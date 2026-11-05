"""Shared helpers for nf-wgs-genotyping bin/ scripts (importable because Python puts a script's
own directory on sys.path). Python >= 3.9; depends on pandas and pyyaml only."""
import csv
import gzip
import json
import re

import yaml

SAMPLESHEET_REQUIRED = ("tumor_id", "normal_id", "caller", "vcf")

# VCF breakend ALT: t[p[  t]p]  ]p]t  [p[t   (t may carry inserted sequence)
_BND_RE = re.compile(r"^([A-Za-z]*)([\[\]])([^:\[\]]+):(\d+)([\[\]])([A-Za-z]*)$")


def open_text(path, mode="rt"):
    """Open a plain or gzip/bgzip text file."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode.replace("t", "") or "r")


def iter_pass_alleles(vcf, tally, contigs=None, missing_filter_is_pass=False):
    """Yield (chrom, pos, REF, ALT) (upper case) for every ALT of every PASS record on an allowed
    contig of a plain or gzipped VCF. FILTER must be exactly PASS ('.' too with
    missing_filter_is_pass); multiallelic records give one item per ALT. Counts go to the
    Counter ``tally``: records, not_pass, off_contig."""
    allowed = set(contigs) if contigs else None
    with open_text(vcf) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            tally["records"] += 1
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alts, filt = f[0], int(f[1]), f[3].upper(), f[4], f[6]
            if not (filt == "PASS" or (filt == "." and missing_filter_is_pass)):
                tally["not_pass"] += 1
                continue
            if allowed is not None and chrom not in allowed:
                tally["off_contig"] += 1
                continue
            for alt in alts.split(","):
                yield chrom, pos, ref, alt.upper()


def parse_contigs(value):
    """'1, 2,X' -> ['1', '2', 'X']; empty/None -> None (meaning: no restriction)."""
    if value is None:
        return None
    items = [c.strip() for c in str(value).split(",") if c.strip()]
    return items or None


def contig_order(fai=None, contigs=None, bam=None):
    """Contig -> rank. Precedence: .fai order, then BAM header order, then the given list."""
    names = []
    if fai:
        with open(fai) as fh:
            names = [line.split("\t", 1)[0] for line in fh if line.strip()]
    elif bam:
        import pysam

        with pysam.AlignmentFile(bam) as fh:
            names = list(fh.references)
    elif contigs:
        names = list(contigs)
    return {name: i for i, name in enumerate(names)}


def site_sort_key(order):
    """Sort key over (chrom, pos): known contigs by rank, unknown ones after, by name."""
    n = len(order)

    def key(chrom, pos):
        rank = order.get(chrom)
        return (rank if rank is not None else n, "" if rank is not None else chrom, int(pos))

    return key


def parse_bnd_alt(alt):
    """Parse a VCF breakend ALT. Returns dict(strand1, chrom2, pos2, strand2, ins_seq) or None.

    strand '+' = the retained sequence lies left of (ends at) the breakend position;
    strand '-' = it lies right of (starts at) the position. Positions are the retained base
    adjacent to the junction (VCF 4.2 spec)."""
    m = _BND_RE.match(alt)
    if not m:
        return None
    pre, b1, chrom2, pos2, b2, post = m.groups()
    if b1 != b2 or bool(pre) == bool(post):
        return None
    if pre:
        strand1, ins_seq = "+", pre[1:]
    else:
        strand1, ins_seq = "-", post[:-1]
    strand2 = "-" if b1 == "[" else "+"
    return {"strand1": strand1, "chrom2": chrom2, "pos2": int(pos2), "strand2": strand2,
            "ins_seq": ins_seq}


def is_single_breakend(alt):
    """GRIDSS-style single breakend ('.SEQ' or 'SEQ.')."""
    return len(alt) > 1 and (alt.startswith(".") or alt.endswith(".")) and "[" not in alt \
        and "]" not in alt


def format_bnd_alt(ref_base, strand1, chrom2, pos2, strand2):
    """Inverse of parse_bnd_alt (no inserted sequence)."""
    b = "[" if strand2 == "-" else "]"
    mate = f"{b}{chrom2}:{int(pos2)}{b}"
    return f"{ref_base}{mate}" if strand1 == "+" else f"{mate}{ref_base}"


def write_table(df, path, dtypes):
    """Write a gzipped CSV with header plus a csverve-style '<path>.yaml' schema."""
    missing = [c for c in df.columns if c not in dtypes]
    if missing:
        raise ValueError(f"no dtype declared for columns: {missing}")
    df.to_csv(path, index=False, sep=",", compression="gzip", na_rep="NaN")
    meta = {"header": True, "sep": ",",
            "columns": [{"name": c, "dtype": dtypes[c]} for c in df.columns]}
    with open(f"{path}.yaml", "w") as fh:
        yaml.safe_dump(meta, fh, sort_keys=False)


def read_samplesheet(path):
    """Read the pipeline samplesheet; require tumor_id, normal_id, caller, vcf columns."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in SAMPLESHEET_REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"samplesheet {path} is missing columns: {missing}")
        return [dict(row) for row in reader]


def write_manifest(path, chunks):
    """chunks: iterable of (chunk_id, n_items)."""
    with open(path, "w") as fh:
        fh.write("chunk_id\tn_items\n")
        for chunk_id, n in chunks:
            fh.write(f"{chunk_id}\t{int(n)}\n")


def read_manifest(path):
    with open(path) as fh:
        next(fh)
        return {cid: int(n) for cid, n in (line.rstrip("\n").split("\t") for line in fh if line.strip())}


def check_gather(manifest_path, received_ids):
    """Fail unless exactly the manifest's chunks were gathered (explicit scatter/gather check)."""
    expected = set(read_manifest(manifest_path))
    received = list(received_ids)
    dupes = sorted({c for c in received if received.count(c) > 1})
    missing = sorted(expected - set(received))
    unexpected = sorted(set(received) - expected)
    if missing or unexpected or dupes:
        raise ValueError(f"gather mismatch: missing={missing} unexpected={unexpected} "
                         f"duplicated={dupes}")


def write_json(obj, path):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True, default=str)
        fh.write("\n")
