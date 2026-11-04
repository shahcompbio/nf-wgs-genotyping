#!/usr/bin/env python3
"""Normalize one caller's SV VCF (BRASS, GRIDSS, SvABA or generic) to canonical BEDPE junctions.

Pass rules
  brass   : FILTER PASS, or '.' (BRASS writes no FILTER) unless --no-brass-missing-filter-is-pass
  svaba   : FILTER PASS
  gridss  : FILTER PASS and the matched normal's FORMAT/VF <= --gridss-max-normal-vf
            (GRIDSS 2.2.2 'somatic.sv.vcf' is a joint normal+tumor call set)
  generic : FILTER PASS or '.'

Breakends are read from the ALT field (VCF 4.2: position = retained base adjacent to the
junction; strand '+' = sequence kept to the left, '-' = kept to the right), so a mate missing from
the file is still recovered. Symbolic <DEL>, <DUP>(:TANDEM) and <INV> are expanded into
junctions. Dropped and counted: single breakends, other alleles (e.g. <INS>), junctions touching
contigs outside --contigs, and intrachromosomal junctions shorter than --min-size. Ends are put in
canonical order (contig rank, position, strand) and mate records collapse to one junction.
"""
import argparse
import csv
import sys
from collections import Counter

import geno_utils as gu

COLUMNS = ["chrom1", "pos1", "strand1", "chrom2", "pos2", "strand2", "ci1_lo", "ci1_hi", "ci2_lo", "ci2_hi",
           "imprecise", "sv_class", "caller", "tumor_id", "normal_id", "record_id", "filter", "qual", "ins_len"]
CALLERS = ("brass", "gridss", "svaba", "generic")
DROP_REASONS = ("not_pass", "germline_normal_support", "single_breakend", "unsupported_allele", "off_contig",
                "small", "mate_duplicate")


def sv_class(chrom1, strand1, chrom2, strand2):
    if chrom1 != chrom2:
        return "TRA"
    return {("+", "-"): "DEL", ("-", "+"): "DUP"}.get((strand1, strand2), "INV")


def _info(field):
    out = {}
    if field in (".", ""):
        return out
    for item in field.split(";"):
        k, _, v = item.partition("=")
        out[k] = v if _ else True
    return out


def _ci(info, key, imprecise, default_ci):
    if key in info and info[key] is not True:
        lo, hi = (int(x) for x in str(info[key]).split(",")[:2])
        return lo, hi
    return (-default_ci, default_ci) if imprecise else (0, 0)


def _end_key(order, chrom, pos, strand):
    rank = order.get(chrom)
    return (rank if rank is not None else len(order), "" if rank is not None else chrom, pos, strand)


def _junctions(chrom, pos, alt, info, imprecise, imprecise_ci):
    """Return list of (end1, end2, ins_len) with end = (chrom, pos, strand, (ci_lo, ci_hi)), or a drop reason."""
    if gu.is_single_breakend(alt):
        return "single_breakend"
    bnd = gu.parse_bnd_alt(alt)
    if bnd:
        own = _ci(info, "CIPOS", imprecise, imprecise_ci)
        remote = _ci(info, "CIRPOS", imprecise, imprecise_ci) if "CIRPOS" in info else own
        return [((chrom, pos, bnd["strand1"], own), (bnd["chrom2"], bnd["pos2"], bnd["strand2"], remote),
                 len(bnd["ins_seq"]))]
    if alt.startswith("<") and "END" in info:
        svtype = alt.strip("<>").split(":")[0]
        end = int(info["END"])
        cip, cie = _ci(info, "CIPOS", imprecise, imprecise_ci), _ci(info, "CIEND", imprecise, imprecise_ci)
        if svtype == "DEL":
            return [((chrom, pos, "+", cip), (chrom, end + 1, "-", cie), 0)]
        if svtype == "DUP":
            return [((chrom, pos + 1, "-", cip), (chrom, end, "+", cie), 0)]
        if svtype == "INV":
            return [((chrom, pos, "+", cip), (chrom, end, "+", cie), 0),
                    ((chrom, pos + 1, "-", cip), (chrom, end + 1, "-", cie), 0)]
    return "unsupported_allele"


def normalize(vcf, caller, tumor_id, normal_id, contigs=None, order=None, min_size=50, imprecise_ci=150,
              gridss_max_normal_vf=0, brass_missing_filter_is_pass=True):
    if caller not in CALLERS:
        raise ValueError(f"unknown SV caller {caller!r}; expected one of {CALLERS}")
    order = order or {}
    allowed = set(contigs) if contigs else None
    dropped, rows, seen = Counter(), [], set()
    records, normal_col = 0, None
    with gu.open_text(vcf) as fh:
        for line in fh:
            if line.startswith("##") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                if caller == "gridss":
                    if normal_id not in samples:
                        raise ValueError(f"{vcf}: GRIDSS normal sample column {normal_id!r} not found in {samples}")
                    normal_col = 9 + samples.index(normal_id)
                continue
            records += 1
            chrom, pos, rid, alts, qual, filt = f[0], int(f[1]), f[2], f[4], f[5], f[6]
            info = _info(f[7])
            if caller == "brass":
                ok = filt == "PASS" or (filt == "." and brass_missing_filter_is_pass)
            elif caller == "generic":
                ok = filt in ("PASS", ".")
            else:
                ok = filt == "PASS"
            if not ok:
                dropped["not_pass"] += 1
                continue
            if caller == "gridss":
                fmt = dict(zip(f[8].split(":"), f[normal_col].split(":")))
                vf = fmt.get("VF", "0")
                if float(0 if vf in (".", "") else vf) > gridss_max_normal_vf:
                    dropped["germline_normal_support"] += 1
                    continue
            imprecise = "IMPRECISE" in info
            for alt in alts.split(","):
                junctions = _junctions(chrom, pos, alt, info, imprecise, imprecise_ci)
                if isinstance(junctions, str):
                    dropped[junctions] += 1
                    continue
                for e1, e2, ins_len in junctions:
                    if allowed is not None and (e1[0] not in allowed or e2[0] not in allowed):
                        dropped["off_contig"] += 1
                        continue
                    if e1[0] == e2[0] and abs(e2[1] - e1[1]) < min_size:
                        dropped["small"] += 1
                        continue
                    if _end_key(order, *e2[:3]) < _end_key(order, *e1[:3]):
                        e1, e2 = e2, e1
                    key = (*e1[:3], *e2[:3])
                    if key in seen:
                        dropped["mate_duplicate"] += 1
                        continue
                    seen.add(key)
                    rows.append({
                        "chrom1": e1[0], "pos1": e1[1], "strand1": e1[2],
                        "chrom2": e2[0], "pos2": e2[1], "strand2": e2[2],
                        "ci1_lo": e1[3][0], "ci1_hi": e1[3][1], "ci2_lo": e2[3][0], "ci2_hi": e2[3][1],
                        "imprecise": imprecise, "sv_class": sv_class(e1[0], e1[2], e2[0], e2[2]),
                        "caller": caller, "tumor_id": tumor_id, "normal_id": normal_id, "record_id": rid,
                        "filter": filt, "qual": qual, "ins_len": ins_len})
    stats = {"input": str(vcf), "caller": caller, "tumor_id": tumor_id, "normal_id": normal_id,
             "records": records, "kept": len(rows), "dropped": {k: dropped.get(k, 0) for k in DROP_REASONS}}
    return rows, stats


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--caller", required=True, choices=CALLERS)
    ap.add_argument("--tumor-id", required=True)
    ap.add_argument("--normal-id", required=True)
    ap.add_argument("--contigs", default=None)
    ap.add_argument("--reference", default=None, help="FASTA whose .fai sets canonical contig order")
    ap.add_argument("--min-size", type=int, default=50)
    ap.add_argument("--imprecise-ci", type=int, default=150)
    ap.add_argument("--gridss-max-normal-vf", type=float, default=0)
    ap.add_argument("--no-brass-missing-filter-is-pass", dest="brass_missing_filter_is_pass",
                    action="store_false")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-stats", required=True)
    a = ap.parse_args(argv)
    contigs = gu.parse_contigs(a.contigs)
    order = gu.contig_order(fai=f"{a.reference}.fai" if a.reference else None, contigs=contigs)
    rows, stats = normalize(a.vcf, a.caller, a.tumor_id, a.normal_id, contigs, order, a.min_size,
                            a.imprecise_ci, a.gridss_max_normal_vf, a.brass_missing_filter_is_pass)
    with open(a.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    gu.write_json(stats, a.out_stats)
    return 0


if __name__ == "__main__":
    sys.exit(main())
