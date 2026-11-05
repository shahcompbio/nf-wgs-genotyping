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
junctions. Ends are put in canonical order (contig rank, position, strand) and mate records
collapse to one junction.

Insertions (sv_class INS) are written like other junctions, (c, pos, +) - (c, pos + del_len + 1, -)
with pos = the reference base left of the insertion, plus ins_len, ins_seq and ins_seq_available
(ins_seq holds the whole inserted sequence; empty and False for length-only calls). Sources:
symbolic <INS> (and subtypes such as <INS:ME:ALU>; length from SVLEN, SVINSLEN or INSLEN, sequence
from SVINSSEQ), sequence-resolved REF/ALT alleles, and breakends whose ALT carries inserted
sequence. A same-contig (+,-) junction with del_len = pos2 - pos1 - 1 deleted bases is a deletion
when del_len >= --min-size, an insertion when del_len < --min-size and ins_len - del_len >=
--min-size, and too small otherwise. With --reference, insertions with a known sequence are
left-aligned (indel_utils.normalize), so callers placing the same insertion at different points
of a repeat agree. IMPRECISE insertions are kept (CIPOS applies to both ends unless CIEND is set).

Size boundary: an event is kept when its net length change is >= --min-size: deleted bases for DEL
(pos2 - pos1 - 1), duplicated bases for DUP (pos2 - pos1 + 1), net inserted bases for INS; INV
uses pos2 - pos1. The indel workflow keeps |len(ALT) - len(REF)| < indel_max_length, so with
equal thresholds (both 50) every insertion/deletion belongs to exactly one workflow.

Dropped and counted: single breakends, other alleles (e.g. <CNV>, explicit-sequence deletions),
<INS> without length or sequence, junctions touching contigs outside --contigs, events below
--min-size, and mate records of an event already kept.
"""
import argparse
import csv
import re
import sys
from collections import Counter

import geno_utils as gu
import indel_utils as iu

COLUMNS = ["chrom1", "pos1", "strand1", "chrom2", "pos2", "strand2", "ci1_lo", "ci1_hi", "ci2_lo", "ci2_hi",
           "imprecise", "sv_class", "caller", "tumor_id", "normal_id", "record_id", "filter", "qual", "ins_len",
           "ins_seq", "ins_seq_available"]
CALLERS = ("brass", "gridss", "svaba", "generic")
DROP_REASONS = ("not_pass", "germline_normal_support", "single_breakend", "unsupported_allele",
                "insertion_length_unknown", "off_contig", "small", "mate_duplicate")
INS_LENGTH_KEYS = ("SVLEN", "SVINSLEN", "INSLEN")
_PLAIN = re.compile(r"^[ACGTN]+$")


def sv_class(chrom1, strand1, chrom2, strand2):
    if chrom1 != chrom2:
        return "TRA"
    return {("+", "-"): "DEL", ("-", "+"): "DUP"}.get((strand1, strand2), "INV")


def event_size(chrom1, pos1, strand1, chrom2, pos2, strand2):
    """Net reference length change of an intrachromosomal junction (None between contigs):
    deleted bases for (+,-), duplicated bases for (-,+), pos2 - pos1 for inversions."""
    if chrom1 != chrom2:
        return None
    if (strand1, strand2) == ("+", "-"):
        return pos2 - pos1 - 1
    if (strand1, strand2) == ("-", "+"):
        return pos2 - pos1 + 1
    return pos2 - pos1


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


def _ins_length(info):
    for key in INS_LENGTH_KEYS:
        v = info.get(key)
        if v not in (None, True, "", "."):
            return abs(int(float(str(v).split(",")[0])))
    return None


def _junctions(chrom, pos, ref, alt, info, imprecise, imprecise_ci):
    """Return list of (end1, end2, ins_len, ins_seq) with end = (chrom, pos, strand, (ci_lo, ci_hi))
    and ins_seq the inserted sequence (None when only its length is known), or a drop reason."""
    if gu.is_single_breakend(alt):
        return "single_breakend"
    cip = _ci(info, "CIPOS", imprecise, imprecise_ci)
    bnd = gu.parse_bnd_alt(alt)
    if bnd:
        remote = _ci(info, "CIRPOS", imprecise, imprecise_ci) if "CIRPOS" in info else cip
        seq = bnd["ins_seq"].upper()
        return [((chrom, pos, bnd["strand1"], cip), (bnd["chrom2"], bnd["pos2"], bnd["strand2"], remote),
                 len(seq), seq)]
    if alt.startswith("<"):
        svtype = alt.strip("<>").split(":")[0].upper()
        if svtype == "INS":
            seq = info.get("SVINSSEQ")
            seq = seq.upper() if isinstance(seq, str) and seq not in ("", ".") else None
            length = len(seq) if seq is not None else _ins_length(info)
            if length is None:
                return "insertion_length_unknown"
            cie = _ci(info, "CIEND", imprecise, imprecise_ci) if "CIEND" in info else cip
            return [((chrom, pos, "+", cip), (chrom, pos + 1, "-", cie), length, seq)]
        if "END" in info:
            end = int(info["END"])
            cie = _ci(info, "CIEND", imprecise, imprecise_ci)
            if svtype == "DEL":
                return [((chrom, pos, "+", cip), (chrom, end + 1, "-", cie), 0, "")]
            if svtype == "DUP":
                return [((chrom, pos + 1, "-", cip), (chrom, end, "+", cie), 0, "")]
            if svtype == "INV":
                return [((chrom, pos, "+", cip), (chrom, end, "+", cie), 0, ""),
                        ((chrom, pos + 1, "-", cip), (chrom, end + 1, "-", cie), 0, "")]
        return "unsupported_allele"
    ref, alt = ref.upper(), alt.upper()
    if _PLAIN.match(ref) and _PLAIN.match(alt) and len(alt) > len(ref):
        # sequence-resolved insertion (or del-ins): trim the shared prefix, then the shared suffix
        p = 0
        while p < len(ref) and ref[p] == alt[p]:
            p += 1
        r, a = ref[p:], alt[p:]
        while r and a and r[-1] == a[-1]:
            r, a = r[:-1], a[:-1]
        anchor = pos + p - 1
        return [((chrom, anchor, "+", cip), (chrom, anchor + len(r) + 1, "-", cip), len(a), a)]
    return "unsupported_allele"


def _left_align(fa, chrom, anchor, del_len, ins_seq):
    """Leftmost equivalent placement of a sequence-resolved insertion (deleting del_len bases after
    `anchor`): returns (anchor, del_len, ins_seq). Unchanged when the contig is not in the FASTA."""
    if chrom not in fa.references:
        return anchor, del_len, ins_seq
    n = fa.get_reference_length(chrom)
    fetch = lambda a, b: fa.fetch(chrom, max(a, 0), min(b, n)).upper()   # noqa: E731
    base = fetch(anchor - 1, anchor + del_len)
    npos, nref, nalt = iu.normalize(fetch, anchor, base, base[0] + ins_seq)
    if nref[0] == nalt[0]:
        return npos, len(nref) - 1, nalt[1:]
    return npos - 1, len(nref), nalt


def normalize(vcf, caller, tumor_id, normal_id, contigs=None, order=None, min_size=50, imprecise_ci=150,
              gridss_max_normal_vf=0, brass_missing_filter_is_pass=True, fasta=None):
    if caller not in CALLERS:
        raise ValueError(f"unknown SV caller {caller!r}; expected one of {CALLERS}")
    order = order or {}
    allowed = set(contigs) if contigs else None
    dropped, rows, seen, ins_tally = Counter(), [], set(), Counter()
    records, normal_col = 0, None
    fa = None
    if fasta:
        import pysam

        fa = pysam.FastaFile(fasta)
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
            chrom, pos, rid, ref, alts, qual, filt = f[0], int(f[1]), f[2], f[3], f[4], f[5], f[6]
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
                junctions = _junctions(chrom, pos, ref, alt, info, imprecise, imprecise_ci)
                if isinstance(junctions, str):
                    dropped[junctions] += 1
                    continue
                for e1, e2, ins_len, ins_seq in junctions:
                    if allowed is not None and (e1[0] not in allowed or e2[0] not in allowed):
                        dropped["off_contig"] += 1
                        continue
                    if _end_key(order, *e2[:3]) < _end_key(order, *e1[:3]):
                        e1, e2 = e2, e1
                    cls = sv_class(e1[0], e1[2], e2[0], e2[2])
                    size = event_size(e1[0], e1[1], e1[2], e2[0], e2[1], e2[2])
                    if cls == "DEL" and 0 <= size < min_size and ins_len - size >= min_size:
                        cls = "INS"
                    elif size is not None and size < min_size:
                        dropped["small"] += 1
                        continue
                    key = (*e1[:3], *e2[:3])
                    if cls == "INS":
                        if ins_seq is not None and fa is not None:
                            anchor, del_len, ins_seq = _left_align(fa, e1[0], e1[1], size, ins_seq)
                            if (anchor, del_len) != (e1[1], size):
                                ins_tally["left_aligned"] += 1
                            e1 = (e1[0], anchor, "+", e1[3])
                            e2 = (e2[0], anchor + del_len + 1, "-", e2[3])
                            ins_len = len(ins_seq)
                        key = (*e1[:3], *e2[:3], ins_seq if ins_seq is not None else f"len:{ins_len}")
                    if key in seen:
                        dropped["mate_duplicate"] += 1
                        continue
                    seen.add(key)
                    if cls == "INS":
                        ins_tally["kept"] += 1
                        ins_tally["sequence_resolved" if ins_seq is not None else "length_only"] += 1
                    rows.append({
                        "chrom1": e1[0], "pos1": e1[1], "strand1": e1[2],
                        "chrom2": e2[0], "pos2": e2[1], "strand2": e2[2],
                        "ci1_lo": e1[3][0], "ci1_hi": e1[3][1], "ci2_lo": e2[3][0], "ci2_hi": e2[3][1],
                        "imprecise": imprecise, "sv_class": cls,
                        "caller": caller, "tumor_id": tumor_id, "normal_id": normal_id, "record_id": rid,
                        "filter": filt, "qual": qual, "ins_len": ins_len, "ins_seq": ins_seq or "",
                        "ins_seq_available": ins_seq is not None and len(ins_seq) == ins_len})
    if fa is not None:
        fa.close()
    stats = {"input": str(vcf), "caller": caller, "tumor_id": tumor_id, "normal_id": normal_id,
             "records": records, "kept": len(rows), "dropped": {k: dropped.get(k, 0) for k in DROP_REASONS},
             "insertions": {k: ins_tally.get(k, 0) for k in ("kept", "sequence_resolved", "length_only",
                                                             "left_aligned")}}
    return rows, stats


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--caller", required=True, choices=CALLERS)
    ap.add_argument("--tumor-id", required=True)
    ap.add_argument("--normal-id", required=True)
    ap.add_argument("--contigs", default=None)
    ap.add_argument("--reference", default=None,
                    help="indexed FASTA: its .fai sets canonical contig order; insertions are left-aligned on it")
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
                            a.imprecise_ci, a.gridss_max_normal_vf, a.brass_missing_filter_is_pass, a.reference)
    with open(a.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    gu.write_json(stats, a.out_stats)
    return 0


if __name__ == "__main__":
    sys.exit(main())
