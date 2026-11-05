#!/usr/bin/env python3
"""Genotype SV insertions (svtyper skips SVTYPE=INS) by realigning reads to the reference and
insertion haplotypes, with the indel machinery of indel_count.py (read filters, soft clips
included, low-BQ bases masked, semi-global BWA-MEM scoring, each fragment counted once,
assignment when best - second best >= --min-score-diff).

Each insertion adjacency is genotyped on its own (as svtyper does per record): haplotypes are the
reference window and the window with the del_len bases after the anchor replaced by the inserted
block. The block is the inserted sequence when known; for length-only calls it is a run of N,
which the aligner scores 0 against any base, so only reads crossing the insertion point decide
(a reference read must carry >= --min-score-diff reference bases past it; an insertion read wins
by every flank base it does not have to mismatch). Blocks longer than 2 x (longest read + 10)
keep their first and last (longest read + 10) bases: no read can span more, and the haplotypes
stay short.

ins_status: genotyped (ref + alt > 0), no_evidence (no informative fragment), skipped_imprecise
(IMPRECISE call: its position is not known well enough to build a haplotype), skipped_no_reference
(no --reference). Skipped rows have empty counts.
"""
import argparse
import sys

import pandas as pd
import pysam

import indel_count as ic
import indel_utils as iu
from sv_union import INS_SITE_COLUMNS as SITE_COLUMNS   # noqa: F401  (input format, written by sv_union)

COLUMNS = ["adj_id", "ref_fragments", "alt_fragments", "uninformative_fragments", "conflict_fragments",
           "ins_status"]
_COUNTS = COLUMNS[1:5]


def _flag(v):
    return v if isinstance(v, bool) else str(v).strip().lower() in ("true", "1")


def alt_block(block, max_read_len):
    """Inserted block bounded to its first and last (max_read_len + FLANK_PAD) bases."""
    keep = max_read_len + ic.FLANK_PAD
    return block if len(block) <= 2 * keep else block[:keep] + block[-keep:]


def _insertion_counts(bam, fa, aligner, chrom, anchor, del_len, block, min_mapq, min_baseq, count_duplicates,
                      min_score_diff):
    """(ref, alt, uninformative, conflict) fragments of one insertion."""
    if chrom not in fa.references:
        return 0, 0, 0, 0
    n = fa.get_reference_length(chrom)
    fetch = lambda a, b: fa.fetch(chrom, max(a, 0), min(b, n)).upper()   # noqa: E731
    ref = fetch(anchor - 1, anchor + del_len)
    vs, ve = iu.variable_span(fetch, n, anchor, ref, ref[0] + block)
    fs, fe = iu.footprint((vs, ve), n)
    reads = ic.locus_reads(bam, chrom, fs, fe, n, min_mapq, count_duplicates)
    if not reads:
        return 0, 0, 0, 0
    max_len = max(a.query_length for a in reads)
    allele = (anchor, ref, ref[0] + alt_block(block, max_len))
    haps = ic.haplotypes(fetch, n, vs, ve, [allele], max_len)
    (r, a), uninformative, conflict = ic.assign_fragments(aligner, reads, haps, min_baseq, min_score_diff)
    return r, a, uninformative, conflict


def count_insertions(bam_path, fasta_path, sites, min_mapq=20, min_baseq=20, count_duplicates=False,
                     min_score_diff=5):
    """sites: iterable of mappings with SITE_COLUMNS. Returns one dict per site (COLUMNS), in order."""
    out = []
    aligner = ic.make_aligner()
    fa = pysam.FastaFile(fasta_path) if fasta_path else None
    try:
        with pysam.AlignmentFile(bam_path) as bam:
            for s in sites:
                row = {"adj_id": str(s["adj_id"]), **{k: None for k in _COUNTS}}
                if _flag(s["imprecise"]):
                    row["ins_status"] = "skipped_imprecise"
                elif fa is None:
                    row["ins_status"] = "skipped_no_reference"
                else:
                    seq = str(s["ins_seq"]) if _flag(s["ins_seq_available"]) else "N" * int(s["ins_len"])
                    counts = _insertion_counts(bam, fa, aligner, str(s["chromosome"]), int(s["position"]),
                                               int(s["del_len"]), seq.upper(), min_mapq, min_baseq,
                                               count_duplicates, min_score_diff)
                    row.update(zip(_COUNTS, counts))
                    row["ins_status"] = "genotyped" if counts[0] + counts[1] > 0 else "no_evidence"
                out.append(row)
    finally:
        if fa is not None:
            fa.close()
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--reference", default=None, help="indexed FASTA (without it nothing is genotyped)")
    ap.add_argument("--sites", required=True, help="insertion chunk TSV from sv_union.py")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--count-duplicates", action="store_true")
    ap.add_argument("--min-score-diff", type=float, default=5)
    a = ap.parse_args(argv)
    if a.min_score_diff <= 0:
        ap.error("--min-score-diff must be > 0")
    sites = pd.read_csv(a.sites, sep="\t", dtype={"adj_id": str, "chromosome": str, "ins_seq": str},
                        keep_default_na=False)
    rows = count_insertions(a.bam, a.reference or None, sites.to_dict("records"), a.min_mapq, a.min_baseq,
                            a.count_duplicates, a.min_score_diff)
    df = pd.DataFrame(rows, columns=COLUMNS)
    for c in _COUNTS:
        df[c] = pd.array(df[c].tolist(), dtype="Int64") if len(df) else pd.array([], dtype="Int64")
    df.to_csv(a.out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
