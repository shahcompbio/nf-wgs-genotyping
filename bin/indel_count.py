#!/usr/bin/env python3
"""Count reference / alternate / other fragments at indel loci in a bulk BAM by realigning each
read to the reference and alternate haplotypes of its locus.

For each locus (alleles A1..Ak sharing a footprint; see indel_union.py):
  * haplotypes: H0 = reference window around the locus; Hi = the same window carrying Ai. The
    window extends (longest read + longest allele + 10) bases beyond every equivalent placement
    of the events, so a read can align anywhere it could have come from;
  * reads: every read whose *unclipped* span (soft clips included) reaches the locus footprint,
    so reads whose indel the aligner soft-clipped are examined. Dropped: unmapped, secondary,
    supplementary, QC-fail, duplicate (unless --count-duplicates), MAPQ < --min-mapq;
  * each read's full sequence (soft clips included; bases with BQ < --min-baseq replaced by N,
    which scores 0 against any base) is aligned to every haplotype semi-globally: the whole read
    must align, haplotype ends are free. Scores are BWA-MEM's defaults: match +1, mismatch -4,
    gap of length k -(6 + k);
  * each fragment (query name) is counted once: its reads' score vectors are summed, so
    overlapping mates count once and a mate that cannot discriminate adds nothing. The fragment
    supports the best haplotype when it beats the second best by >= --min-score-diff (default 5
    = one high-quality discriminating base: match - mismatch). Otherwise it is a conflict when
    its reads, taken alone, support different haplotypes, and uninformative when they do not
    (e.g. a read that ends inside the repeat).
Per allele Ai: ref_count (H0), alt_count (Hi), other_count (another allele of the locus);
uninformative_count and conflict_count are per locus.
"""
import argparse
import sys

import pandas as pd
import pysam
from Bio.Align import PairwiseAligner

import indel_utils as iu

KEYS = ["chromosome", "position", "ref", "alt"]
COUNT_COLUMNS = KEYS + ["ref_count", "alt_count", "other_count", "uninformative_count", "conflict_count"]
FETCH_MARGIN = 1000          # bp around the footprint searched for reads (aligned coordinates)
FLANK_PAD = 10


def make_aligner():
    """Semi-global (read end-to-end, haplotype ends free), BWA-MEM scoring; Biopython's gap open
    score is the cost of the first gap position (O + E)."""
    aligner = PairwiseAligner(mode="global", match_score=1, mismatch_score=-4, open_gap_score=-7,
                              extend_gap_score=-1, wildcard="N")
    aligner.query_end_gap_score = 0       # score(haplotype, read): leading/trailing haplotype bases
    return aligner


def _clips(cigar):
    def lead(ops):
        n = 0
        for op, length in ops:
            if op == 5:
                continue
            if op != 4:
                break
            n += length
        return n
    return lead(cigar), lead(reversed(cigar))


def _keep(a, min_mapq, count_duplicates):
    return not (a.is_unmapped or a.is_secondary or a.is_supplementary or a.is_qcfail
                or (a.is_duplicate and not count_duplicates) or a.mapping_quality < min_mapq
                or a.query_sequence is None or a.reference_end is None)


def _masked(a, min_baseq):
    seq, quals = a.query_sequence.upper(), a.query_qualities
    if quals is None or min_baseq <= 0:
        return seq
    return "".join(b if q >= min_baseq else "N" for b, q in zip(seq, quals))


def _best(scores):
    order = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    return order[0], scores[order[0]] - scores[order[1]]


def locus_reads(bam, chrom, fs, fe, contig_len, min_mapq, count_duplicates):
    """Reads passing the filters whose *unclipped* span reaches [fs, fe) (0-based)."""
    try:
        candidates = bam.fetch(chrom, max(0, fs - FETCH_MARGIN), min(contig_len, fe + FETCH_MARGIN))
    except ValueError:                    # contig absent from the BAM header
        return []
    reads = []
    for a in candidates:
        if not _keep(a, min_mapq, count_duplicates):
            continue
        lead, trail = _clips(a.cigartuples)
        if a.reference_start - lead < fe and a.reference_end + trail > fs:
            reads.append(a)
    return reads


def haplotypes(fetch, contig_len, vs, ve, alleles, max_read_len):
    """[H0 = reference window, H1..Hk = window carrying each allele]; the window extends
    (longest read + longest allele + FLANK_PAD) beyond the variable span [vs, ve)."""
    flank = max_read_len + max(max(len(r), len(t)) for _, r, t in alleles) + FLANK_PAD
    ws, we = max(0, vs - flank), min(contig_len, ve + flank)
    window = fetch(ws, we)
    return [window] + [iu.apply_allele(window, ws, pos, ref, alt) for pos, ref, alt in alleles]


def assign_fragments(aligner, reads, haps, min_baseq, min_score_diff):
    """Per-haplotype fragment counts, uninformative and conflict fragments (module docstring)."""
    fragments = {}
    for a in reads:
        query = _masked(a, min_baseq)
        fragments.setdefault(a.query_name, []).append([aligner.score(h, query) for h in haps])
    counts, uninformative, conflict = [0] * len(haps), 0, 0
    for per_read in fragments.values():
        total = [sum(col) for col in zip(*per_read)]
        best, margin = _best(total)
        if margin >= min_score_diff:
            counts[best] += 1
            continue
        calls = {b for b, m in map(_best, per_read) if m >= min_score_diff}
        if len(calls) > 1:
            conflict += 1
        else:
            uninformative += 1
    return counts, uninformative, conflict


def _locus_counts(bam, fa, aligner, chrom, alleles, min_mapq, min_baseq, count_duplicates, min_score_diff):
    """alleles: [(pos, ref, alt)] -> (per-haplotype fragment counts [H0..Hk], uninformative, conflict)."""
    zero = ([0] * (len(alleles) + 1), 0, 0)
    contig_len = fa.get_reference_length(chrom)
    fetch = lambda a, b: fa.fetch(chrom, max(a, 0), min(b, contig_len)).upper()   # noqa: E731
    spans = [iu.variable_span(fetch, contig_len, pos, ref, alt) for pos, ref, alt in alleles]
    vs, ve = min(s[0] for s in spans), max(s[1] for s in spans)
    fs, fe = iu.footprint((vs, ve), contig_len)
    reads = locus_reads(bam, chrom, fs, fe, contig_len, min_mapq, count_duplicates)
    if not reads:
        return zero
    haps = haplotypes(fetch, contig_len, vs, ve, alleles, max(a.query_length for a in reads))
    return assign_fragments(aligner, reads, haps, min_baseq, min_score_diff)


def count_sites(bam_path, fasta_path, sites, min_mapq=20, min_baseq=20, count_duplicates=False,
                min_score_diff=5):
    """sites: iterable of (chromosome, position, ref, alt, locus_id); alleles sharing a locus_id
    are genotyped jointly. Returns one dict per site (COUNT_COLUMNS), in input order."""
    sites = [(str(c), int(p), r.upper(), t.upper(), str(l)) for c, p, r, t, l in sites]
    loci = {}
    for i, site in enumerate(sites):
        loci.setdefault(site[4], []).append(i)
    out = [None] * len(sites)
    aligner = make_aligner()
    with pysam.AlignmentFile(bam_path) as bam, pysam.FastaFile(fasta_path) as fa:
        for idx in loci.values():
            chrom = sites[idx[0]][0]
            if any(sites[i][0] != chrom for i in idx):
                raise ValueError(f"locus {sites[idx[0]][4]} spans several contigs")
            alleles = [sites[i][1:4] for i in idx]
            counts, uninformative, conflict = _locus_counts(bam, fa, aligner, chrom, alleles, min_mapq,
                                                            min_baseq, count_duplicates, min_score_diff)
            for j, i in enumerate(idx):
                c, p, r, t, _ = sites[i]
                alt_n = counts[j + 1]
                out[i] = {"chromosome": c, "position": p, "ref": r, "alt": t, "ref_count": counts[0],
                          "alt_count": alt_n, "other_count": sum(counts[1:]) - alt_n,
                          "uninformative_count": uninformative, "conflict_count": conflict}
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--reference", required=True, help="indexed FASTA")
    ap.add_argument("--sites", required=True, help="TSV with chromosome, position, ref, alt, locus_id")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--count-duplicates", action="store_true")
    ap.add_argument("--min-score-diff", type=float, default=5)
    a = ap.parse_args(argv)
    if a.min_score_diff <= 0:
        ap.error("--min-score-diff must be > 0")
    sites = pd.read_csv(a.sites, sep="\t", dtype={"chromosome": str, "ref": str, "alt": str, "locus_id": str})
    rows = count_sites(a.bam, a.reference, sites[KEYS + ["locus_id"]].itertuples(index=False), a.min_mapq,
                       a.min_baseq, a.count_duplicates, a.min_score_diff)
    pd.DataFrame(rows, columns=COUNT_COLUMNS).to_csv(a.out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
