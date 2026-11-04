#!/usr/bin/env python3
"""Count reference / alternate / other fragments at SNV sites in a bulk BAM.

Bulk port of mondrianutils SnvGenotyper with every filter explicit (no pysam defaults):
  * reads dropped: unmapped, secondary, supplementary, QC-fail, duplicate (unless
    --count-duplicates), MAPQ < --min-mapq;
  * pileup positions dropped: deletions, reference skips, base quality < --min-baseq;
  * improper pairs / orphans are kept (pysam ignore_orphans=False);
  * each fragment (query name) is counted once: overlapping mates that agree count once; if they
    disagree the higher-quality base wins; a quality tie is a conflict and is not counted.
"""
import argparse
import sys

import pandas as pd
import pysam

COUNT_COLUMNS = ["chromosome", "position", "ref", "alt", "ref_count", "alt_count", "other_count",
                 "conflict_count"]


def _site_counts(bam, chrom, pos, ref, alt, min_mapq, min_baseq, count_duplicates):
    obs = {}
    try:
        columns = bam.pileup(chrom, pos - 1, pos, truncate=True, stepper="nofilter", ignore_orphans=False,
                             ignore_overlaps=False, min_base_quality=0, min_mapping_quality=0,
                             max_depth=10 ** 6)
        for col in columns:
            if col.reference_pos != pos - 1:
                continue
            for pr in col.pileups:
                a = pr.alignment
                if (a.is_unmapped or a.is_secondary or a.is_supplementary or a.is_qcfail
                        or (a.is_duplicate and not count_duplicates) or a.mapping_quality < min_mapq):
                    continue
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                qp = pr.query_position
                bq = a.query_qualities[qp] if a.query_qualities is not None else 0
                if bq < min_baseq:
                    continue
                obs.setdefault(a.query_name, []).append((a.query_sequence[qp].upper(), bq))
    except ValueError:
        # contig absent from the BAM header: report zero counts
        pass
    ref_n = alt_n = other_n = conflict_n = 0
    for reads in obs.values():
        bases = {b for b, _ in reads}
        if len(bases) == 1:
            base = reads[0][0]
        else:
            best = max(q for _, q in reads)
            top = {b for b, q in reads if q == best}
            if len(top) != 1:
                conflict_n += 1
                continue
            base = top.pop()
        if base == ref:
            ref_n += 1
        elif base == alt:
            alt_n += 1
        else:
            other_n += 1
    return ref_n, alt_n, other_n, conflict_n


def count_sites(bam_path, sites, min_mapq=20, min_baseq=20, count_duplicates=False):
    out = []
    with pysam.AlignmentFile(bam_path) as bam:
        for chrom, pos, ref, alt in sites:
            r, a, o, c = _site_counts(bam, str(chrom), int(pos), ref.upper(), alt.upper(), min_mapq,
                                      min_baseq, count_duplicates)
            out.append({"chromosome": str(chrom), "position": int(pos), "ref": ref, "alt": alt,
                        "ref_count": r, "alt_count": a, "other_count": o, "conflict_count": c})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sites", required=True, help="TSV with chromosome, position, ref, alt")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--count-duplicates", action="store_true")
    a = ap.parse_args(argv)
    sites = pd.read_csv(a.sites, sep="\t", dtype={"chromosome": str, "ref": str, "alt": str})
    rows = count_sites(a.bam, sites[["chromosome", "position", "ref", "alt"]].itertuples(index=False),
                       a.min_mapq, a.min_baseq, a.count_duplicates)
    pd.DataFrame(rows, columns=COUNT_COLUMNS).to_csv(a.out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
