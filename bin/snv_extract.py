#!/usr/bin/env python3
"""Extract PASS, biallelic single-nucleotide variants from one caller VCF.

Rules (identical for every caller; the caller name is provenance only):
  * FILTER must be exactly PASS (a missing FILTER '.' counts only with --missing-filter-is-pass);
  * multiallelic records are split into one row per ALT;
  * an allele is kept only if REF and ALT are single bases in ACGT and differ
    (indels, MNVs, symbolic and breakend alleles are dropped);
  * records on contigs outside --contigs are dropped;
  * with --reference, REF must match the FASTA; mismatches are dropped and the run fails if they
    exceed --max-ref-mismatch-frac of SNV alleles (a sign the VCF is on another assembly).
"""
import argparse
import csv
import sys
from collections import Counter

import geno_utils as gu

BASES = set("ACGT")
COLUMNS = ["chrom", "pos", "ref", "alt", "tumor_id", "normal_id", "caller"]


def extract(vcf, caller, tumor_id, normal_id, contigs=None, fasta=None, missing_filter_is_pass=False,
            max_ref_mismatch_frac=0.01):
    fa = None
    if fasta:
        import pysam

        fa = pysam.FastaFile(fasta)
    dropped = Counter()
    rows, seen = [], set()
    for chrom, pos, ref, alt in gu.iter_pass_alleles(vcf, dropped, contigs, missing_filter_is_pass):
        if len(ref) != 1 or len(alt) != 1 or ref not in BASES or alt not in BASES or ref == alt:
            dropped["not_snv"] += 1
            continue
        if fa is not None:
            try:
                fa_base = fa.fetch(chrom, pos - 1, pos).upper()
            except (KeyError, ValueError):
                fa_base = None
            if fa_base != ref:
                dropped["ref_mismatch"] += 1
                continue
        key = (chrom, pos, ref, alt)
        if key in seen:
            dropped["duplicate"] += 1
            continue
        seen.add(key)
        rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "tumor_id": tumor_id,
                     "normal_id": normal_id, "caller": caller})
    records = dropped["records"]
    snv_alleles = len(rows) + dropped["duplicate"] + dropped["ref_mismatch"]
    if fa is not None and snv_alleles and dropped["ref_mismatch"] / snv_alleles > max_ref_mismatch_frac:
        raise ValueError(
            f"{vcf}: REF mismatch against {fasta} for {dropped['ref_mismatch']}/{snv_alleles} SNV alleles "
            f"(> {max_ref_mismatch_frac}); wrong assembly?")
    stats = {"input": str(vcf), "caller": caller, "tumor_id": tumor_id, "normal_id": normal_id,
             "records": records, "kept": len(rows),
             "dropped": {k: dropped.get(k, 0) for k in
                         ("not_pass", "not_snv", "off_contig", "ref_mismatch", "duplicate")}}
    return rows, stats


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--caller", required=True)
    ap.add_argument("--tumor-id", required=True)
    ap.add_argument("--normal-id", required=True)
    ap.add_argument("--contigs", default=None, help="comma-separated contigs to keep (default: all)")
    ap.add_argument("--reference", default=None, help="indexed FASTA for the REF check")
    ap.add_argument("--missing-filter-is-pass", action="store_true")
    ap.add_argument("--max-ref-mismatch-frac", type=float, default=0.01)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-stats", required=True)
    a = ap.parse_args(argv)
    rows, stats = extract(a.vcf, a.caller, a.tumor_id, a.normal_id, gu.parse_contigs(a.contigs), a.reference,
                          a.missing_filter_is_pass, a.max_ref_mismatch_frac)
    with open(a.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    gu.write_json(stats, a.out_stats)
    return 0


if __name__ == "__main__":
    sys.exit(main())
