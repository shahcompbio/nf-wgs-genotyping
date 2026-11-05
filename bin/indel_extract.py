#!/usr/bin/env python3
"""Extract PASS indels from one caller VCF (Pindel, Strelka, Mutect or generic), normalized
against the reference so that every caller's representation of an event is identical.

Rules (identical for every caller; the caller name is provenance only):
  * FILTER must be exactly PASS (a missing FILTER '.' counts only with --missing-filter-is-pass);
  * multiallelic records are split into one row per ALT;
  * an allele is kept only if REF and ALT are plain ACGT and differ in length (insertions,
    deletions and complex del-ins); SNVs, equal-length MNVs, symbolic and breakend alleles are
    dropped;
  * records on contigs outside --contigs are dropped;
  * REF must match the FASTA; mismatches are dropped and the run fails if they exceed
    --max-ref-mismatch-frac of candidate indel alleles (a sign the VCF is on another assembly);
  * alleles are left-aligned and trimmed (vt/bcftools-norm rule) before deduplication;
  * alleles whose length change is >= --max-length are dropped (structural variants).
"""
import argparse
import csv
import sys
from collections import Counter

import pysam

import geno_utils as gu
import indel_utils as iu

COLUMNS = ["chrom", "pos", "ref", "alt", "tumor_id", "normal_id", "caller"]
DROP_KEYS = ("not_pass", "not_indel", "mnv", "too_long", "off_contig", "ref_mismatch", "duplicate")


def extract(vcf, caller, tumor_id, normal_id, fasta, contigs=None, missing_filter_is_pass=False,
            max_ref_mismatch_frac=0.01, max_length=50):
    fa = pysam.FastaFile(fasta)
    dropped = Counter()
    rows, seen = [], set()
    candidates = normalized = 0
    for chrom, pos, ref, alt in gu.iter_pass_alleles(vcf, dropped, contigs, missing_filter_is_pass):
        kind = iu.allele_kind(ref, alt)
        if kind == "mnv":
            dropped["mnv"] += 1
            continue
        if kind != "indel":
            dropped["not_indel"] += 1
            continue
        candidates += 1
        try:
            fa_ref = fa.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper()
        except (KeyError, ValueError):
            fa_ref = None
        if fa_ref != ref:
            dropped["ref_mismatch"] += 1
            continue
        if abs(len(alt) - len(ref)) >= max_length:
            dropped["too_long"] += 1
            continue
        npos, nref, nalt = iu.normalize(lambda a, b, c=chrom: fa.fetch(c, max(a, 0), b), pos, ref, alt)
        normalized += (npos, nref, nalt) != (pos, ref, alt)
        key = (chrom, npos, nref, nalt)
        if key in seen:
            dropped["duplicate"] += 1
            continue
        seen.add(key)
        rows.append({"chrom": chrom, "pos": npos, "ref": nref, "alt": nalt, "tumor_id": tumor_id,
                     "normal_id": normal_id, "caller": caller})
    if candidates and dropped["ref_mismatch"] / candidates > max_ref_mismatch_frac:
        raise ValueError(
            f"{vcf}: REF mismatch against {fasta} for {dropped['ref_mismatch']}/{candidates} indel alleles "
            f"(> {max_ref_mismatch_frac}); wrong assembly?")
    stats = {"input": str(vcf), "caller": caller, "tumor_id": tumor_id, "normal_id": normal_id,
             "variant_type": "indel", "records": dropped["records"], "kept": len(rows),
             "normalized": normalized, "dropped": {k: dropped.get(k, 0) for k in DROP_KEYS}}
    return rows, stats


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--caller", required=True)
    ap.add_argument("--tumor-id", required=True)
    ap.add_argument("--normal-id", required=True)
    ap.add_argument("--reference", required=True, help="indexed FASTA (REF check and normalization)")
    ap.add_argument("--contigs", default=None, help="comma-separated contigs to keep (default: all)")
    ap.add_argument("--missing-filter-is-pass", action="store_true")
    ap.add_argument("--max-ref-mismatch-frac", type=float, default=0.01)
    ap.add_argument("--max-length", type=int, default=50,
                    help="drop alleles whose length change is >= this (default 50: structural variants)")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-stats", required=True)
    a = ap.parse_args(argv)
    if a.max_length < 1:
        ap.error("--max-length must be >= 1")
    rows, stats = extract(a.vcf, a.caller, a.tumor_id, a.normal_id, a.reference, gu.parse_contigs(a.contigs),
                          a.missing_filter_is_pass, a.max_ref_mismatch_frac, a.max_length)
    with open(a.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    gu.write_json(stats, a.out_stats)
    return 0


if __name__ == "__main__":
    sys.exit(main())
