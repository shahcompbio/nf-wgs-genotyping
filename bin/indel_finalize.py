#!/usr/bin/env python3
"""Gather indel counts (explicitly checked against the chunk manifest), join them onto the indel
union, and write the final indel table, its csverve-style YAML, a QC JSON, and a single-sample
VCF of the same rows (bgzipped + tabix index; FORMAT counts only, no GT)."""
import argparse
import json
import os
import sys

import numpy as np

import geno_utils as gu
import snv_finalize

COUNT_COLUMNS = ("ref_count", "alt_count", "other_count", "uninformative_count", "conflict_count")
FINAL_COLUMNS = ["sample_id", "chromosome", "position", "ref", "alt", "indel_length", "ref_count",
                 "alt_count", "other_count", "depth", "vaf", "uninformative_count", "conflict_count",
                 "locus_id", "var_start", "var_end", "n_callers", "callers", "n_source_samples",
                 "source_samples", "sources", "called_in_target"]
DTYPES = {"sample_id": "category", "chromosome": "category", "position": "int", "ref": "str", "alt": "str",
          "indel_length": "int", "ref_count": "int", "alt_count": "int", "other_count": "int", "depth": "int",
          "vaf": "float", "uninformative_count": "int", "conflict_count": "int", "locus_id": "str",
          "var_start": "int", "var_end": "int", "n_callers": "int", "callers": "str", "n_source_samples": "int",
          "source_samples": "str", "sources": "str", "called_in_target": "bool"}
# Added to the union VCF's header (contigs + INFO CALLERS/SAMPLES/NCALLERS/NSAMPLES/INTARGET).
VCF_HEADER_LINES = (
    '##INFO=<ID=SOURCES,Number=.,Type=String,Description="Source calls as <tumor>_vs_<normal>:<caller>">',
    '##INFO=<ID=LOCUS,Number=1,Type=String,Description="Locus genotyped jointly: union alleles whose '
    'repeat-aware footprints overlap (chrom:start-end, 1-based)">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Fragments supporting ref and alt haplotypes '
    '(realignment; each fragment counted once; filters in docs/params.md)">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Informative fragments: ref + alt + other">',
    '##FORMAT=<ID=AF,Number=A,Type=Float,Description="alt / (ref + alt); missing when ref + alt = 0">',
    '##FORMAT=<ID=OC,Number=1,Type=Integer,Description="Fragments supporting another union allele of the locus">',
    '##FORMAT=<ID=UF,Number=1,Type=Integer,Description="Fragments reaching the locus that support no '
    'haplotype (e.g. ending inside the repeat); per locus">',
)


def finalize(union_tsv, counts_files, manifest, sample_id, stats_files=()):
    df = snv_finalize.join_counts(union_tsv, counts_files, manifest, COUNT_COLUMNS, {"locus_id": str})
    df["depth"] = df["ref_count"] + df["alt_count"] + df["other_count"]
    denom = df["ref_count"] + df["alt_count"]
    df["vaf"] = np.where(denom > 0, df["alt_count"] / denom.where(denom > 0, 1), np.nan)
    df["sample_id"] = sample_id
    df["called_in_target"] = df["called_in_target"].astype(bool)
    in_target = df[df["called_in_target"]]
    loci = df.drop_duplicates("locus_id")              # uninformative/conflict are per locus
    qc = {
        "type": "indel_genotyping",
        "sample_id": sample_id,
        "union_size": int(len(df)),
        "loci": int(len(loci)),
        "called_in_target": int(len(in_target)),
        "positive_control_fraction": (float((in_target["alt_count"] >= 1).mean()) if len(in_target) else None),
        "sites_with_alt_support": int((df["alt_count"] >= 1).sum()),
        "median_depth": float(df["depth"].median()) if len(df) else None,
        "uninformative_fragments": int(loci["uninformative_count"].sum()),
        "conflict_fragments": int(loci["conflict_count"].sum()),
        "chunks": len(counts_files),
        "inputs": [json.load(open(p)) for p in stats_files],
    }
    return df[FINAL_COLUMNS], qc


def _indel_fields(rec, row):
    rec.info["LOCUS"] = str(row.locus_id)
    rec.samples[row.sample_id]["UF"] = int(row.uninformative_count)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--union", required=True)
    ap.add_argument("--union-vcf", required=True, help="indel_union VCF (its header is reused)")
    ap.add_argument("--counts", nargs="+", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--stats", nargs="*", default=[])
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)
    df, qc = finalize(a.union, a.counts, a.manifest, a.sample_id, a.stats)
    gu.write_table(df, os.path.join(a.outdir, f"{a.sample_id}_indel_genotyping.csv.gz"), DTYPES)
    gu.write_json(qc, os.path.join(a.outdir, f"{a.sample_id}_indel_qc.json"))
    snv_finalize.write_vcf(df, a.union_vcf, os.path.join(a.outdir, f"{a.sample_id}_indel_genotyping.vcf.gz"),
                           a.sample_id, header_lines=VCF_HEADER_LINES, extra=_indel_fields)
    return 0


if __name__ == "__main__":
    sys.exit(main())
