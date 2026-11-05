#!/usr/bin/env python3
"""Gather SNV counts (explicitly checked against the chunk manifest), join them onto the union,
and write the final genotyping table, its csverve-style YAML, a QC JSON, and a single-sample VCF
of the same rows (bgzipped + tabix index; FORMAT counts only, no GT)."""
import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd
import pysam

import geno_utils as gu

KEYS = ["chromosome", "position", "ref", "alt"]
FINAL_COLUMNS = ["sample_id", "chromosome", "position", "ref", "alt", "ref_count", "alt_count",
                 "other_count", "depth", "vaf", "n_callers", "callers", "n_source_samples",
                 "source_samples", "sources", "called_in_target"]
DTYPES = {"sample_id": "category", "chromosome": "category", "position": "int", "ref": "category",
          "alt": "category", "ref_count": "int", "alt_count": "int", "other_count": "int", "depth": "int",
          "vaf": "float", "n_callers": "int", "callers": "str", "n_source_samples": "int",
          "source_samples": "str", "sources": "str", "called_in_target": "bool"}
_CHUNK_RE = re.compile(r"counts_(chunk_\d+)\.tsv$")
# Added to the union VCF's header (contigs + INFO CALLERS/SAMPLES/NCALLERS/NSAMPLES/INTARGET).
# No GT: the pipeline counts fragments and makes no genotype call.
VCF_HEADER_LINES = (
    '##INFO=<ID=SOURCES,Number=.,Type=String,Description="Source calls as <tumor>_vs_<normal>:<caller>">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Fragments supporting ref and alt '
    '(each fragment counted once; filters in docs/params.md)">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Fragments counted: ref + alt + other">',
    '##FORMAT=<ID=AF,Number=A,Type=Float,Description="alt / (ref + alt); missing when ref + alt = 0">',
    '##FORMAT=<ID=OC,Number=1,Type=Integer,Description="Fragments with a base other than ref or alt">',
)


def _chunk_id(path):
    m = _CHUNK_RE.search(os.path.basename(path))
    if not m:
        raise ValueError(f"cannot derive chunk id from counts file name: {path}")
    return m.group(1)


def join_counts(union_tsv, counts_files, manifest, count_columns, union_dtypes=None):
    """Explicit gather: fail unless the counts files are exactly the manifest's chunks and every
    union row has exactly one counts row (joined on chromosome, position, ref, alt)."""
    gu.check_gather(manifest, [_chunk_id(p) for p in counts_files])
    dt = {"chromosome": str, "ref": str, "alt": str}
    union = pd.read_csv(union_tsv, sep="\t", dtype={**dt, "callers": str, "source_samples": str, "sources": str,
                                                    **(union_dtypes or {})})
    counts = pd.concat([pd.read_csv(p, sep="\t", dtype=dt) for p in counts_files], ignore_index=True)
    if counts.duplicated(KEYS).any():
        raise ValueError("duplicated counts rows for some union sites")
    merged = union.merge(counts, on=KEYS, how="outer", indicator=True)
    bad = merged[merged["_merge"] != "both"]
    if len(bad):
        raise ValueError(f"union sites and counts rows disagree for {len(bad)} sites, e.g. "
                         f"{bad[KEYS].head(3).values.tolist()}")
    df = union.merge(counts, on=KEYS, how="left")
    for c in count_columns:
        df[c] = df[c].astype(int)
    return df


def finalize(union_tsv, counts_files, manifest, sample_id, stats_files=()):
    df = join_counts(union_tsv, counts_files, manifest, ("ref_count", "alt_count", "other_count", "conflict_count"))
    df["depth"] = df["ref_count"] + df["alt_count"] + df["other_count"]
    denom = df["ref_count"] + df["alt_count"]
    df["vaf"] = np.where(denom > 0, df["alt_count"] / denom.where(denom > 0, 1), np.nan)
    df["sample_id"] = sample_id
    df["called_in_target"] = df["called_in_target"].astype(bool)
    in_target = df[df["called_in_target"]]
    qc = {
        "type": "snv_genotyping",
        "sample_id": sample_id,
        "union_size": int(len(df)),
        "called_in_target": int(len(in_target)),
        "positive_control_fraction": (float((in_target["alt_count"] >= 1).mean()) if len(in_target) else None),
        "sites_with_alt_support": int((df["alt_count"] >= 1).sum()),
        "median_depth": float(df["depth"].median()) if len(df) else None,
        "conflict_fragments": int(df["conflict_count"].sum()),
        "chunks": len(counts_files),
        "inputs": [json.load(open(p)) for p in stats_files],
    }
    return df[FINAL_COLUMNS], qc


def write_vcf(df, union_vcf, path, sample_id, header_lines=VCF_HEADER_LINES, extra=None):
    """Write the final table `df` as a bgzipped, tabix-indexed VCF at `path` (must end in .vcf.gz):
    one record per row in the same (union, contig-sorted) order, header from `union_vcf` plus
    `header_lines`. `extra(record, row)` sets further INFO/FORMAT fields."""
    if not path.endswith(".vcf.gz"):
        raise ValueError(f"VCF output must end in .vcf.gz: {path}")
    with pysam.VariantFile(union_vcf) as src:
        header = src.header.copy()
    for line in header_lines:
        header.add_line(line)
    header.add_sample(sample_id)
    plain = path[:-3]
    with pysam.VariantFile(plain, "w", header=header) as out:
        for r in df.itertuples(index=False):
            rec = out.new_record(contig=str(r.chromosome), start=int(r.position) - 1,
                                 alleles=(str(r.ref), str(r.alt)), filter="PASS")
            rec.info["CALLERS"] = tuple(r.callers.split(";"))
            rec.info["SAMPLES"] = tuple(r.source_samples.split(";"))
            rec.info["NCALLERS"] = int(r.n_callers)
            rec.info["NSAMPLES"] = int(r.n_source_samples)
            rec.info["SOURCES"] = tuple(r.sources.split(";"))
            if r.called_in_target:
                rec.info["INTARGET"] = True
            s = rec.samples[sample_id]
            s["AD"] = (int(r.ref_count), int(r.alt_count))
            s["DP"] = int(r.depth)
            s["AF"] = (None,) if np.isnan(r.vaf) else (float(r.vaf),)
            s["OC"] = int(r.other_count)
            if extra is not None:
                extra(rec, r)
            out.write(rec)
    pysam.tabix_index(plain, preset="vcf", force=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--union", required=True)
    ap.add_argument("--union-vcf", required=True, help="snv_union VCF (its header is reused)")
    ap.add_argument("--counts", nargs="+", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--stats", nargs="*", default=[])
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)
    df, qc = finalize(a.union, a.counts, a.manifest, a.sample_id, a.stats)
    gu.write_table(df, os.path.join(a.outdir, f"{a.sample_id}_snv_genotyping.csv.gz"), DTYPES)
    gu.write_json(qc, os.path.join(a.outdir, f"{a.sample_id}_genotyping_qc.json"))
    write_vcf(df, a.union_vcf, os.path.join(a.outdir, f"{a.sample_id}_snv_genotyping.vcf.gz"), a.sample_id)
    return 0


if __name__ == "__main__":
    sys.exit(main())
