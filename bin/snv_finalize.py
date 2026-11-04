#!/usr/bin/env python3
"""Gather SNV counts (explicitly checked against the chunk manifest), join them onto the union,
and write the final genotyping table, its csverve-style YAML, and a QC JSON."""
import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd

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


def _chunk_id(path):
    m = _CHUNK_RE.search(os.path.basename(path))
    if not m:
        raise ValueError(f"cannot derive chunk id from counts file name: {path}")
    return m.group(1)


def finalize(union_tsv, counts_files, manifest, sample_id, stats_files=()):
    gu.check_gather(manifest, [_chunk_id(p) for p in counts_files])
    dt = {"chromosome": str, "ref": str, "alt": str}
    union = pd.read_csv(union_tsv, sep="\t", dtype={**dt, "callers": str, "source_samples": str, "sources": str})
    counts = pd.concat([pd.read_csv(p, sep="\t", dtype=dt) for p in counts_files], ignore_index=True)
    if counts.duplicated(KEYS).any():
        raise ValueError("duplicated counts rows for some union sites")
    merged = union.merge(counts, on=KEYS, how="outer", indicator=True)
    bad = merged[merged["_merge"] != "both"]
    if len(bad):
        raise ValueError(f"union sites and counts rows disagree for {len(bad)} sites, e.g. "
                         f"{bad[KEYS].head(3).values.tolist()}")
    df = union.merge(counts, on=KEYS, how="left")
    for c in ("ref_count", "alt_count", "other_count", "conflict_count"):
        df[c] = df[c].astype(int)
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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--union", required=True)
    ap.add_argument("--counts", nargs="+", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--stats", nargs="*", default=[])
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)
    df, qc = finalize(a.union, a.counts, a.manifest, a.sample_id, a.stats)
    gu.write_table(df, os.path.join(a.outdir, f"{a.sample_id}_snv_genotyping.csv.gz"), DTYPES)
    gu.write_json(qc, os.path.join(a.outdir, f"{a.sample_id}_genotyping_qc.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
