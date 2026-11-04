#!/usr/bin/env python3
"""Gather svtyper chunk outputs (checked against the chunk manifest), join them onto the SV union,
and write the final table + YAML, a sorted/indexed svtyper VCF and a QC JSON.

BND pairs are genotyped by svtyper from the same fragments on both records, so only the '_1'
record is read (never the sum of both ends). genotype_status:
  genotyped           counts present and some evidence
  no_evidence         counts present but no informative fragments/reads
  not_genotyped       svtyper emitted the record without counts (e.g. GT ./.)
  missing_from_output the union SV is absent from svtyper's output
"""
import argparse
import json
import os
import re
import sys
from collections import Counter

import numpy as np
import pandas as pd
import pysam

import geno_utils as gu

FIELD_MAP = {"AOF": "alt_fragments", "ROF": "ref_fragments", "FS": "alt_split_fragments",
             "FP": "alt_span_fragments", "FC": "alt_clip_fragments", "FX": "conflict_fragments",
             "GT": "svt_gt", "GQ": "svt_gq", "SQ": "svt_sq", "DP": "svt_dp", "RO": "svt_ro", "AO": "svt_ao",
             "RS": "svt_rs", "AS": "svt_as", "ASC": "svt_asc", "RP": "svt_rp", "AP": "svt_ap", "AB": "svt_ab"}
FLOAT_FIELDS = {"svt_sq", "svt_ab"}
COUNT_COLUMNS = [v for k, v in FIELD_MAP.items() if v != "svt_gt" and v not in FLOAT_FIELDS]
FINAL_COLUMNS = ["sample_id", "sv_id", "chromosome_1", "position_1", "strand_1", "chromosome_2", "position_2",
                 "strand_2", "sv_class", "svtyper_type", "sv_length",
                 "alt_fragments", "ref_fragments", "alt_split_fragments", "alt_span_fragments",
                 "alt_clip_fragments", "conflict_fragments", "frag_vaf",
                 "svt_gt", "svt_gq", "svt_sq", "svt_dp", "svt_ro", "svt_ao", "svt_rs", "svt_as", "svt_asc",
                 "svt_rp", "svt_ap", "svt_ab", "genotype_status",
                 "n_members", "n_callers", "callers", "n_source_samples", "source_samples", "sources",
                 "called_in_target", "spread_1", "spread_2", "rep_caller", "rep_source"]
DTYPES = {c: "Int64" for c in COUNT_COLUMNS}
DTYPES.update({"sample_id": "category", "sv_id": "str", "chromosome_1": "category", "position_1": "int",
               "strand_1": "category", "chromosome_2": "category", "position_2": "int", "strand_2": "category",
               "sv_class": "category", "svtyper_type": "category", "sv_length": "float", "frag_vaf": "float",
               "svt_gt": "str", "svt_sq": "float", "svt_ab": "float", "genotype_status": "category",
               "n_members": "int", "n_callers": "int", "callers": "str", "n_source_samples": "int",
               "source_samples": "str", "sources": "str", "called_in_target": "bool", "spread_1": "int",
               "spread_2": "int", "rep_caller": "category", "rep_source": "str"})
_CHUNK_RE = re.compile(r"(chunk_\d+)\.gt\.vcf$")


def _chunk_id(path):
    m = _CHUNK_RE.search(os.path.basename(path))
    if not m:
        raise ValueError(f"cannot derive chunk id from svtyper output name: {path}")
    return m.group(1)


def _num(v, as_float=False):
    if v in (None, "", "."):
        return None
    try:
        return float(v) if as_float else int(float(v))
    except ValueError:
        return None


def parse_svtyper_vcf(path):
    """{sv_id: {final_column: value}} from one svtyper output (0-byte file allowed)."""
    out = {}
    if os.path.getsize(path) == 0:
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            t = line.rstrip("\n").split("\t")
            rid = t[2]
            if rid.endswith("_2"):
                continue
            sv_id = rid[:-2] if rid.endswith("_1") else rid
            fmt = dict(zip(t[8].split(":"), t[9].split(":"))) if len(t) > 9 else {}
            vals = {}
            for k, col in FIELD_MAP.items():
                v = fmt.get(k)
                vals[col] = (v if v not in (None, "") else None) if col == "svt_gt" else _num(v, col in FLOAT_FIELDS)
            out[sv_id] = vals
    return out


def _status(v):
    if v is None:
        return "missing_from_output"
    if v["alt_fragments"] is not None and v["ref_fragments"] is not None:
        return "genotyped" if (v["alt_fragments"] + v["ref_fragments"]) > 0 else "no_evidence"
    if v["svt_dp"] is not None:
        return "genotyped" if v["svt_dp"] > 0 else "no_evidence"
    return "not_genotyped"


def finalize(union_bedpe, gt_files, manifest, sample_id, stats_files=()):
    gu.check_gather(manifest, [_chunk_id(p) for p in gt_files])
    union = pd.read_csv(union_bedpe, sep="\t", dtype={"chromosome_1": str, "chromosome_2": str, "callers": str,
                                                      "source_samples": str, "sources": str, "sv_id": str})
    parsed = {}
    for p in gt_files:
        parsed.update(parse_svtyper_vcf(p))
    rows = []
    for u in union.to_dict("records"):
        v = parsed.get(u["sv_id"])
        row = {**u, **(v or {c: None for c in FIELD_MAP.values()}), "genotype_status": _status(v)}
        rows.append(row)
    df = pd.DataFrame(rows, columns=list(union.columns) + list(FIELD_MAP.values()) + ["genotype_status"])
    for c in COUNT_COLUMNS:
        df[c] = pd.array(df[c].tolist(), dtype="Int64") if len(df) else pd.array([], dtype="Int64")
    for c in FLOAT_FIELDS:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype(float)
    alt = df["alt_fragments"].astype("Float64")
    tot = alt + df["ref_fragments"].astype("Float64")
    df["frag_vaf"] = (alt / tot.where(tot > 0)).astype(float) if len(df) else pd.Series(dtype=float)
    df["sample_id"] = sample_id
    df["called_in_target"] = df["called_in_target"].astype(bool)
    in_target = df[df["called_in_target"]]
    supported = (in_target["alt_fragments"].fillna(0) >= 1) if len(in_target) else pd.Series(dtype=bool)
    qc = {
        "type": "sv_genotyping",
        "sample_id": sample_id,
        "union_size": int(len(df)),
        "called_in_target": int(len(in_target)),
        "positive_control_fraction": float(supported.mean()) if len(in_target) else None,
        "svs_with_alt_support": int((df["alt_fragments"].fillna(0) >= 1).sum()) if len(df) else 0,
        "status_counts": dict(sorted(Counter(df["genotype_status"]).items())),
        "chunks": len(gt_files),
        "inputs": [json.load(open(p)) for p in stats_files],
    }
    return df[FINAL_COLUMNS], qc


def concat_vcfs(gt_files, out_vcf, order, sample_id):
    header, recs = None, []
    for p in gt_files:
        if os.path.getsize(p) == 0:
            continue
        with open(p) as fh:
            lines = fh.read().splitlines()
        if header is None:
            header = [l for l in lines if l.startswith("#")]
        recs.extend(l.split("\t") for l in lines if l and not l.startswith("#"))
    if header is None:
        header = ["##fileformat=VCFv4.2",
                  "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_id])]
    key = gu.site_sort_key(order)
    recs.sort(key=lambda r: (key(r[0], r[1]), r[2]))
    with open(out_vcf, "w") as fh:
        fh.write("\n".join(header) + "\n")
        for r in recs:
            fh.write("\t".join(r) + "\n")
    pysam.tabix_index(out_vcf, preset="vcf", force=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--union", required=True)
    ap.add_argument("--genotyped", nargs="+", required=True, help="svtyper outputs named <chunk_id>.gt.vcf")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--stats", nargs="*", default=[])
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--reference", default=None)
    ap.add_argument("--contigs", default=None)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)
    df, qc = finalize(a.union, a.genotyped, a.manifest, a.sample_id, a.stats)
    gu.write_table(df, os.path.join(a.outdir, f"{a.sample_id}_sv_genotyping.csv.gz"), DTYPES)
    gu.write_json(qc, os.path.join(a.outdir, f"{a.sample_id}_genotyping_qc.json"))
    order = gu.contig_order(fai=f"{a.reference}.fai" if a.reference else None, contigs=gu.parse_contigs(a.contigs))
    concat_vcfs(sorted(a.genotyped), os.path.join(a.outdir, f"{a.sample_id}_svtyper.vcf"), order, a.sample_id)
    return 0


if __name__ == "__main__":
    sys.exit(main())
