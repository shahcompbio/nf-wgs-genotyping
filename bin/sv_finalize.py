#!/usr/bin/env python3
"""Gather svtyper and insertion-realignment chunk outputs (each checked against its chunk
manifest), join them onto the union adjacencies and clusters, and write the adjacency table, the
cluster table (+ csverve YAMLs), the adjacency table as BEDPE, a sorted/indexed svtyper VCF and a
QC JSON.

Every union adjacency is genotyped: by svtyper (DEL/DUP/BND records) or by sv_ins_count.py
(insertions, genotyper 'realign'). <S>_sv_adjacency_genotyping.csv.gz has one row per adjacency
(adj_id -> sv_id and its source calls). <S>_sv_genotyping.csv.gz keeps one row per cluster with
the v1.2 columns, taken from the cluster's representative adjacency (rep_adj_id), followed by a
summary of all its adjacencies.

<S>_sv_adjacency_genotyping.bedpe: the adjacency table in BEDPE layout, one row per adjacency.
Columns 1-10: chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2, where each end is
the 0-based half-open junction base [position - 1, position) (the CI is in ci*_lo/ci*_hi, not in
the interval), name = adj_id, score = alt_fragments, strands as in the tables (DEL +/-, DUP -/+).
Then every other adjacency-table column. Header line starts with '#'; missing values are '.';
rows in reference order of (chrom1, start1), then (chrom2, start2), name.

BND pairs are genotyped by svtyper from the same fragments on both records, so only the '_1'
record is read (never the sum of both ends). genotype_status:
  genotyped             counts present and some evidence
  no_evidence           counts present but no informative fragments/reads
  not_genotyped         svtyper emitted the record without counts (e.g. GT ./.)
  missing_from_output   the adjacency is absent from svtyper's output
  skipped_imprecise     IMPRECISE insertion (position too uncertain to build a haplotype)
  skipped_no_reference  insertion, but no reference FASTA was given
"""
import argparse
import json
import os
import re
import sys
from collections import Counter

import pandas as pd
import pysam

import geno_utils as gu

FIELD_MAP = {"AOF": "alt_fragments", "ROF": "ref_fragments", "FS": "alt_split_fragments",
             "FP": "alt_span_fragments", "FC": "alt_clip_fragments", "FX": "conflict_fragments",
             "GT": "svt_gt", "GQ": "svt_gq", "SQ": "svt_sq", "DP": "svt_dp", "RO": "svt_ro", "AO": "svt_ao",
             "RS": "svt_rs", "AS": "svt_as", "ASC": "svt_asc", "RP": "svt_rp", "AP": "svt_ap", "AB": "svt_ab"}
FLOAT_FIELDS = {"svt_sq", "svt_ab"}
COUNT_COLUMNS = [v for k, v in FIELD_MAP.items() if v != "svt_gt" and v not in FLOAT_FIELDS]
INS_COUNTS = ("ref_fragments", "alt_fragments", "uninformative_fragments", "conflict_fragments")
GENOTYPE_COLUMNS = ["alt_fragments", "ref_fragments", "alt_split_fragments", "alt_span_fragments",
                    "alt_clip_fragments", "conflict_fragments", "frag_vaf",
                    "svt_gt", "svt_gq", "svt_sq", "svt_dp", "svt_ro", "svt_ao", "svt_rs", "svt_as", "svt_asc",
                    "svt_rp", "svt_ap", "svt_ab", "genotype_status"]
PROVENANCE = ["n_members", "n_callers", "callers", "n_source_samples", "source_samples", "sources",
              "called_in_target"]
FINAL_COLUMNS = (["sample_id", "sv_id", "chromosome_1", "position_1", "strand_1", "chromosome_2", "position_2",
                  "strand_2", "sv_class", "svtyper_type", "sv_length"] + GENOTYPE_COLUMNS + PROVENANCE
                 + ["spread_1", "spread_2", "rep_caller", "rep_source",
                    # v1.3: the cluster's adjacencies
                    "rep_adj_id", "genotyper", "n_adjacencies", "n_adjacencies_genotyped", "n_adjacencies_with_alt",
                    "max_alt_fragments", "max_frag_vaf", "n_linked_clusters"])
ADJ_FINAL_COLUMNS = (["sample_id", "adj_id", "sv_id", "chromosome_1", "position_1", "strand_1", "chromosome_2",
                      "position_2", "strand_2", "sv_class", "genotyper", "svtyper_type", "sv_length", "ins_len",
                      "ins_seq_available", "imprecise", "ci1_lo", "ci1_hi", "ci2_lo", "ci2_hi"]
                     + GENOTYPE_COLUMNS[:6] + ["uninformative_fragments"] + GENOTYPE_COLUMNS[6:]
                     + ["is_cluster_rep"] + PROVENANCE)
_COORDINATES = ["adj_id", "chromosome_1", "position_1", "strand_1", "chromosome_2", "position_2", "strand_2"]
BEDPE_COLUMNS = (["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"]
                 + [c for c in ADJ_FINAL_COLUMNS if c not in _COORDINATES])
_COMMON_DTYPES = {c: "Int64" for c in COUNT_COLUMNS}
_COMMON_DTYPES.update({"sample_id": "category", "sv_id": "str", "chromosome_1": "category", "position_1": "int",
                       "strand_1": "category", "chromosome_2": "category", "position_2": "int",
                       "strand_2": "category", "sv_class": "category", "svtyper_type": "category",
                       "sv_length": "float", "frag_vaf": "float", "svt_gt": "str", "svt_sq": "float",
                       "svt_ab": "float", "genotype_status": "category", "n_members": "int", "n_callers": "int",
                       "callers": "str", "n_source_samples": "int", "source_samples": "str", "sources": "str",
                       "called_in_target": "bool", "genotyper": "category"})
DTYPES = {**_COMMON_DTYPES, "spread_1": "int", "spread_2": "int", "rep_caller": "category", "rep_source": "str",
          "rep_adj_id": "str", "n_adjacencies": "int", "n_adjacencies_genotyped": "int",
          "n_adjacencies_with_alt": "int", "max_alt_fragments": "Int64", "max_frag_vaf": "float",
          "n_linked_clusters": "int"}
ADJ_DTYPES = {**_COMMON_DTYPES, "adj_id": "str", "ins_len": "int", "ins_seq_available": "bool", "imprecise": "bool",
              "ci1_lo": "int", "ci1_hi": "int", "ci2_lo": "int", "ci2_hi": "int",
              "uninformative_fragments": "Int64", "is_cluster_rep": "bool"}
_CHUNK_RE = re.compile(r"(chunk_\d+)\.gt\.vcf$")
_INS_CHUNK_RE = re.compile(r"(ins_chunk_\d+)\.tsv$")
_STR = {"chromosome_1": str, "chromosome_2": str, "callers": str, "source_samples": str, "sources": str,
        "sv_id": str, "adj_id": str, "rep_adj_id": str, "ins_seq": str}


def _chunk_id(path, regex=_CHUNK_RE):
    m = regex.search(os.path.basename(path))
    if not m:
        raise ValueError(f"cannot derive chunk id from chunk output name: {path}")
    return m.group(1)


def _num(v, as_float=False):
    if v in (None, "", "."):
        return None
    try:
        return float(v) if as_float else int(float(v))
    except ValueError:
        return None


def parse_svtyper_vcf(path):
    """{adj_id: {final_column: value}} from one svtyper output (0-byte file allowed)."""
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
            adj_id = rid[:-2] if rid.endswith("_1") else rid
            fmt = dict(zip(t[8].split(":"), t[9].split(":"))) if len(t) > 9 else {}
            vals = {}
            for k, col in FIELD_MAP.items():
                v = fmt.get(k)
                vals[col] = (v if v not in (None, "") else None) if col == "svt_gt" else _num(v, col in FLOAT_FIELDS)
            out[adj_id] = vals
    return out


def _status(v):
    if v is None:
        return "missing_from_output"
    if v["alt_fragments"] is not None and v["ref_fragments"] is not None:
        return "genotyped" if (v["alt_fragments"] + v["ref_fragments"]) > 0 else "no_evidence"
    if v["svt_dp"] is not None:
        return "genotyped" if v["svt_dp"] > 0 else "no_evidence"
    return "not_genotyped"


def read_insertion_counts(ins_files, ins_manifest, expected_ids):
    """{adj_id: row} from sv_ins_count outputs; fails unless the files are exactly the manifest's
    chunks and every insertion adjacency has exactly one row."""
    if ins_manifest is None:
        if expected_ids:
            raise ValueError(f"{len(expected_ids)} insertion adjacencies but no insertion chunk manifest")
        return {}
    gu.check_gather(ins_manifest, [_chunk_id(p, _INS_CHUNK_RE) for p in ins_files])
    frames = [pd.read_csv(p, sep="\t", dtype={"adj_id": str, "ins_status": str}) for p in ins_files]
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["adj_id"])
    if df["adj_id"].duplicated().any():
        raise ValueError("duplicated insertion counts rows for some adjacencies")
    got = set(df["adj_id"])
    if got != set(expected_ids):
        raise ValueError(f"insertion adjacencies and counts rows disagree: missing="
                         f"{sorted(set(expected_ids) - got)[:3]} unexpected={sorted(got - set(expected_ids))[:3]}")
    return {r["adj_id"]: r for r in df.to_dict("records")}


def _int_or_none(v):
    return None if v is None or pd.isna(v) else int(v)


def _as_bool(s):
    return s.map(lambda v: v if isinstance(v, bool) else str(v).strip().lower() in ("true", "1")).astype(bool)


def genotype_adjacencies(adjs, parsed, ins):
    rows = []
    empty = {c: None for c in FIELD_MAP.values()}
    for a in adjs.to_dict("records"):
        if a["genotyper"] == "realign":
            r = ins[a["adj_id"]]
            g = {**empty, **{c: _int_or_none(r[c]) for c in INS_COUNTS}, "genotype_status": r["ins_status"]}
        else:
            v = parsed.get(a["adj_id"])
            g = {**(v or empty), "uninformative_fragments": None, "genotype_status": _status(v)}
        rows.append({**a, **g})
    cols = list(adjs.columns) + list(FIELD_MAP.values()) + ["uninformative_fragments", "genotype_status"]
    df = pd.DataFrame(rows, columns=cols)
    for c in COUNT_COLUMNS + ["uninformative_fragments"]:
        df[c] = pd.array(df[c].tolist(), dtype="Int64") if len(df) else pd.array([], dtype="Int64")
    for c in FLOAT_FIELDS:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype(float)
    df["frag_vaf"] = _vaf(df)
    for c in ("called_in_target", "is_cluster_rep", "ins_seq_available", "imprecise"):
        df[c] = _as_bool(df[c]) if len(df) else df[c].astype(bool)
    return df


def _vaf(df):
    if not len(df):
        return pd.Series(dtype=float)
    alt = df["alt_fragments"].astype("Float64")
    tot = alt + df["ref_fragments"].astype("Float64")
    return (alt / tot.where(tot > 0)).astype(float)


def summarize_clusters(union, adj_df):
    geno = adj_df[["adj_id"] + [c for c in GENOTYPE_COLUMNS]].rename(columns={"adj_id": "rep_adj_id"})
    df = union.merge(geno, on="rep_adj_id", how="left", validate="one_to_one")
    a = adj_df.assign(_genotyped=adj_df["genotype_status"] == "genotyped",
                      _with_alt=adj_df["alt_fragments"].fillna(0) >= 1)
    summary = a.groupby("sv_id").agg(n_adjacencies_genotyped=("_genotyped", "sum"),
                                     n_adjacencies_with_alt=("_with_alt", "sum"),
                                     max_alt_fragments=("alt_fragments", "max"),
                                     max_frag_vaf=("frag_vaf", "max")).reset_index()
    df = df.merge(summary, on="sv_id", how="left", validate="one_to_one")
    for c in ("n_adjacencies_genotyped", "n_adjacencies_with_alt"):
        df[c] = df[c].fillna(0).astype(int)
    df["max_alt_fragments"] = pd.array(df["max_alt_fragments"].tolist(), dtype="Int64") if len(df) else \
        pd.array([], dtype="Int64")
    df["max_frag_vaf"] = df["max_frag_vaf"].astype(float)
    for c in COUNT_COLUMNS:
        df[c] = df[c].astype("Int64")
    df["called_in_target"] = _as_bool(df["called_in_target"]) if len(df) else df["called_in_target"].astype(bool)
    return df


def finalize(union_bedpe, adjacencies_tsv, gt_files, manifest, sample_id, stats_files=(), ins_files=(),
             ins_manifest=None):
    """-> (cluster table, adjacency table, qc)."""
    gu.check_gather(manifest, [_chunk_id(p) for p in gt_files])
    union = pd.read_csv(union_bedpe, sep="\t", dtype=_STR)
    adjs = pd.read_csv(adjacencies_tsv, sep="\t", dtype=_STR)
    ins_ids = adjs.loc[adjs["genotyper"] == "realign", "adj_id"].tolist()
    ins = read_insertion_counts(list(ins_files), ins_manifest, ins_ids)
    parsed = {}
    for p in gt_files:
        parsed.update(parse_svtyper_vcf(p))
    adj_df = genotype_adjacencies(adjs, parsed, ins)
    adj_df["sample_id"] = sample_id
    clusters = summarize_clusters(union, adj_df)
    clusters["sample_id"] = sample_id

    in_target = clusters[clusters["called_in_target"]]
    ins_adj = adj_df[adj_df["genotyper"] == "realign"]
    qc = {
        "type": "sv_genotyping",
        "sample_id": sample_id,
        "union_size": int(len(clusters)),
        "called_in_target": int(len(in_target)),
        # representative adjacency (as in v1.2) and any adjacency of the cluster
        "positive_control_fraction": (float((in_target["alt_fragments"].fillna(0) >= 1).mean())
                                      if len(in_target) else None),
        "positive_control_fraction_any_adjacency": (float((in_target["n_adjacencies_with_alt"] >= 1).mean())
                                                    if len(in_target) else None),
        "svs_with_alt_support": int((clusters["alt_fragments"].fillna(0) >= 1).sum()) if len(clusters) else 0,
        "status_counts": dict(sorted(Counter(clusters["genotype_status"]).items())),
        "n_adjacencies": int(len(adj_df)),
        "adjacency_status_counts": dict(sorted(Counter(adj_df["genotype_status"]).items())),
        "clusters_with_multiple_adjacencies": int((clusters["n_adjacencies"] > 1).sum()),
        "clusters_linked_to_other_clusters": int((clusters["n_linked_clusters"] > 0).sum()),
        "insertions": {"adjacencies": int(len(ins_adj)),
                       "sequence_resolved": int(ins_adj["ins_seq_available"].sum()),
                       "length_only": int((~ins_adj["ins_seq_available"]).sum()),
                       "status_counts": dict(sorted(Counter(ins_adj["genotype_status"]).items()))},
        "chunks": len(gt_files),
        "insertion_chunks": len(ins_files),
        "inputs": [json.load(open(p)) for p in stats_files],
    }
    return clusters[FINAL_COLUMNS], adj_df[ADJ_FINAL_COLUMNS], qc


def bedpe_table(adj_df, order):
    """The final adjacency table in BEDPE layout (BEDPE_COLUMNS), sorted by reference order."""
    pos1, pos2 = adj_df["position_1"].astype(int), adj_df["position_2"].astype(int)
    head = pd.DataFrame({"chrom1": adj_df["chromosome_1"], "start1": pos1 - 1, "end1": pos1,
                         "chrom2": adj_df["chromosome_2"], "start2": pos2 - 1, "end2": pos2,
                         "name": adj_df["adj_id"], "score": adj_df["alt_fragments"],
                         "strand1": adj_df["strand_1"], "strand2": adj_df["strand_2"]})
    df = pd.concat([head, adj_df[BEDPE_COLUMNS[10:]]], axis=1)[BEDPE_COLUMNS]
    key = gu.site_sort_key(order)
    rows = list(zip(df["chrom1"], df["end1"], df["chrom2"], df["end2"], df["name"]))
    idx = sorted(range(len(df)), key=lambda i: (key(rows[i][0], rows[i][1]), key(rows[i][2], rows[i][3]),
                                                 rows[i][4]))
    return df.iloc[idx].reset_index(drop=True)


def write_bedpe(df, path):
    """Header line '#chrom1<TAB>start1...' then tab-separated rows, missing values as '.'."""
    with open(path, "w") as fh:
        fh.write("#" + "\t".join(df.columns) + "\n")
        df.to_csv(fh, sep="\t", header=False, index=False, na_rep=".")


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
    ap.add_argument("--union", required=True, help="cluster union (<S>_union_svs.bedpe.gz)")
    ap.add_argument("--adjacencies", required=True, help="<S>_union_sv_adjacencies.tsv.gz")
    ap.add_argument("--genotyped", nargs="+", required=True, help="svtyper outputs named <chunk_id>.gt.vcf")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--ins-counts", nargs="*", default=[], help="sv_ins_count outputs named *<ins_chunk_id>.tsv")
    ap.add_argument("--ins-manifest", required=True)
    ap.add_argument("--stats", nargs="*", default=[])
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--reference", default=None)
    ap.add_argument("--contigs", default=None)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)
    clusters, adjs, qc = finalize(a.union, a.adjacencies, a.genotyped, a.manifest, a.sample_id, a.stats,
                                  a.ins_counts, a.ins_manifest)
    order = gu.contig_order(fai=f"{a.reference}.fai" if a.reference else None, contigs=gu.parse_contigs(a.contigs))
    gu.write_table(clusters, os.path.join(a.outdir, f"{a.sample_id}_sv_genotyping.csv.gz"), DTYPES)
    gu.write_table(adjs, os.path.join(a.outdir, f"{a.sample_id}_sv_adjacency_genotyping.csv.gz"), ADJ_DTYPES)
    write_bedpe(bedpe_table(adjs, order), os.path.join(a.outdir, f"{a.sample_id}_sv_adjacency_genotyping.bedpe"))
    gu.write_json(qc, os.path.join(a.outdir, f"{a.sample_id}_genotyping_qc.json"))
    concat_vcfs(sorted(a.genotyped), os.path.join(a.outdir, f"{a.sample_id}_svtyper.vcf"), order, a.sample_id)
    return 0


if __name__ == "__main__":
    sys.exit(main())
