#!/usr/bin/env python3
"""Cluster normalized SV junctions across callers and source samples, write the union with
provenance, an svtyper-compatible VCF, and chunk VCFs for scatter/gather.

Clustering: junctions are bucketed by (chrom1, chrom2, strand1, strand2); within a bucket two
junctions are linked when |pos1 - pos1'| <= slop and |pos2 - pos2'| <= slop; clusters are the
connected components (union-find), so the result does not depend on input order.
Representative: precise before imprecise, then --rep-priority caller order, then narrowest CI,
then (tumor_id, record_id). Cluster ids SV000001... follow genomic order of the representative.

svtyper VCF: intrachromosomal (+,-) -> <DEL> (POS=pos1, END=pos2-1); intrachromosomal (-,+) ->
<DUP> (POS=pos1-1, END=pos2); everything else -> a BND pair <id>_1/<id>_2 with MATEID. svtyper adds
1 to the POS of reverse ('-') ends, so '-' ends are written one base left of the junction. Every
record carries CIPOS (and CIEND for DEL/DUP), which svtyper requires.
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

import geno_utils as gu
from sv_normalize import COLUMNS as MEMBER_COLUMNS

UNION_COLUMNS = ["sv_id", "chromosome_1", "position_1", "strand_1", "chromosome_2", "position_2", "strand_2",
                 "ci1_lo", "ci1_hi", "ci2_lo", "ci2_hi", "imprecise", "sv_class", "svtyper_type", "sv_length",
                 "n_members", "n_callers", "callers", "n_source_samples", "source_samples", "sources",
                 "called_in_target", "spread_1", "spread_2", "rep_caller", "rep_source"]


def svtyper_type(chrom1, strand1, chrom2, strand2):
    if chrom1 == chrom2 and (strand1, strand2) == ("+", "-"):
        return "DEL"
    if chrom1 == chrom2 and (strand1, strand2) == ("-", "+"):
        return "DUP"
    return "BND"


def _find(parent, i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


def cluster_union(rows, sample_id, order, slop, rep_priority):
    rows = [dict(r) for r in rows]
    for r in rows:
        r["pos1"], r["pos2"] = int(r["pos1"]), int(r["pos2"])
        for k in ("ci1_lo", "ci1_hi", "ci2_lo", "ci2_hi"):
            r[k] = int(r[k])
        r["imprecise"] = str(r["imprecise"]).lower() in ("true", "1")
    n = len(rows)
    parent = list(range(n))
    buckets = defaultdict(list)
    for i, r in enumerate(rows):
        buckets[(r["chrom1"], r["chrom2"], r["strand1"], r["strand2"])].append(i)
    for idxs in buckets.values():
        idxs.sort(key=lambda i: (rows[i]["pos1"], rows[i]["pos2"]))
        for ai, a in enumerate(idxs):
            for b in idxs[ai + 1:]:
                if rows[b]["pos1"] - rows[a]["pos1"] > slop:
                    break
                if abs(rows[b]["pos2"] - rows[a]["pos2"]) <= slop:
                    pa, pb = _find(parent, a), _find(parent, b)
                    if pa != pb:
                        parent[max(pa, pb)] = min(pa, pb)
    comps = defaultdict(list)
    for i in range(n):
        comps[_find(parent, i)].append(i)

    prio = {c: k for k, c in enumerate(rep_priority)}

    def rep_key(r):
        width = (r["ci1_hi"] - r["ci1_lo"]) + (r["ci2_hi"] - r["ci2_lo"])
        return (r["imprecise"], prio.get(r["caller"], len(prio)), width, r["tumor_id"], str(r["record_id"]),
                r["pos1"], r["pos2"], r["normal_id"])

    def end_key(chrom, pos):
        return gu.site_sort_key(order)(chrom, pos)

    clusters = []
    for members in comps.values():
        ms = [rows[i] for i in members]
        rep = min(ms, key=rep_key)
        clusters.append((rep, ms))
    clusters.sort(key=lambda c: (end_key(c[0]["chrom1"], c[0]["pos1"]), end_key(c[0]["chrom2"], c[0]["pos2"]),
                                 c[0]["strand1"], c[0]["strand2"]))
    out, member_rows = [], []
    for k, (rep, ms) in enumerate(clusters, start=1):
        sv_id = f"SV{k:06d}"
        callers = sorted({m["caller"] for m in ms})
        samples = sorted({m["tumor_id"] for m in ms})
        same = rep["chrom1"] == rep["chrom2"]
        out.append({
            "sv_id": sv_id, "chromosome_1": rep["chrom1"], "position_1": rep["pos1"], "strand_1": rep["strand1"],
            "chromosome_2": rep["chrom2"], "position_2": rep["pos2"], "strand_2": rep["strand2"],
            "ci1_lo": rep["ci1_lo"], "ci1_hi": rep["ci1_hi"], "ci2_lo": rep["ci2_lo"], "ci2_hi": rep["ci2_hi"],
            "imprecise": rep["imprecise"], "sv_class": rep["sv_class"],
            "svtyper_type": svtyper_type(rep["chrom1"], rep["strand1"], rep["chrom2"], rep["strand2"]),
            "sv_length": float(rep["pos2"] - rep["pos1"]) if same else np.nan,
            "n_members": len(ms), "n_callers": len(callers), "callers": ";".join(callers),
            "n_source_samples": len(samples), "source_samples": ";".join(samples),
            "sources": ";".join(sorted({f"{m['tumor_id']}_vs_{m['normal_id']}:{m['caller']}:{m['record_id']}"
                                        for m in ms})),
            "called_in_target": sample_id in samples,
            "spread_1": max(m["pos1"] for m in ms) - min(m["pos1"] for m in ms),
            "spread_2": max(m["pos2"] for m in ms) - min(m["pos2"] for m in ms),
            "rep_caller": rep["caller"], "rep_source": f"{rep['tumor_id']}_vs_{rep['normal_id']}:{rep['record_id']}",
        })
        for m in sorted(ms, key=lambda m: (m["caller"], m["tumor_id"], str(m["record_id"]))):
            member_rows.append({"sv_id": sv_id, **{c: m[c] for c in MEMBER_COLUMNS}})
    union = pd.DataFrame(out, columns=UNION_COLUMNS)
    members = pd.DataFrame(member_rows, columns=["sv_id"] + MEMBER_COLUMNS)
    return union, members


def build_union(tsvs, sample_id, order, slop, rep_priority):
    rows = []
    for t in tsvs:
        df = pd.read_csv(t, sep="\t", dtype={"chrom1": str, "chrom2": str, "tumor_id": str, "normal_id": str,
                                             "record_id": str, "filter": str, "qual": str}, keep_default_na=False)
        rows.extend(df.to_dict("records"))
    return cluster_union(rows, sample_id, order, slop, rep_priority)


def _vcf_records(union, fasta=None):
    fa = pysam.FastaFile(fasta) if fasta else None

    def ref(chrom, pos):
        if fa is None or pos < 1:
            return "N"
        try:
            return (fa.fetch(chrom, pos - 1, pos) or "N").upper()
        except (KeyError, ValueError):
            return "N"

    recs = []
    for u in union.itertuples(index=False):
        imp = ";IMPRECISE" if u.imprecise else ""
        ci1, ci2 = f"{u.ci1_lo},{u.ci1_hi}", f"{u.ci2_lo},{u.ci2_hi}"
        if u.svtyper_type in ("DEL", "DUP"):
            pos, end = (u.position_1, u.position_2 - 1) if u.svtyper_type == "DEL" else (u.position_1 - 1, u.position_2)
            info = f"SVTYPE={u.svtyper_type};END={end};CIPOS={ci1};CIEND={ci2};EVENT={u.sv_id}{imp}"
            recs.append((u.chromosome_1, pos, u.sv_id, ref(u.chromosome_1, pos), f"<{u.svtyper_type}>", info))
            continue
        p1 = u.position_1 if u.strand_1 == "+" else u.position_1 - 1
        p2 = u.position_2 if u.strand_2 == "+" else u.position_2 - 1
        r1, r2 = ref(u.chromosome_1, p1), ref(u.chromosome_2, p2)
        id1, id2 = f"{u.sv_id}_1", f"{u.sv_id}_2"
        recs.append((u.chromosome_1, p1, id1, r1, gu.format_bnd_alt(r1, u.strand_1, u.chromosome_2, p2, u.strand_2),
                     f"SVTYPE=BND;MATEID={id2};CIPOS={ci1};EVENT={u.sv_id}{imp}"))
        recs.append((u.chromosome_2, p2, id2, r2, gu.format_bnd_alt(r2, u.strand_2, u.chromosome_1, p1, u.strand_1),
                     f"SVTYPE=BND;MATEID={id1};CIPOS={ci2};EVENT={u.sv_id}{imp}"))
    return recs


def write_svtyper_vcf(union, path, order, fai=None, fasta=None):
    if fasta and not fai:
        fai = f"{fasta}.fai"
    lengths = {}
    if fai:
        with open(fai) as fh:
            lengths = {l.split("\t")[0]: l.split("\t")[1] for l in fh if l.strip()}
    recs = _vcf_records(union, fasta)
    key = gu.site_sort_key(order)
    recs.sort(key=lambda r: (key(r[0], r[1]), r[2]))
    contigs = list(order) + sorted({r[0] for r in recs} - set(order))
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##source=nf-wgs-genotyping sv_union\n")
        for c in contigs:
            fh.write(f"##contig=<ID={c}" + (f",length={lengths[c]}" if c in lengths else "") + ">\n")
        fh.write('##ALT=<ID=DEL,Description="Deletion">\n##ALT=<ID=DUP,Description="Duplication">\n'
                 '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
                 '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n'
                 '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">\n'
                 '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">\n'
                 '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">\n'
                 '##INFO=<ID=EVENT,Number=1,Type=String,Description="Union SV id">\n'
                 '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise representative call">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, rid, ref, alt, info in recs:
            fh.write(f"{chrom}\t{pos}\t{rid}\t{ref}\t{alt}\t.\tPASS\t{info}\n")
    return path


def write_outputs(union, members, sample_id, outdir, chunk_size, order, fai=None, fasta=None):
    os.makedirs(os.path.join(outdir, "chunks"), exist_ok=True)
    union.to_csv(os.path.join(outdir, f"{sample_id}_union_svs.bedpe.gz"), sep="\t", index=False,
                 compression="gzip")
    members.to_csv(os.path.join(outdir, f"{sample_id}_union_svs.members.tsv.gz"), sep="\t", index=False,
                   compression="gzip")
    full = write_svtyper_vcf(union, os.path.join(outdir, f"{sample_id}_union_svs.vcf"), order, fai, fasta)
    pysam.tabix_index(full, preset="vcf", force=True)
    chunks = []
    n_chunks = max(1, -(-len(union) // chunk_size))
    for i in range(n_chunks):
        part = union.iloc[i * chunk_size:(i + 1) * chunk_size]
        cid = f"chunk_{i:05d}"
        write_svtyper_vcf(part, os.path.join(outdir, "chunks", f"{cid}.vcf"), order, fai, fasta)
        chunks.append((cid, len(part)))
    gu.write_manifest(os.path.join(outdir, "chunks.manifest"), chunks)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", nargs="+", required=True, help="sv_normalize TSVs")
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--reference", default=None, help="indexed FASTA (contig order, lengths, REF bases)")
    ap.add_argument("--contigs", default=None, help="contig order when no reference is given")
    ap.add_argument("--slop", type=int, default=100)
    ap.add_argument("--rep-priority", default="gridss,brass,svaba")
    ap.add_argument("--chunk-size", type=int, default=500)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)
    if a.chunk_size < 1 or a.slop < 0:
        ap.error("--chunk-size must be >= 1 and --slop >= 0")
    fai = f"{a.reference}.fai" if a.reference else None
    order = gu.contig_order(fai=fai, contigs=gu.parse_contigs(a.contigs))
    union, members = build_union(a.inputs, a.sample_id, order, a.slop, gu.parse_contigs(a.rep_priority) or [])
    write_outputs(union, members, a.sample_id, a.outdir, a.chunk_size, order, fai, a.reference)
    return 0


if __name__ == "__main__":
    sys.exit(main())
