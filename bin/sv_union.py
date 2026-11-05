#!/usr/bin/env python3
"""Union of normalized SV junctions across callers and source samples: every distinct adjacency is
a genotyping candidate; clusters of nearby adjacencies summarize them. Writes the cluster union
(BEDPE-like, with provenance), the adjacency table, the members table, an svtyper-compatible VCF,
and chunk files for scatter/gather (svtyper VCFs; insertion TSVs for sv_ins_count.py).

Adjacency: a distinct normalized junction, i.e. identical (chrom, pos, strand) at both ends, plus,
for insertions, the inserted sequence (or 'len:<n>' when only the length is known). Calls of the
same adjacency by several callers/samples are its members; its CI and IMPRECISE flag are those of
its best member (representative order below). Each adjacency is genotyped on its own, so nearby
but different breakpoints in a cluster all get a genotype (adj_id = <sv_id>.<k>).

Clusters: adjacencies are bucketed by (chrom1, chrom2, strand1, strand2, insertion or not). Two
adjacencies are linked when |pos1 - pos1'| <= slop and |pos2 - pos2'| <= slop. Links are taken in
increasing distance max(|dpos1|, |dpos2|) (ties by the adjacencies' coordinates, so the result does
not depend on input order) and merge two clusters only if every member of the merged cluster stays
within slop of every other at both ends (spread_1, spread_2 <= slop): bounded diameter, no
transitive chaining of distant events. n_linked_clusters counts other clusters linked to a cluster
(neighborhoods to review). Representative member: precise before imprecise, sequence-resolved
before length-only insertions, then --rep-priority caller order, narrowest CI, (tumor_id,
record_id). Cluster ids SV000001... follow genomic order of the representative; the cluster table
reports the representative's coordinates (rep_adj_id) as before.

svtyper VCF (non-insertion adjacencies; svtyper skips SVTYPE=INS): intrachromosomal (+,-) ->
<DEL> (POS=pos1, END=pos2-1); intrachromosomal (-,+) -> <DUP> (POS=pos1-1, END=pos2); everything
else -> a BND pair <adj_id>_1/<adj_id>_2 with MATEID. svtyper adds 1 to the POS of reverse ('-')
ends, so '-' ends are written one base left of the junction. Every record carries CIPOS (and
CIEND for DEL/DUP), which svtyper requires, EVENT=<adj_id> and CLUSTER=<sv_id>. Chunks never split
a BND pair. Insertions are written to ins_chunks/ (INS_SITE_COLUMNS) for sv_ins_count.py.
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
                 "called_in_target", "spread_1", "spread_2", "rep_caller", "rep_source",
                 "genotyper", "rep_adj_id", "n_adjacencies", "n_linked_clusters"]
ADJ_COLUMNS = ["adj_id", "sv_id", "chromosome_1", "position_1", "strand_1", "chromosome_2", "position_2",
               "strand_2", "ci1_lo", "ci1_hi", "ci2_lo", "ci2_hi", "imprecise", "sv_class", "genotyper",
               "svtyper_type", "sv_length", "ins_len", "ins_seq", "ins_seq_available", "n_members", "n_callers",
               "callers", "n_source_samples", "source_samples", "sources", "called_in_target",
               "is_cluster_rep"]
INS_SITE_COLUMNS = ["adj_id", "chromosome", "position", "del_len", "ins_len", "ins_seq", "ins_seq_available",
                    "imprecise"]


def svtyper_type(chrom1, strand1, chrom2, strand2):
    if chrom1 == chrom2 and (strand1, strand2) == ("+", "-"):
        return "DEL"
    if chrom1 == chrom2 and (strand1, strand2) == ("-", "+"):
        return "DUP"
    return "BND"


def _as_bool(v):
    return v if isinstance(v, bool) else str(v).strip().lower() in ("true", "1")


def _prep(rows):
    out = []
    for r in rows:
        r = dict(r)
        r["pos1"], r["pos2"] = int(r["pos1"]), int(r["pos2"])
        for k in ("ci1_lo", "ci1_hi", "ci2_lo", "ci2_hi"):
            r[k] = int(r[k])
        r["imprecise"] = _as_bool(r["imprecise"])
        r["ins_len"] = int(r.get("ins_len") or 0)
        seq = r.get("ins_seq")
        r["ins_seq"] = "" if seq is None or (isinstance(seq, float) and np.isnan(seq)) else str(seq)
        r["ins_seq_available"] = _as_bool(r.get("ins_seq_available", True))
        out.append(r)
    return out


def _is_ins(r):
    return r["sv_class"] == "INS"


def _adj_key(r):
    ins = (r["ins_seq"] if r["ins_seq_available"] else f"len:{r['ins_len']}") if _is_ins(r) else ""
    return (r["chrom1"], r["pos1"], r["strand1"], r["chrom2"], r["pos2"], r["strand2"], ins)


def _sv_length(r):
    if _is_ins(r):
        return float(r["ins_len"] - (r["pos2"] - r["pos1"] - 1))
    return float(r["pos2"] - r["pos1"]) if r["chrom1"] == r["chrom2"] else np.nan


def _provenance(ms, sample_id):
    callers = sorted({m["caller"] for m in ms})
    samples = sorted({m["tumor_id"] for m in ms})
    return {"n_members": len(ms), "n_callers": len(callers), "callers": ";".join(callers),
            "n_source_samples": len(samples), "source_samples": ";".join(samples),
            "sources": ";".join(sorted({f"{m['tumor_id']}_vs_{m['normal_id']}:{m['caller']}:{m['record_id']}"
                                        for m in ms})),
            "called_in_target": sample_id in samples}


def _clusters(keys, slop):
    """Bounded-diameter clustering of adjacency keys -> (root per key, linked root pairs)."""
    buckets = defaultdict(list)
    for k in keys:
        buckets[(k[0], k[3], k[2], k[5], k[6] != "")].append(k)
    links = []
    for ks in buckets.values():
        ks.sort(key=lambda k: (k[1], k[4], k))
        for i, a in enumerate(ks):
            for b in ks[i + 1:]:
                if b[1] - a[1] > slop:
                    break
                d2 = abs(b[4] - a[4])
                if d2 <= slop:
                    links.append((max(b[1] - a[1], d2), min(a, b), max(a, b)))
    links.sort()
    parent = {k: k for k in keys}
    bounds = {k: [k[1], k[1], k[4], k[4]] for k in keys}

    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    for _, a, b in links:
        ra, rb = find(a), find(b)
        if ra == rb:
            continue
        ba, bb = bounds[ra], bounds[rb]
        merged = [min(ba[0], bb[0]), max(ba[1], bb[1]), min(ba[2], bb[2]), max(ba[3], bb[3])]
        if merged[1] - merged[0] <= slop and merged[3] - merged[2] <= slop:
            root, other = min(ra, rb), max(ra, rb)
            parent[other] = root
            bounds[root] = merged
    roots = {k: find(k) for k in keys}
    linked = defaultdict(set)
    for _, a, b in links:
        if roots[a] != roots[b]:
            linked[roots[a]].add(roots[b])
            linked[roots[b]].add(roots[a])
    return roots, linked


def cluster_union(rows, sample_id, order, slop, rep_priority):
    """-> (cluster union, adjacencies, members) DataFrames (UNION_COLUMNS, ADJ_COLUMNS,
    ['sv_id', 'adj_id'] + MEMBER_COLUMNS)."""
    rows = _prep(rows)
    prio = {c: k for k, c in enumerate(rep_priority)}

    def rep_key(r):
        width = (r["ci1_hi"] - r["ci1_lo"]) + (r["ci2_hi"] - r["ci2_lo"])
        return (r["imprecise"], not r["ins_seq_available"], prio.get(r["caller"], len(prio)), width,
                r["tumor_id"], str(r["record_id"]), r["pos1"], r["pos2"], r["normal_id"])

    site = gu.site_sort_key(order)
    groups = defaultdict(list)
    for r in rows:
        groups[_adj_key(r)].append(r)
    keys = sorted(groups)
    roots, linked = _clusters(keys, slop)
    comps = defaultdict(list)
    for k in keys:
        comps[roots[k]].append(k)

    clusters = []
    for root, ks in comps.items():
        rep = min((m for k in ks for m in groups[k]), key=rep_key)
        clusters.append((rep, root, ks))
    clusters.sort(key=lambda c: (site(c[0]["chrom1"], c[0]["pos1"]), site(c[0]["chrom2"], c[0]["pos2"]),
                                 c[0]["strand1"], c[0]["strand2"], _adj_key(c[0])))
    out, adj_rows, member_rows = [], [], []
    for n, (rep, root, ks) in enumerate(clusters, start=1):
        sv_id = f"SV{n:06d}"
        ks = sorted(ks, key=lambda k: (site(k[0], k[1]), site(k[3], k[4]), k[2], k[5], k[6]))
        rep_key_adj = _adj_key(rep)
        ms = [m for k in ks for m in groups[k]]
        ins = _is_ins(rep)
        rep_adj_id = None
        for i, k in enumerate(ks, start=1):
            adj_id = f"{sv_id}.{i}"
            best = min(groups[k], key=rep_key)
            is_rep = k == rep_key_adj
            rep_adj_id = adj_id if is_rep else rep_adj_id
            adj_rows.append({
                "adj_id": adj_id, "sv_id": sv_id, "chromosome_1": k[0], "position_1": k[1], "strand_1": k[2],
                "chromosome_2": k[3], "position_2": k[4], "strand_2": k[5],
                "ci1_lo": best["ci1_lo"], "ci1_hi": best["ci1_hi"], "ci2_lo": best["ci2_lo"],
                "ci2_hi": best["ci2_hi"], "imprecise": best["imprecise"], "sv_class": best["sv_class"],
                "genotyper": "realign" if ins else "svtyper",
                "svtyper_type": "" if ins else svtyper_type(k[0], k[2], k[3], k[5]),
                "sv_length": _sv_length(best), "ins_len": best["ins_len"], "ins_seq": best["ins_seq"],
                "ins_seq_available": best["ins_seq_available"], **_provenance(groups[k], sample_id),
                "is_cluster_rep": is_rep})
            for m in sorted(groups[k], key=lambda m: (m["caller"], m["tumor_id"], str(m["record_id"]))):
                member_rows.append({"sv_id": sv_id, "adj_id": adj_id, **{c: m.get(c) for c in MEMBER_COLUMNS}})
        out.append({
            "sv_id": sv_id, "chromosome_1": rep["chrom1"], "position_1": rep["pos1"], "strand_1": rep["strand1"],
            "chromosome_2": rep["chrom2"], "position_2": rep["pos2"], "strand_2": rep["strand2"],
            "ci1_lo": rep["ci1_lo"], "ci1_hi": rep["ci1_hi"], "ci2_lo": rep["ci2_lo"], "ci2_hi": rep["ci2_hi"],
            "imprecise": rep["imprecise"], "sv_class": rep["sv_class"],
            "svtyper_type": "" if ins else svtyper_type(rep["chrom1"], rep["strand1"], rep["chrom2"], rep["strand2"]),
            "sv_length": _sv_length(rep), **_provenance(ms, sample_id),
            "spread_1": max(m["pos1"] for m in ms) - min(m["pos1"] for m in ms),
            "spread_2": max(m["pos2"] for m in ms) - min(m["pos2"] for m in ms),
            "rep_caller": rep["caller"], "rep_source": f"{rep['tumor_id']}_vs_{rep['normal_id']}:{rep['record_id']}",
            "genotyper": "realign" if ins else "svtyper", "rep_adj_id": rep_adj_id, "n_adjacencies": len(ks),
            "n_linked_clusters": len(linked.get(root, ())),
        })
    union = pd.DataFrame(out, columns=UNION_COLUMNS)
    adjs = pd.DataFrame(adj_rows, columns=ADJ_COLUMNS)
    members = pd.DataFrame(member_rows, columns=["sv_id", "adj_id"] + MEMBER_COLUMNS)
    return union, adjs, members


def build_union(tsvs, sample_id, order, slop, rep_priority):
    rows = []
    for t in tsvs:
        df = pd.read_csv(t, sep="\t", dtype={"chrom1": str, "chrom2": str, "tumor_id": str, "normal_id": str,
                                             "record_id": str, "filter": str, "qual": str, "ins_seq": str},
                         keep_default_na=False)
        rows.extend(df.to_dict("records"))
    return cluster_union(rows, sample_id, order, slop, rep_priority)


def _vcf_records(adjs, fasta=None):
    fa = pysam.FastaFile(fasta) if fasta else None

    def ref(chrom, pos):
        if fa is None or pos < 1:
            return "N"
        try:
            return (fa.fetch(chrom, pos - 1, pos) or "N").upper()
        except (KeyError, ValueError):
            return "N"

    recs = []
    for u in adjs.itertuples(index=False):
        if u.genotyper != "svtyper":
            continue
        imp = ";IMPRECISE" if _as_bool(u.imprecise) else ""
        ci1, ci2 = f"{u.ci1_lo},{u.ci1_hi}", f"{u.ci2_lo},{u.ci2_hi}"
        tags = f"EVENT={u.adj_id};CLUSTER={u.sv_id}{imp}"
        if u.svtyper_type in ("DEL", "DUP"):
            pos, end = (u.position_1, u.position_2 - 1) if u.svtyper_type == "DEL" else (u.position_1 - 1, u.position_2)
            info = f"SVTYPE={u.svtyper_type};END={end};CIPOS={ci1};CIEND={ci2};{tags}"
            recs.append((u.chromosome_1, pos, u.adj_id, ref(u.chromosome_1, pos), f"<{u.svtyper_type}>", info))
            continue
        p1 = u.position_1 if u.strand_1 == "+" else u.position_1 - 1
        p2 = u.position_2 if u.strand_2 == "+" else u.position_2 - 1
        r1, r2 = ref(u.chromosome_1, p1), ref(u.chromosome_2, p2)
        id1, id2 = f"{u.adj_id}_1", f"{u.adj_id}_2"
        recs.append((u.chromosome_1, p1, id1, r1, gu.format_bnd_alt(r1, u.strand_1, u.chromosome_2, p2, u.strand_2),
                     f"SVTYPE=BND;MATEID={id2};CIPOS={ci1};{tags}"))
        recs.append((u.chromosome_2, p2, id2, r2, gu.format_bnd_alt(r2, u.strand_2, u.chromosome_1, p1, u.strand_1),
                     f"SVTYPE=BND;MATEID={id1};CIPOS={ci2};{tags}"))
    if fa is not None:
        fa.close()
    return recs


def write_svtyper_vcf(adjs, path, order, fai=None, fasta=None):
    """svtyper input VCF of the svtyper-genotyped rows of an adjacency table."""
    if fasta and not fai:
        fai = f"{fasta}.fai"
    lengths = {}
    if fai:
        with open(fai) as fh:
            lengths = {l.split("\t")[0]: l.split("\t")[1] for l in fh if l.strip()}
    recs = _vcf_records(adjs, fasta)
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
                 '##INFO=<ID=EVENT,Number=1,Type=String,Description="Union adjacency id (adj_id)">\n'
                 '##INFO=<ID=CLUSTER,Number=1,Type=String,Description="Union SV cluster id (sv_id)">\n'
                 '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise call (best member of the '
                 'adjacency)">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, rid, ref, alt, info in recs:
            fh.write(f"{chrom}\t{pos}\t{rid}\t{ref}\t{alt}\t.\tPASS\t{info}\n")
    return path


def insertion_sites(adjs):
    ins = adjs[adjs["genotyper"] == "realign"]
    return pd.DataFrame({
        "adj_id": ins["adj_id"], "chromosome": ins["chromosome_1"], "position": ins["position_1"],
        "del_len": ins["position_2"] - ins["position_1"] - 1, "ins_len": ins["ins_len"], "ins_seq": ins["ins_seq"],
        "ins_seq_available": ins["ins_seq_available"], "imprecise": ins["imprecise"]}, columns=INS_SITE_COLUMNS)


def write_outputs(union, adjs, members, sample_id, outdir, chunk_size, order, fai=None, fasta=None,
                  ins_chunk_size=200):
    for sub in ("chunks", "ins_chunks"):
        os.makedirs(os.path.join(outdir, sub), exist_ok=True)
    union.to_csv(os.path.join(outdir, f"{sample_id}_union_svs.bedpe.gz"), sep="\t", index=False,
                 compression="gzip")
    adjs.to_csv(os.path.join(outdir, f"{sample_id}_union_sv_adjacencies.tsv.gz"), sep="\t", index=False,
                compression="gzip")
    members.to_csv(os.path.join(outdir, f"{sample_id}_union_svs.members.tsv.gz"), sep="\t", index=False,
                   compression="gzip")
    sv = adjs[adjs["genotyper"] == "svtyper"]
    full = write_svtyper_vcf(sv, os.path.join(outdir, f"{sample_id}_union_svs.vcf"), order, fai, fasta)
    pysam.tabix_index(full, preset="vcf", force=True)
    chunks = []
    for i in range(max(1, -(-len(sv) // chunk_size))):
        part = sv.iloc[i * chunk_size:(i + 1) * chunk_size]
        cid = f"chunk_{i:05d}"
        write_svtyper_vcf(part, os.path.join(outdir, "chunks", f"{cid}.vcf"), order, fai, fasta)
        chunks.append((cid, len(part)))
    gu.write_manifest(os.path.join(outdir, "chunks.manifest"), chunks)
    sites = insertion_sites(adjs)
    ins_chunks = []
    for i in range(max(1, -(-len(sites) // ins_chunk_size))):
        part = sites.iloc[i * ins_chunk_size:(i + 1) * ins_chunk_size]
        cid = f"ins_chunk_{i:05d}"
        part.to_csv(os.path.join(outdir, "ins_chunks", f"{cid}.tsv"), sep="\t", index=False)
        ins_chunks.append((cid, len(part)))
    gu.write_manifest(os.path.join(outdir, "ins_chunks.manifest"), ins_chunks)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", nargs="+", required=True, help="sv_normalize TSVs")
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--reference", default=None, help="indexed FASTA (contig order, lengths, REF bases)")
    ap.add_argument("--contigs", default=None, help="contig order when no reference is given")
    ap.add_argument("--slop", type=int, default=100)
    ap.add_argument("--rep-priority", default="gridss,brass,svaba")
    ap.add_argument("--chunk-size", type=int, default=500, help="adjacencies per svtyper chunk")
    ap.add_argument("--ins-chunk-size", type=int, default=200, help="insertions per sv_ins_count chunk")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)
    if a.chunk_size < 1 or a.ins_chunk_size < 1 or a.slop < 0:
        ap.error("--chunk-size and --ins-chunk-size must be >= 1 and --slop >= 0")
    fai = f"{a.reference}.fai" if a.reference else None
    order = gu.contig_order(fai=fai, contigs=gu.parse_contigs(a.contigs))
    union, adjs, members = build_union(a.inputs, a.sample_id, order, a.slop,
                                       gu.parse_contigs(a.rep_priority) or [])
    write_outputs(union, adjs, members, a.sample_id, a.outdir, a.chunk_size, order, fai, a.reference,
                  a.ins_chunk_size)
    return 0


if __name__ == "__main__":
    sys.exit(main())
