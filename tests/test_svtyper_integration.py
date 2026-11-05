"""Integration tests against the real svtyper fork (skipped if it cannot be imported)."""
import os
import subprocess

import pandas as pd
import pytest

import svtyper_helpers as sh
import geno_utils as gu
import sv_finalize
import sv_ins_count
import sv_normalize
import sv_union
import svtyper_run

pytestmark = pytest.mark.skipif(not sh.available(), reason="svtyper fork not importable (set SVTYPER_SRC)")

HERE = os.path.dirname(os.path.abspath(__file__))
SVD = os.path.join(HERE, "data", "svtyper")
NA_BAM = os.path.join(SVD, "NA12878.target_loci.sorted.bam")
NA_LIB = os.path.join(SVD, "NA12878.bam.json")
CLIP = os.path.join(SVD, "clipmatch")


def _svtyper(vcf, bam, lib, out, reference=None, extra=""):
    cmd, _ = svtyper_run.build_command(sh.svtyper_cmd(), vcf, bam, lib, out, reference=reference,
                                       extra_args=extra)
    subprocess.run(cmd, env=sh.svtyper_env(), check=True, capture_output=True)
    return out


def _union_from_vcf(vcf, tmp_path, sample, order, contigs=None, slop=0, fasta=None, chunk_size=10000):
    rows, _ = sv_normalize.normalize(vcf, "generic", sample, "N", contigs=contigs, order=order, min_size=1)
    union, adjs, members = sv_union.cluster_union(rows, sample, order, slop, ["generic"])
    od = tmp_path / "union"
    od.mkdir(exist_ok=True)
    sv_union.write_outputs(union, adjs, members, sample, str(od), chunk_size, order, fasta=fasta)
    return od


def _ins_counts(od, bam, tmp_path, fasta=None):
    """Run sv_ins_count.py on every insertion chunk of a union directory."""
    out = []
    for chunk in sorted((od / "ins_chunks").glob("ins_chunk_*.tsv")):
        p = tmp_path / f"ins_counts_{chunk.name}"
        sv_ins_count.main(["--bam", bam, "--sites", str(chunk), "--out", str(p),
                           *(["--reference", fasta] if fasta else [])])
        out.append(str(p))
    return out


def _finalize(od, sample, gt_files, ins_files):
    return sv_finalize.finalize(str(od / f"{sample}_union_svs.bedpe.gz"),
                                str(od / f"{sample}_union_sv_adjacencies.tsv.gz"), gt_files,
                                str(od / "chunks.manifest"), sample, [], ins_files, str(od / "ins_chunks.manifest"))


def test_writer_positions_match_svtyper_breakpoint_parser(tmp_path):
    parsers = sh.import_parsers()
    order = {"1": 0, "2": 1, "3": 2}
    rows = []
    for (c1, p1, s1, c2, p2, s2) in [("1", 100, "+", "1", 900, "-"), ("1", 2000, "-", "1", 2600, "+"),
                                     ("2", 100, "+", "2", 700, "+"), ("2", 1000, "-", "2", 4000, "-"),
                                     ("2", 5000, "-", "3", 400, "+"), ("2", 6000, "+", "3", 900, "-")]:
        rows.append({"chrom1": c1, "pos1": p1, "strand1": s1, "chrom2": c2, "pos2": p2, "strand2": s2,
                     "ci1_lo": 0, "ci1_hi": 0, "ci2_lo": 0, "ci2_hi": 0, "imprecise": False,
                     "sv_class": sv_normalize.sv_class(c1, s1, c2, s2), "caller": "generic", "tumor_id": "T",
                     "normal_id": "N", "record_id": f"r{len(rows)}", "filter": "PASS", "qual": ".", "ins_len": 0})
    union, adjs, _ = sv_union.cluster_union(rows, "T", order, 0, ["generic"])
    path = str(tmp_path / "u.vcf")
    sv_union.write_svtyper_vcf(adjs, path, order)
    lines = open(path).read().splitlines(keepends=True)
    vcf = parsers.Vcf()
    vcf.add_header([l for l in lines if l.startswith("#")])
    got = {}
    for line in lines:
        if line.startswith("#"):
            continue
        var = parsers.Variant(line.rstrip("\n").split("\t"), vcf)
        bp = vcf.get_variant_breakpoints(var, 1e10)
        if bp is not None:
            got[bp["id"].split("_")[0]] = ((bp["A"]["chrom"], bp["A"]["pos"], bp["A"]["is_reverse"]),
                                           (bp["B"]["chrom"], bp["B"]["pos"], bp["B"]["is_reverse"]))
    assert len(got) == len(adjs) == 6
    for u in adjs.itertuples(index=False):
        assert got[u.adj_id] == ((u.chromosome_1, u.position_1, u.strand_1 == "-"),
                                 (u.chromosome_2, u.position_2, u.strand_2 == "-")), u.adj_id


def _parse_counts(path, fields=("AO", "RO", "AS", "AP", "RS", "RP")):
    out = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        t = line.rstrip("\n").split("\t")
        fmt = dict(zip(t[8].split(":"), t[9].split(":")))
        end = [kv.split("=")[1] for kv in t[7].split(";") if kv.startswith("END=")]
        out[(t[0], int(t[1]), int(end[0]))] = tuple(fmt[f] for f in fields)
    return out


def test_del_dup_round_trip_matches_svtyper_on_original_records(tmp_path):
    lines = open(os.path.join(SVD, "example.vcf")).read().splitlines(keepends=True)
    head = [l for l in lines if l.startswith("#")]
    recs = [l for l in lines if not l.startswith("#")]
    dels = [l for l in recs if "SVTYPE=DEL" in l][:25]
    dups = [l for l in recs if "SVTYPE=DUP" in l]
    subset = tmp_path / "subset.vcf"
    subset.write_text("".join(head + dels + dups))
    orig = _parse_counts(_svtyper(str(subset), NA_BAM, NA_LIB, str(tmp_path / "orig.gt.vcf")))
    od = _union_from_vcf(str(subset), tmp_path, "NA12878", gu.contig_order(bam=NA_BAM))
    ours = _parse_counts(_svtyper(str(od / "chunks" / "chunk_00000.vcf"), NA_BAM, NA_LIB,
                                  str(tmp_path / "ours.gt.vcf")))
    assert len(orig) == len(dels) + len(dups)
    assert ours == orig


def test_clipmatch_fixture_through_pipeline_scripts(tmp_path):
    fasta = os.path.join(CLIP, "ref.fa")
    order = gu.contig_order(fai=fasta + ".fai")
    od = _union_from_vcf(os.path.join(CLIP, "svs.vcf"), tmp_path, "clip", order, slop=0, fasta=fasta)
    gt = _svtyper(str(od / "chunks" / "chunk_00000.vcf"), os.path.join(CLIP, "reads.bam"),
                  os.path.join(CLIP, "lib.json"), str(tmp_path / "chunk_00000.gt.vcf"), reference=fasta,
                  extra=sh.DEFAULT_SVTYPER_ARGS + " -m 25")
    df, _, _ = _finalize(od, "clip", [gt], _ins_counts(od, os.path.join(CLIP, "reads.bam"), tmp_path, fasta))
    by_ctg = df.set_index("chromosome_1")
    sv1 = by_ctg.loc["ctg1"]
    assert (sv1["alt_fragments"], sv1["alt_split_fragments"], sv1["alt_span_fragments"],
            sv1["alt_clip_fragments"]) == (19, 8, 7, 4)
    assert by_ctg.loc["ctg3", "alt_fragments"] == 0
    assert by_ctg.loc["ctg5", "alt_fragments"] == 0


def test_synthetic_sv_fixture_end_to_end_on_na12878(synth, tmp_path):
    contigs = [str(i) for i in range(1, 23)]
    order = gu.contig_order(bam=NA_BAM)
    tsvs = []
    for (who, caller), vcf in sorted(synth["sv_vcfs"].items()):
        tumor = synth["sv_target"] if who == "T" else synth["other_tumor"]
        rows, _ = sv_normalize.normalize(vcf, caller, tumor, "N1", contigs=contigs, order=order)
        p = tmp_path / f"{who}.{caller}.tsv"
        pd.DataFrame(rows, columns=sv_normalize.COLUMNS).to_csv(p, sep="\t", index=False)
        tsvs.append(str(p))
    union, adjs, members = sv_union.build_union(tsvs, synth["sv_target"], order, 100, ["gridss", "brass", "svaba"])
    od = tmp_path / "u"
    od.mkdir()
    sv_union.write_outputs(union, adjs, members, synth["sv_target"], str(od), 500, order)
    lib = str(tmp_path / "lib.json")
    cmd, _ = svtyper_run.build_command(sh.svtyper_cmd(), None, NA_BAM, lib, None, lib_only=True)
    subprocess.run(cmd, env=sh.svtyper_env(), check=True, capture_output=True)
    assert os.path.getsize(lib) > 0
    gt = _svtyper(str(od / "chunks" / "chunk_00000.vcf"), NA_BAM, lib, str(tmp_path / "chunk_00000.gt.vcf"),
                  extra=sh.DEFAULT_SVTYPER_ARGS)
    df, adj, qc = _finalize(od, synth["sv_target"], [gt], _ins_counts(od, NA_BAM, tmp_path))
    assert len(df) == synth["sv_expected"]["n_union"]
    assert qc["positive_control_fraction"] == 1.0
    assert (df[~df["called_in_target"]]["alt_fragments"] == 0).all()
    assert (df["genotype_status"] != "missing_from_output").all()
    # every adjacency has its own svtyper record and genotype: the homozygous deletion 72354 was
    # called by BRASS at the true breakpoints and by SvABA 15 bp to the right
    assert len(adj) == len(df) + 2 and (adj["genotype_status"] != "missing_from_output").all()
    pos, end = synth["sv_expected"]["shifted_del"]
    pair = adj[(adj["chromosome_1"] == "2") & adj["position_1"].isin([pos, pos + 15])].sort_values("position_1")
    assert pair["sv_id"].nunique() == 1 and len(pair) == 2 and pair["is_cluster_rep"].tolist() == [True, False]
    assert (pair["alt_fragments"] >= 1).all()
    assert pair.iloc[0]["callers"] == "brass" and pair.iloc[1]["callers"] == "svaba"
    with open(gt) as fh:
        ids = {l.split("\t")[2] for l in fh if not l.startswith("#")}
    assert set(adj["adj_id"]) == {i[:-2] if i.endswith(("_1", "_2")) else i for i in ids}
