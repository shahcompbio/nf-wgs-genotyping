import glob
import json
import math
import os
import shutil
import subprocess

import pandas as pd
import pysam
import pytest
import yaml

import snv_count
import snv_extract
import snv_finalize
import snv_union

FINAL_COLUMNS = ["sample_id", "chromosome", "position", "ref", "alt", "ref_count", "alt_count",
                 "other_count", "depth", "vaf", "n_callers", "callers", "n_source_samples",
                 "source_samples", "sources", "called_in_target"]


def _run_upstream(synth, tmp_path, chunk_size=3):
    stats = []
    tsvs = []
    for caller, tumor in (("caveman", "T1"), ("mutect", "T1"), ("strelka", "T2")):
        out = tmp_path / f"{tumor}.{caller}.tsv"
        st = tmp_path / f"{tumor}.{caller}.json"
        snv_extract.main(["--vcf", synth["snv_vcfs"][caller], "--caller", caller, "--tumor-id", tumor,
                          "--normal-id", "N1", "--contigs", "1,2", "--reference", synth["fasta"],
                          "--out-tsv", str(out), "--out-stats", str(st)])
        tsvs.append(str(out))
        stats.append(str(st))
    udir = tmp_path / "union"
    udir.mkdir()
    snv_union.main(["--inputs", *tsvs, "--sample-id", "T1", "--reference", synth["fasta"],
                    "--chunk-size", str(chunk_size), "--outdir", str(udir)])
    counts = []
    for chunk in sorted(glob.glob(str(udir / "chunks" / "chunk_*.tsv"))):
        cid = os.path.basename(chunk)[:-4]
        out = tmp_path / f"counts_{cid}.tsv"
        snv_count.main(["--bam", synth["bam"], "--sites", chunk, "--out", str(out)])
        counts.append(str(out))
    return udir, counts, stats


def test_finalize_joins_counts_onto_union(synth, tmp_path):
    udir, counts, stats = _run_upstream(synth, tmp_path)
    df, qc = snv_finalize.finalize(str(udir / "T1_union_snvs.tsv.gz"), counts, str(udir / "chunks.manifest"),
                                   "T1", stats)
    assert list(df.columns) == FINAL_COLUMNS
    assert len(df) == 4
    s1 = df.iloc[0]
    t = synth["truth_S1"]
    assert (s1["ref_count"], s1["alt_count"], s1["other_count"]) == (t["ref_count"], t["alt_count"], t["other_count"])
    assert s1["depth"] == t["ref_count"] + t["alt_count"] + t["other_count"]
    assert s1["vaf"] == pytest.approx(t["alt_count"] / (t["ref_count"] + t["alt_count"]))
    assert math.isnan(df.iloc[3]["vaf"])
    assert (df["sample_id"] == "T1").all()
    assert qc["union_size"] == 4
    assert qc["called_in_target"] == 2
    assert qc["positive_control_fraction"] == pytest.approx(0.5)
    assert qc["conflict_fragments"] == t["conflict_count"]
    assert qc["inputs"][0]["kept"] >= 1


def test_finalize_fails_on_missing_chunk(synth, tmp_path):
    udir, counts, stats = _run_upstream(synth, tmp_path)
    with pytest.raises(ValueError, match="missing"):
        snv_finalize.finalize(str(udir / "T1_union_snvs.tsv.gz"), counts[:-1], str(udir / "chunks.manifest"),
                              "T1", stats)


def test_finalize_fails_when_a_site_has_no_counts_row(synth, tmp_path):
    udir, counts, stats = _run_upstream(synth, tmp_path)
    df = pd.read_csv(counts[0], sep="\t", dtype={"chromosome": str})
    df.iloc[:-1].to_csv(counts[0], sep="\t", index=False)
    with pytest.raises(ValueError, match="counts rows"):
        snv_finalize.finalize(str(udir / "T1_union_snvs.tsv.gz"), counts, str(udir / "chunks.manifest"),
                              "T1", stats)


def test_cli_writes_table_yaml_and_qc(synth, tmp_path):
    udir, counts, stats = _run_upstream(synth, tmp_path)
    out = tmp_path / "final"
    out.mkdir()
    _finalize_cli(udir, counts, stats, out)
    table = pd.read_csv(out / "T1_snv_genotyping.csv.gz", dtype={"chromosome": str})
    assert list(table.columns) == FINAL_COLUMNS
    meta = yaml.safe_load(open(out / "T1_snv_genotyping.csv.gz.yaml"))
    assert [c["name"] for c in meta["columns"]] == FINAL_COLUMNS
    qc = json.load(open(out / "T1_genotyping_qc.json"))
    assert qc["sample_id"] == "T1" and qc["type"] == "snv_genotyping"


def _finalize_cli(udir, counts, stats, out, sample="T1"):
    snv_finalize.main(["--union", str(udir / f"{sample}_union_snvs.tsv.gz"),
                       "--union-vcf", str(udir / f"{sample}_union_snvs.vcf.gz"), "--counts", *counts,
                       "--manifest", str(udir / "chunks.manifest"), "--stats", *stats,
                       "--sample-id", sample, "--outdir", str(out)])


def test_cli_writes_vcf_matching_table(synth, tmp_path):
    udir, counts, stats = _run_upstream(synth, tmp_path)
    out = tmp_path / "final"
    out.mkdir()
    _finalize_cli(udir, counts, stats, out)
    table = pd.read_csv(out / "T1_snv_genotyping.csv.gz", dtype={"chromosome": str})
    path = out / "T1_snv_genotyping.vcf.gz"
    assert (out / "T1_snv_genotyping.vcf.gz.tbi").exists()
    fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(synth["fasta"] + ".fai") if l.strip()}
    with pysam.VariantFile(str(path)) as vf:
        assert tuple(vf.header.samples) == ("T1",)
        assert {c: vf.header.contigs[c].length for c in ("1", "2")} == {c: fai[c] for c in ("1", "2")}
        records = list(vf)
    keys = [(r.chrom, r.pos, r.ref, r.alts[0]) for r in records]
    assert keys == list(table[["chromosome", "position", "ref", "alt"]].itertuples(index=False, name=None))
    assert any(math.isnan(v) for v in table["vaf"])   # the missing-AF branch is exercised
    for rec, row in zip(records, table.itertuples(index=False)):
        s = rec.samples["T1"]
        assert "GT" not in rec.format                  # counts only: the pipeline calls no genotype
        assert s["AD"] == (row.ref_count, row.alt_count)
        assert (s["DP"], s["OC"]) == (row.depth, row.other_count)
        if math.isnan(row.vaf):
            assert s.get("AF") in (None, (None,))
        else:
            assert s["AF"][0] == pytest.approx(row.vaf, rel=1e-5)
        assert list(rec.filter) == ["PASS"]
        assert rec.info["CALLERS"] == tuple(row.callers.split(";"))
        assert rec.info["SAMPLES"] == tuple(row.source_samples.split(";"))
        assert rec.info["SOURCES"] == tuple(row.sources.split(";"))
        assert ("INTARGET" in rec.info) == bool(row.called_in_target)
    with pysam.VariantFile(str(path)) as vf:          # tabix index usable
        assert len(list(vf.fetch("1"))) == int((table["chromosome"] == "1").sum())
    if shutil.which("bcftools"):
        res = subprocess.run(["bcftools", "view", "-H", str(path)], capture_output=True, text=True)
        assert res.returncode == 0, res.stderr
        assert len(res.stdout.splitlines()) == len(table)


def test_empty_union_end_to_end_gives_header_only_table(synth, tmp_path):
    vcf = tmp_path / "none.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    tsv, st = tmp_path / "x.tsv", tmp_path / "x.json"
    snv_extract.main(["--vcf", str(vcf), "--caller", "mutect", "--tumor-id", "T1", "--normal-id", "N1",
                      "--out-tsv", str(tsv), "--out-stats", str(st)])
    udir = tmp_path / "union"
    udir.mkdir()
    snv_union.main(["--inputs", str(tsv), "--sample-id", "T1", "--outdir", str(udir)])
    counts = tmp_path / "counts_chunk_00000.tsv"
    snv_count.main(["--bam", synth["bam"], "--sites", str(udir / "chunks" / "chunk_00000.tsv"), "--out", str(counts)])
    out = tmp_path / "final"
    out.mkdir()
    _finalize_cli(udir, [str(counts)], [str(st)], out)
    table = pd.read_csv(out / "T1_snv_genotyping.csv.gz")
    assert len(table) == 0 and list(table.columns) == FINAL_COLUMNS
    assert json.load(open(out / "T1_genotyping_qc.json"))["positive_control_fraction"] is None
    with pysam.VariantFile(str(out / "T1_snv_genotyping.vcf.gz")) as vf:
        assert tuple(vf.header.samples) == ("T1",) and len(list(vf)) == 0
