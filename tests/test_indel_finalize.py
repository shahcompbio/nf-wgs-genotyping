import glob
import json
import os
import shutil
import subprocess

import pandas as pd
import pysam
import pytest
import yaml

import indel_count
import indel_extract
import indel_finalize
import indel_union
import presence_gt

FINAL_COLUMNS = ["sample_id", "chromosome", "position", "ref", "alt", "indel_length", "ref_count",
                 "alt_count", "other_count", "depth", "vaf", "uninformative_count", "conflict_count",
                 "locus_id", "var_start", "var_end", "n_callers", "callers", "n_source_samples",
                 "source_samples", "sources", "called_in_target", "gt", "gq", "pl"]
NAMES = ["D1", "HP1", "HP2", "I3", "C4"]


def _run_upstream(synth, tmp_path, chunk_size=2):
    tsvs, stats = [], []
    for caller, tumor in (("pindel", "T1"), ("strelka", "T1"), ("mutect", "T2")):
        out, st = tmp_path / f"{tumor}.{caller}.indels.tsv", tmp_path / f"{tumor}.{caller}.indels.stats.json"
        indel_extract.main(["--vcf", synth["indel_vcfs"][caller], "--caller", caller, "--tumor-id", tumor,
                            "--normal-id", "N1", "--contigs", "1,2,3", "--reference", synth["fasta"],
                            "--out-tsv", str(out), "--out-stats", str(st)])
        tsvs.append(str(out))
        stats.append(str(st))
    udir = tmp_path / "union"
    indel_union.main(["--inputs", *tsvs, "--sample-id", "T1", "--reference", synth["fasta"],
                      "--chunk-size", str(chunk_size), "--outdir", str(udir)])
    counts = []
    for chunk in sorted(glob.glob(str(udir / "chunks" / "chunk_*.tsv"))):
        out = tmp_path / f"counts_{os.path.basename(chunk)}"
        indel_count.main(["--bam", synth["bam"], "--reference", synth["fasta"], "--sites", chunk,
                          "--out", str(out)])
        counts.append(str(out))
    return udir, counts, stats


def test_finalize_joins_counts_onto_union(synth, tmp_path):
    udir, counts, stats = _run_upstream(synth, tmp_path)
    df, qc = indel_finalize.finalize(str(udir / "T1_union_indels.tsv.gz"), counts,
                                     str(udir / "chunks.manifest"), "T1", stats)
    assert list(df.columns) == FINAL_COLUMNS
    assert len(df) == 5
    for name, row in zip(NAMES, df.itertuples(index=False)):
        t = synth["truth_indels"][name]
        assert {k: getattr(row, k) for k in t} == t, name
        assert row.depth == t["ref_count"] + t["alt_count"] + t["other_count"]
        assert row.vaf == pytest.approx(t["alt_count"] / (t["ref_count"] + t["alt_count"]))
    assert (df["sample_id"] == "T1").all()
    assert qc["type"] == "indel_genotyping" and qc["union_size"] == 5 and qc["loci"] == 4
    assert qc["called_in_target"] == 4 and qc["positive_control_fraction"] == pytest.approx(1.0)
    # per-locus counts are summed once per locus, not once per allele
    assert qc["uninformative_fragments"] == 1 + 2 + 1 + 0
    assert qc["conflict_fragments"] == 1
    assert qc["chunks"] == 3 and len(qc["inputs"]) == 3
    # every designed allele has ref and alt fragments: both alleles present (indel error rate 0.01)
    assert df["gt"].tolist() == ["0/1"] * 5 and (df["gq"] >= 20).all()
    gts, gqs, pls = presence_gt.call(df["ref_count"], df["alt_count"], 0.01)
    assert df["gq"].tolist() == gqs.tolist() and df["pl"].tolist() == [presence_gt.format_pl(p) for p in pls]
    assert qc["gt_counts"] == {"0/0": 0, "0/1": 5, "1/1": 0, "./.": 0}
    assert qc["gt_model"]["error_rate"] == 0.01 and "empirical_error_rate" not in qc


def test_finalize_fails_on_missing_chunk_or_missing_rows(synth, tmp_path):
    udir, counts, stats = _run_upstream(synth, tmp_path)
    args = (str(udir / "T1_union_indels.tsv.gz"),)
    with pytest.raises(ValueError, match="missing"):
        indel_finalize.finalize(*args, counts[:-1], str(udir / "chunks.manifest"), "T1", stats)
    df = pd.read_csv(counts[0], sep="\t", dtype={"chromosome": str})
    df.iloc[:-1].to_csv(counts[0], sep="\t", index=False)
    with pytest.raises(ValueError, match="counts rows"):
        indel_finalize.finalize(*args, counts, str(udir / "chunks.manifest"), "T1", stats)


def _finalize_cli(udir, counts, stats, out, sample="T1", extra=()):
    indel_finalize.main(["--union", str(udir / f"{sample}_union_indels.tsv.gz"),
                         "--union-vcf", str(udir / f"{sample}_union_indels.vcf.gz"), "--counts", *counts,
                         "--manifest", str(udir / "chunks.manifest"), "--stats", *stats,
                         "--sample-id", sample, "--outdir", str(out), *extra])


def test_cli_writes_table_yaml_qc_and_matching_vcf(synth, tmp_path):
    udir, counts, stats = _run_upstream(synth, tmp_path)
    out = tmp_path / "final"
    out.mkdir()
    _finalize_cli(udir, counts, stats, out)
    table = pd.read_csv(out / "T1_indel_genotyping.csv.gz", dtype={"chromosome": str})
    assert list(table.columns) == FINAL_COLUMNS
    meta = yaml.safe_load(open(out / "T1_indel_genotyping.csv.gz.yaml"))
    assert [c["name"] for c in meta["columns"]] == FINAL_COLUMNS
    assert json.load(open(out / "T1_indel_qc.json"))["sample_id"] == "T1"
    path = out / "T1_indel_genotyping.vcf.gz"
    assert (out / "T1_indel_genotyping.vcf.gz.tbi").exists()
    with pysam.VariantFile(str(path)) as vf:
        assert tuple(vf.header.samples) == ("T1",)
        records = list(vf)
    keys = [(r.chrom, r.pos, r.ref, r.alts[0]) for r in records]
    assert keys == list(table[["chromosome", "position", "ref", "alt"]].itertuples(index=False, name=None))
    for rec, row in zip(records, table.itertuples(index=False)):
        s = rec.samples["T1"]
        assert list(rec.format)[0] == "GT"
        assert s["GT"] == presence_gt.gt_tuple(row.gt) and s["GQ"] == row.gq
        assert s["PL"] == tuple(int(x) for x in row.pl.split(","))
        assert s["AD"] == (row.ref_count, row.alt_count)
        assert (s["DP"], s["OC"], s["UF"]) == (row.depth, row.other_count, row.uninformative_count)
        assert s["AF"][0] == pytest.approx(row.vaf, rel=1e-5)
        assert rec.info["LOCUS"] == row.locus_id
        assert rec.info["SOURCES"] == tuple(row.sources.split(";"))
        assert ("INTARGET" in rec.info) == bool(row.called_in_target)
    if shutil.which("bcftools"):
        res = subprocess.run(["bcftools", "view", "-H", str(path)], capture_output=True, text=True)
        assert res.returncode == 0, res.stderr
        assert len(res.stdout.splitlines()) == len(table)


def test_cli_indel_gt_parameters(synth, tmp_path):
    udir, counts, stats = _run_upstream(synth, tmp_path)
    out = tmp_path / "final"
    out.mkdir()
    _finalize_cli(udir, counts, stats, out, extra=["--gt-error-rate", "0.03", "--gt-min-gq", "40"])
    qc = json.load(open(out / "T1_indel_qc.json"))
    assert (qc["gt_model"]["error_rate"], qc["gt_model"]["min_gq"]) == (0.03, 40)
    table = pd.read_csv(out / "T1_indel_genotyping.csv.gz", dtype={"chromosome": str})
    assert "./." in set(table["gt"])          # the 3/3 complex allele is not called at GQ >= 40
    with pysam.VariantFile(str(out / "T1_indel_genotyping.vcf.gz")) as vf:
        assert "error_rate=0.03," in str(vf.header)


def test_empty_union_end_to_end(synth, tmp_path):
    vcf = tmp_path / "none.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    tsv, st = tmp_path / "x.tsv", tmp_path / "x.json"
    indel_extract.main(["--vcf", str(vcf), "--caller", "pindel", "--tumor-id", "T1", "--normal-id", "N1",
                        "--reference", synth["fasta"], "--out-tsv", str(tsv), "--out-stats", str(st)])
    udir = tmp_path / "union"
    indel_union.main(["--inputs", str(tsv), "--sample-id", "T1", "--reference", synth["fasta"],
                      "--outdir", str(udir)])
    counts = tmp_path / "counts_chunk_00000.tsv"
    indel_count.main(["--bam", synth["bam"], "--reference", synth["fasta"], "--sites",
                      str(udir / "chunks" / "chunk_00000.tsv"), "--out", str(counts)])
    out = tmp_path / "final"
    out.mkdir()
    _finalize_cli(udir, [str(counts)], [str(st)], out)
    table = pd.read_csv(out / "T1_indel_genotyping.csv.gz")
    assert len(table) == 0 and list(table.columns) == FINAL_COLUMNS
    qc = json.load(open(out / "T1_indel_qc.json"))
    assert qc["positive_control_fraction"] is None and qc["loci"] == 0
    with pysam.VariantFile(str(out / "T1_indel_genotyping.vcf.gz")) as vf:
        assert len(list(vf)) == 0
