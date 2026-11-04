import gzip
import os

import pandas as pd
import pysam

import geno_utils as gu
import snv_extract
import snv_union


def _extracted(synth, tmp_path):
    paths = []
    for caller, tumor in (("caveman", "T1"), ("mutect", "T1"), ("strelka", "T2")):
        out = tmp_path / f"{tumor}.{caller}.tsv"
        snv_extract.main(["--vcf", synth["snv_vcfs"][caller], "--caller", caller, "--tumor-id", tumor,
                          "--normal-id", "N1", "--contigs", "1,2", "--out-tsv", str(out),
                          "--out-stats", str(tmp_path / f"{tumor}.{caller}.json")])
        paths.append(str(out))
    return paths


def test_union_groups_sites_with_provenance(synth, tmp_path):
    df = snv_union.build_union(_extracted(synth, tmp_path), "T1", gu.contig_order(fai=synth["fasta"] + ".fai"))
    keys = [tuple(x) for x in df[["chromosome", "position", "ref", "alt"]].itertuples(index=False)]
    assert keys == sorted(synth["snv_union"], key=lambda k: (k[0], k[1], k[3]))
    s1 = df.iloc[0]
    assert s1["callers"] == "caveman;mutect" and s1["n_callers"] == 2
    assert s1["source_samples"] == "T1" and s1["n_source_samples"] == 1
    assert s1["sources"] == "T1_vs_N1:caveman;T1_vs_N1:mutect"
    assert bool(s1["called_in_target"]) is True
    assert df[df.chromosome == "2"]["called_in_target"].tolist() == [False, False]


def test_union_sorts_by_fai_order_not_lexically(tmp_path):
    t = tmp_path / "a.tsv"
    t.write_text("chrom\tpos\tref\talt\ttumor_id\tnormal_id\tcaller\n"
                 "10\t5\tA\tG\tT\tN\tmutect\n2\t9\tA\tG\tT\tN\tmutect\n2\t3\tC\tT\tT\tN\tmutect\n")
    df = snv_union.build_union([str(t)], "X", {"2": 0, "10": 1})
    assert df[["chromosome", "position"]].values.tolist() == [["2", 3], ["2", 9], ["10", 5]]


def test_write_outputs_chunks_manifest_and_indexed_vcf(synth, tmp_path):
    df = snv_union.build_union(_extracted(synth, tmp_path), "T1", gu.contig_order(fai=synth["fasta"] + ".fai"))
    out = tmp_path / "out"
    out.mkdir()
    snv_union.write_outputs(df, "T1", str(out), chunk_size=3, contigs_order=gu.contig_order(
        fai=synth["fasta"] + ".fai"))
    assert gu.read_manifest(str(out / "chunks.manifest")) == {"chunk_00000": 3, "chunk_00001": 1}
    chunk = pd.read_csv(out / "chunks" / "chunk_00001.tsv", sep="\t", dtype={"chromosome": str})
    assert list(chunk.columns) == ["chromosome", "position", "ref", "alt"]
    union_tsv = pd.read_csv(out / "T1_union_snvs.tsv.gz", sep="\t", dtype={"chromosome": str})
    assert len(union_tsv) == 4
    vcf = pysam.VariantFile(str(out / "T1_union_snvs.vcf.gz"))
    assert os.path.exists(str(out / "T1_union_snvs.vcf.gz.tbi"))
    recs = list(vcf.fetch("1"))
    assert len(recs) == 2
    assert recs[0].info["NCALLERS"] == 2 and recs[0].info["INTARGET"] is True
    assert list(recs[0].info["CALLERS"]) == ["caveman", "mutect"]


def test_empty_union_still_writes_one_header_only_chunk(tmp_path):
    t = tmp_path / "a.tsv"
    t.write_text("chrom\tpos\tref\talt\ttumor_id\tnormal_id\tcaller\n")
    df = snv_union.build_union([str(t)], "X", {"1": 0})
    assert len(df) == 0
    snv_union.write_outputs(df, "X", str(tmp_path), chunk_size=10, contigs_order={"1": 0})
    assert gu.read_manifest(str(tmp_path / "chunks.manifest")) == {"chunk_00000": 0}
    assert (tmp_path / "chunks" / "chunk_00000.tsv").read_text().strip() == "chromosome\tposition\tref\talt"
    with gzip.open(tmp_path / "X_union_snvs.vcf.gz", "rt") as fh:
        assert any(l.startswith("#CHROM") for l in fh)
