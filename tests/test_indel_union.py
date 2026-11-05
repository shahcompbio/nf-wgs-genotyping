import os

import pandas as pd
import pysam

import geno_utils as gu
import indel_extract
import indel_union

KEYS = ["chromosome", "position", "ref", "alt"]


def _extracted(synth, tmp_path):
    paths = []
    for caller, tumor in (("pindel", "T1"), ("strelka", "T1"), ("mutect", "T2")):
        out = tmp_path / f"{tumor}.{caller}.indels.tsv"
        indel_extract.main(["--vcf", synth["indel_vcfs"][caller], "--caller", caller, "--tumor-id", tumor,
                            "--normal-id", "N1", "--contigs", "1,2,3", "--reference", synth["fasta"],
                            "--out-tsv", str(out), "--out-stats", str(tmp_path / f"{tumor}.{caller}.json")])
        paths.append(str(out))
    return paths


def _union(synth, tmp_path):
    order = gu.contig_order(fai=synth["fasta"] + ".fai")
    df = indel_union.build_union(_extracted(synth, tmp_path), "T1", order)
    return indel_union.annotate_loci(df, synth["fasta"]), order


def test_union_merges_callers_and_normalized_representations(synth, tmp_path):
    df, _ = _union(synth, tmp_path)
    s = synth["indel_sites"]
    names = ["D1", "HP1", "HP2", "I3", "C4"]
    assert [tuple(x) for x in df[KEYS].itertuples(index=False)] == [tuple(s[n]) for n in names]
    for name, row in zip(names, df.itertuples(index=False)):
        exp = synth["indel_union"][name]
        assert (row.callers, row.source_samples) == (exp["callers"], exp["source_samples"]), name
    assert df["called_in_target"].tolist() == [True, True, False, True, True]
    assert df.iloc[3]["sources"] == "T1_vs_N1:pindel;T1_vs_N1:strelka"


def test_loci_footprints_cover_repeats_and_group_overlapping_alleles(synth, tmp_path):
    df, _ = _union(synth, tmp_path)
    assert df["indel_length"].tolist() == [-3, -1, -2, 2, -3]
    assert list(zip(df["var_start"], df["var_end"])) == [(500, 504), (1000, 1009), (1000, 1009),
                                                       (1500, 1513), (2000, 2005)]
    assert df["locus_id"].tolist() == ["3:500-504", "3:1000-1009", "3:1000-1009", "3:1500-1513",
                                       "3:2000-2005"]


def test_alleles_sharing_one_flanking_base_are_separate_loci(synth, tmp_path):
    seq = synth["seqs"]["3"]
    assert seq[503:505] == "GG"                     # 504 G, 505 G: +T or +A after 504 is unique
    t = tmp_path / "x.tsv"
    t.write_text("chrom\tpos\tref\talt\ttumor_id\tnormal_id\tcaller\n"
                 "3\t500\tGTCA\tG\tT1\tN1\tpindel\n"
                 "3\t504\tG\tGT\tT1\tN1\tpindel\n"
                 "3\t504\tG\tGA\tT1\tN1\tpindel\n")
    df = indel_union.annotate_loci(indel_union.build_union([str(t)], "T1", {"3": 0}), synth["fasta"])
    assert df["locus_id"].tolist() == ["3:500-504", "3:504-505", "3:504-505"]


def test_write_outputs_keeps_loci_whole_in_chunks(synth, tmp_path):
    df, order = _union(synth, tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    indel_union.write_outputs(df, "T1", str(out), chunk_size=2, contigs_order=order,
                              fai=synth["fasta"] + ".fai")
    assert gu.read_manifest(str(out / "chunks.manifest")) == {"chunk_00000": 1, "chunk_00001": 2,
                                                             "chunk_00002": 2}
    chunk = pd.read_csv(out / "chunks" / "chunk_00001.tsv", sep="\t", dtype={"chromosome": str})
    assert list(chunk.columns) == KEYS + ["locus_id"]
    assert chunk["locus_id"].tolist() == ["3:1000-1009", "3:1000-1009"]
    union = pd.read_csv(out / "T1_union_indels.tsv.gz", sep="\t", dtype={"chromosome": str})
    assert len(union) == 5 and "locus_id" in union.columns
    assert os.path.exists(out / "T1_union_indels.vcf.gz.tbi")
    with pysam.VariantFile(str(out / "T1_union_indels.vcf.gz")) as vf:
        recs = list(vf.fetch("3"))
        assert "nf-wgs-genotyping indel_union" in str(vf.header)
    assert [(r.pos, r.ref, r.alts[0]) for r in recs] == [tuple(x[1:]) for x in df[KEYS].itertuples(index=False)]


def test_empty_union_writes_one_header_only_chunk(synth, tmp_path):
    t = tmp_path / "a.tsv"
    t.write_text("chrom\tpos\tref\talt\ttumor_id\tnormal_id\tcaller\n")
    out = tmp_path / "out"
    indel_union.main(["--inputs", str(t), "--sample-id", "X", "--reference", synth["fasta"],
                      "--outdir", str(out)])
    assert gu.read_manifest(str(out / "chunks.manifest")) == {"chunk_00000": 0}
    assert (out / "chunks" / "chunk_00000.tsv").read_text().strip() == "\t".join(KEYS + ["locus_id"])
    assert len(pd.read_csv(out / "X_union_indels.tsv.gz", sep="\t")) == 0
