import pandas as pd

import snv_count


def _sites(synth):
    s3_ref = synth["seqs"]["2"][synth["S3"][1] - 1]
    s3_alt = [b for b in "ACGT" if b != s3_ref][0]
    s2_ref = synth["seqs"]["1"][synth["S2"][1] - 1]
    return [
        ("1", synth["S1"][1], synth["S1_ref"], synth["S1_alt"]),
        ("1", synth["S2"][1], s2_ref, [b for b in "ACGT" if b != s2_ref][0]),
        ("2", synth["S3"][1], s3_ref, s3_alt),
    ]


def _counts(synth, **kw):
    return {(r["chromosome"], r["position"]): r for r in snv_count.count_sites(synth["bam"], _sites(synth), **kw)}


def test_s1_counts_fragments_once_with_explicit_filters(synth):
    c = _counts(synth)[("1", synth["S1"][1])]
    assert {k: c[k] for k in synth["truth_S1"]} == synth["truth_S1"]


def test_reference_only_and_uncovered_sites(synth):
    c = _counts(synth)
    assert {k: c[("1", synth["S2"][1])][k] for k in synth["truth_S2"]} == synth["truth_S2"]
    s3 = c[("2", synth["S3"][1])]
    assert (s3["ref_count"], s3["alt_count"], s3["other_count"], s3["conflict_count"]) == (0, 0, 0, 0)


def test_count_duplicates_and_baseq_thresholds_are_parameters(synth):
    s1 = ("1", synth["S1"][1])
    assert _counts(synth, count_duplicates=True)[s1]["alt_count"] == synth["truth_S1"]["alt_count"] + 3
    assert _counts(synth, min_baseq=10)[s1]["alt_count"] == synth["truth_S1"]["alt_count"] + 2
    assert _counts(synth, min_mapq=10)[s1]["alt_count"] == synth["truth_S1"]["alt_count"] + 2


def test_cli_reads_chunk_and_writes_counts(synth, tmp_path):
    chunk = tmp_path / "chunk_00007.tsv"
    pd.DataFrame(_sites(synth), columns=["chromosome", "position", "ref", "alt"]).to_csv(chunk, sep="\t", index=False)
    out = tmp_path / "counts_chunk_00007.tsv"
    snv_count.main(["--bam", synth["bam"], "--sites", str(chunk), "--out", str(out)])
    df = pd.read_csv(out, sep="\t", dtype={"chromosome": str})
    assert list(df.columns) == ["chromosome", "position", "ref", "alt", "ref_count", "alt_count",
                                "other_count", "conflict_count"]
    assert df["alt_count"].tolist()[0] == synth["truth_S1"]["alt_count"]


def test_cli_header_only_chunk_gives_header_only_counts(synth, tmp_path):
    chunk = tmp_path / "chunk_00000.tsv"
    chunk.write_text("chromosome\tposition\tref\talt\n")
    out = tmp_path / "counts_chunk_00000.tsv"
    snv_count.main(["--bam", synth["bam"], "--sites", str(chunk), "--out", str(out)])
    assert out.read_text().strip().split("\t")[0] == "chromosome"
    assert len(out.read_text().strip().splitlines()) == 1
