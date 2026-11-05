"""Fragment counting at indel loci by realignment to ref/alt haplotypes (truth: make_fixtures.py)."""
import pandas as pd
import pytest

import indel_count

NAMES = ["D1", "HP1", "HP2", "I3", "C4"]
LOCI = {"D1": "3:500-504", "HP1": "3:1000-1009", "HP2": "3:1000-1009", "I3": "3:1500-1513",
        "C4": "3:2000-2005"}
COUNTS = ["ref_count", "alt_count", "other_count", "uninformative_count", "conflict_count"]


def _sites(synth, names=NAMES):
    return [(*synth["indel_sites"][n], LOCI[n]) for n in names]


def _counts(synth, names=NAMES, **kw):
    rows = indel_count.count_sites(synth["bam"], synth["fasta"], _sites(synth, names), **kw)
    return dict(zip(names, rows))


def test_counts_reproduce_designed_truth(synth):
    got = _counts(synth)
    for name in NAMES:
        assert {k: got[name][k] for k in COUNTS} == synth["truth_indels"][name], name
        chrom, pos, ref, alt = synth["indel_sites"][name]
        assert (got[name]["chromosome"], got[name]["position"], got[name]["ref"], got[name]["alt"]) == \
            (chrom, pos, ref, alt)


def test_read_filters_and_thresholds_are_parameters(synth):
    t = synth["truth_indels"]["D1"]
    assert _counts(synth, ["D1"], count_duplicates=True)["D1"]["alt_count"] == t["alt_count"] + 2
    assert _counts(synth, ["D1"], min_mapq=10)["D1"]["alt_count"] == t["alt_count"] + 1
    lowbq = _counts(synth, ["D1"], min_baseq=0)["D1"]          # the BQ-2 read is no longer masked
    assert (lowbq["alt_count"], lowbq["uninformative_count"]) == (t["alt_count"] + 1, 0)
    strict = _counts(synth, ["D1"], min_score_diff=1000)["D1"]
    assert (strict["ref_count"], strict["alt_count"], strict["conflict_count"]) == (0, 0, 0)


def test_alleles_of_a_locus_are_genotyped_jointly(synth):
    # alone, HP1 cannot claim the -AA fragments (they are closer to -A than to the reference, but
    # by less than one discriminating base): they become uninformative, never HP1 alt
    alone = indel_count.count_sites(synth["bam"], synth["fasta"], [(*synth["indel_sites"]["HP1"], "solo")])[0]
    t = synth["truth_indels"]["HP1"]
    assert (alone["alt_count"], alone["other_count"]) == (t["alt_count"], 0)
    assert alone["uninformative_count"] == t["uninformative_count"] + t["other_count"]


def test_uncovered_locus_gives_zero_counts(synth):
    seq = synth["seqs"]["2"]
    site = ("2", 2000, seq[1999:2001], seq[1999], "2:2000-2002")
    row = indel_count.count_sites(synth["bam"], synth["fasta"], [site])[0]
    assert [row[k] for k in COUNTS] == [0, 0, 0, 0, 0]


def test_cli_reads_chunk_and_writes_counts(synth, tmp_path):
    chunk = tmp_path / "chunk_00003.tsv"
    pd.DataFrame(_sites(synth), columns=["chromosome", "position", "ref", "alt", "locus_id"]).to_csv(
        chunk, sep="\t", index=False)
    out = tmp_path / "counts_chunk_00003.tsv"
    indel_count.main(["--bam", synth["bam"], "--reference", synth["fasta"], "--sites", str(chunk),
                      "--out", str(out)])
    df = pd.read_csv(out, sep="\t", dtype={"chromosome": str})
    assert list(df.columns) == ["chromosome", "position", "ref", "alt"] + COUNTS
    assert df["alt_count"].tolist() == [synth["truth_indels"][n]["alt_count"] for n in NAMES]


def test_cli_header_only_chunk_gives_header_only_counts(synth, tmp_path):
    chunk = tmp_path / "chunk_00000.tsv"
    chunk.write_text("chromosome\tposition\tref\talt\tlocus_id\n")
    out = tmp_path / "counts_chunk_00000.tsv"
    indel_count.main(["--bam", synth["bam"], "--reference", synth["fasta"], "--sites", str(chunk),
                      "--out", str(out)])
    assert out.read_text().strip() == "\t".join(["chromosome", "position", "ref", "alt"] + COUNTS)


@pytest.mark.parametrize("bad", ["--min-score-diff=0", "--min-score-diff=-1"])
def test_cli_rejects_non_positive_score_threshold(synth, tmp_path, bad):
    with pytest.raises(SystemExit):
        indel_count.main(["--bam", synth["bam"], "--reference", synth["fasta"], "--sites", "x",
                          "--out", str(tmp_path / "o.tsv"), bad])
