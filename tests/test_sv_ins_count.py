"""SV insertion genotyping by realignment (bin/sv_ins_count.py) on the contig-4 fixtures."""
import pandas as pd
import pytest

import sv_ins_count

SITE_COLUMNS = ["adj_id", "chromosome", "position", "del_len", "ins_len", "ins_seq", "ins_seq_available",
                "imprecise"]


def _site(adj_id, chrom, pos, ins_len, ins_seq="", imprecise=False, del_len=0):
    return {"adj_id": adj_id, "chromosome": chrom, "position": pos, "del_len": del_len, "ins_len": ins_len,
            "ins_seq": ins_seq, "ins_seq_available": bool(ins_seq), "imprecise": imprecise}


def _fixture_sites(synth):
    s = synth["ins_sites"]
    c1, a1, _, seq1 = s["INS1"]
    c2, a2, _, seq2 = s["INS2"]
    return [_site("SV1.1", c1, a1, 60, seq1),                       # sequence-resolved
            _site("SV1.2", c1, a1, 60),                             # length only: N block
            _site("SV2.1", c2, a2, 50, seq2),                       # in a repeat
            _site("SV3.1", "4", s["IMPRECISE"][1], 300, imprecise=True),
            _site("SV4.1", "4", s["NO_EVIDENCE"][1], 100)]


def _by_id(rows):
    return {r["adj_id"]: r for r in rows}


COUNTS = ("ref_fragments", "alt_fragments", "uninformative_fragments", "conflict_fragments")


def test_fixture_truth_for_sequence_and_length_only_haplotypes(synth):
    out = _by_id(sv_ins_count.count_insertions(synth["bam"], synth["fasta"], _fixture_sites(synth)))
    truth = synth["truth_ins"]
    for adj in ("SV1.1", "SV1.2"):
        assert {k: out[adj][k] for k in COUNTS} == truth["INS1"], adj
        assert out[adj]["ins_status"] == "genotyped"
    assert {k: out["SV2.1"][k] for k in COUNTS} == truth["INS2"]
    assert {k: out["SV4.1"][k] for k in COUNTS} == truth["NO_EVIDENCE"]
    assert out["SV4.1"]["ins_status"] == "no_evidence"
    assert out["SV3.1"]["ins_status"] == "skipped_imprecise"
    assert all(out["SV3.1"][k] is None for k in COUNTS)


def test_without_reference_nothing_is_genotyped(synth):
    out = _by_id(sv_ins_count.count_insertions(synth["bam"], None, _fixture_sites(synth)))
    assert out["SV3.1"]["ins_status"] == "skipped_imprecise"
    assert {out[a]["ins_status"] for a in ("SV1.1", "SV1.2", "SV2.1", "SV4.1")} == {"skipped_no_reference"}
    assert all(out["SV1.1"][k] is None for k in COUNTS)


def test_long_insertions_use_a_bounded_haplotype(synth):
    c1, a1, _, _ = synth["ins_sites"]["INS1"]
    # a 5-kb length-only insertion at INS1: the reads cannot tell it from the 60-bp one
    out = sv_ins_count.count_insertions(synth["bam"], synth["fasta"], [_site("L", c1, a1, 5000)])[0]
    assert {k: out[k] for k in COUNTS} == synth["truth_ins"]["INS1"]
    block = sv_ins_count.alt_block("A" * 1000 + "C" * 3000 + "G" * 1000, max_read_len=50)
    assert block == "A" * 60 + "G" * 60
    assert sv_ins_count.alt_block("ACGT" * 10, max_read_len=50) == "ACGT" * 10


def test_site_without_reads_gives_no_evidence(synth):
    out = sv_ins_count.count_insertions(synth["bam"], synth["fasta"], [_site("X", "4", 3800, 80)])[0]
    assert out["ins_status"] == "no_evidence" and out["alt_fragments"] == 0


def test_cli_writes_counts_in_site_order(synth, tmp_path):
    sites = tmp_path / "ins_chunk_00000.tsv"
    pd.DataFrame(_fixture_sites(synth), columns=SITE_COLUMNS).to_csv(sites, sep="\t", index=False)
    out = tmp_path / "ins_counts_ins_chunk_00000.tsv"
    sv_ins_count.main(["--bam", synth["bam"], "--reference", synth["fasta"], "--sites", str(sites),
                       "--out", str(out)])
    df = pd.read_csv(out, sep="\t", keep_default_na=False)
    assert list(df.columns) == sv_ins_count.COLUMNS
    assert df["adj_id"].tolist() == ["SV1.1", "SV1.2", "SV2.1", "SV3.1", "SV4.1"]
    assert df.loc[3, "ins_status"] == "skipped_imprecise" and df.loc[3, "alt_fragments"] == ""
    sv_ins_count.main(["--bam", synth["bam"], "--sites", str(sites), "--out", str(out)])   # no FASTA
    assert set(pd.read_csv(out, sep="\t")["ins_status"]) == {"skipped_imprecise", "skipped_no_reference"}


def test_empty_chunk(synth, tmp_path):
    sites = tmp_path / "ins_chunk_00000.tsv"
    pd.DataFrame(columns=SITE_COLUMNS).to_csv(sites, sep="\t", index=False)
    out = tmp_path / "o.tsv"
    sv_ins_count.main(["--bam", synth["bam"], "--reference", synth["fasta"], "--sites", str(sites),
                       "--out", str(out)])
    assert out.read_text().splitlines() == ["\t".join(sv_ins_count.COLUMNS)]


@pytest.mark.parametrize("bad", ["0", "-1"])
def test_min_score_diff_must_be_positive(synth, tmp_path, bad):
    with pytest.raises(SystemExit):
        sv_ins_count.main(["--bam", synth["bam"], "--sites", "x", "--out", "y", "--min-score-diff", bad])
