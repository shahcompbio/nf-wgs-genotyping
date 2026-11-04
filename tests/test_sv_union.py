import gzip
import os
import random

import pandas as pd
import pytest

import geno_utils as gu
import sv_normalize
import sv_union

ORDER = {"1": 0, "2": 1, "3": 2}
PRIORITY = ["gridss", "brass", "svaba"]


def _row(c1, p1, s1, c2, p2, s2, caller="brass", tumor="T1", rid="r", imprecise=False, ci=0):
    return {"chrom1": c1, "pos1": p1, "strand1": s1, "chrom2": c2, "pos2": p2, "strand2": s2,
            "ci1_lo": -ci, "ci1_hi": ci, "ci2_lo": -ci, "ci2_hi": ci, "imprecise": imprecise,
            "sv_class": sv_normalize.sv_class(c1, s1, c2, s2), "caller": caller, "tumor_id": tumor,
            "normal_id": "N1", "record_id": rid, "filter": "PASS", "qual": ".", "ins_len": 0}


def _union(rows, slop=100, sample="T1"):
    return sv_union.cluster_union(rows, sample, ORDER, slop, PRIORITY)


def test_calls_within_slop_on_same_strands_merge_with_provenance():
    rows = [_row("1", 1000, "+", "1", 5000, "-", "brass", "T1", "b1"),
            _row("1", 1015, "+", "1", 5015, "-", "svaba", "T1", "s1", imprecise=True, ci=150),
            _row("1", 1030, "+", "1", 4990, "-", "brass", "T2", "b9")]
    union, members = _union(rows)
    assert len(union) == 1
    u = union.iloc[0]
    assert u["sv_id"] == "SV000001"
    assert (u["position_1"], u["position_2"], u["rep_caller"]) == (1000, 5000, "brass")
    assert u["callers"] == "brass;svaba" and u["n_callers"] == 2
    assert u["source_samples"] == "T1;T2" and bool(u["called_in_target"]) is True
    assert (u["spread_1"], u["spread_2"]) == (30, 25)
    assert u["n_members"] == 3
    assert set(members["sv_id"]) == {"SV000001"} and len(members) == 3


def test_slop_boundary_is_inclusive():
    a = _row("1", 1000, "+", "1", 5000, "-", rid="a")
    assert len(_union([a, _row("1", 1100, "+", "1", 5000, "-", rid="b")], slop=100)[0]) == 1
    assert len(_union([a, _row("1", 1101, "+", "1", 5000, "-", rid="b")], slop=100)[0]) == 2
    assert len(_union([a, _row("1", 1000, "+", "1", 5101, "-", rid="b")], slop=100)[0]) == 2


def test_opposite_strands_never_merge():
    rows = [_row("1", 1000, "+", "1", 5000, "-"), _row("1", 1000, "-", "1", 5000, "+", rid="x"),
            _row("1", 1000, "+", "1", 5000, "+", rid="y")]
    assert len(_union(rows)[0]) == 3


def test_representative_prefers_precise_then_priority():
    rows = [_row("1", 1000, "+", "1", 5000, "-", "gridss", rid="g", imprecise=True, ci=50),
            _row("1", 1010, "+", "1", 5010, "-", "svaba", rid="s"),
            _row("1", 1020, "+", "1", 5020, "-", "brass", rid="b")]
    u = _union(rows)[0].iloc[0]
    assert u["rep_caller"] == "brass" and u["position_1"] == 1020


def test_clustering_is_independent_of_input_order():
    rng = random.Random(1)
    rows = [_row("1", 1000 + 60 * i, "+", "2", 7000 + 60 * i, "-", rid=f"r{i}") for i in range(8)]
    rows += [_row("3", 50 + 500 * i, "-", "3", 9000, "+", rid=f"q{i}") for i in range(4)]
    base = _union(rows)[0]
    for _ in range(5):
        rng.shuffle(rows)
        pd.testing.assert_frame_equal(_union(rows)[0], base)


def test_ids_follow_genomic_order():
    rows = [_row("2", 10, "+", "2", 900, "-", rid="late"), _row("1", 500, "+", "3", 10, "-", rid="mid"),
            _row("1", 100, "+", "1", 900, "-", rid="early")]
    u = _union(rows)[0]
    assert u[["sv_id", "chromosome_1", "position_1"]].values.tolist() == [
        ["SV000001", "1", 100], ["SV000002", "1", 500], ["SV000003", "2", 10]]
    assert u["svtyper_type"].tolist() == ["DEL", "BND", "DEL"]
    assert u["sv_length"].tolist()[0] == 800 and pd.isna(u["sv_length"].tolist()[1])


def _records(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = fh.read().splitlines()
    head = [l for l in lines if l.startswith("##")]
    recs = [l.split("\t") for l in lines if not l.startswith("#")]
    return head, recs


def test_svtyper_vcf_conventions(tmp_path):
    rows = [_row("1", 100, "+", "1", 900, "-", rid="del", ci=5),
            _row("1", 2000, "-", "1", 2600, "+", rid="dup"),
            _row("2", 100, "+", "2", 700, "+", rid="inv"),
            _row("2", 1000, "-", "3", 400, "-", rid="tra")]
    union, members = _union(rows)
    path = str(tmp_path / "u.vcf")
    sv_union.write_svtyper_vcf(union, path, ORDER)
    head, recs = _records(path)
    by_id = {r[2]: r for r in recs}
    d = by_id["SV000001"]
    assert (d[1], d[4]) == ("100", "<DEL>") and "END=899" in d[7] and "CIPOS=-5,5" in d[7] and "CIEND=-5,5" in d[7]
    u = by_id["SV000002"]
    assert (u[1], u[4]) == ("1999", "<DUP>") and "END=2600" in u[7]
    i1, i2 = by_id["SV000003_1"], by_id["SV000003_2"]
    assert (i1[0], i1[1], i1[4]) == ("2", "100", "N]2:700]")
    assert (i2[0], i2[1], i2[4]) == ("2", "700", "N]2:100]")
    t1, t2 = by_id["SV000004_1"], by_id["SV000004_2"]
    # '-' breakends are written one base left of the junction (svtyper adds 1 back)
    assert (t1[0], t1[1], t1[4]) == ("2", "999", "[3:399[N")
    assert (t2[0], t2[1], t2[4]) == ("3", "399", "[2:999[N")
    for r in recs:
        assert "CIPOS=" in r[7]
        if "SVTYPE=BND" in r[7]:
            assert "MATEID=" in r[7]
    assert any(h.startswith("##INFO=<ID=CIPOS") for h in head)


def test_chunks_keep_bnd_mates_together(tmp_path):
    rows = [_row("1", 100 + 1000 * i, "-", "2", 400, "-", rid=f"t{i}") for i in range(3)]
    union, members = _union(rows)
    sv_union.write_outputs(union, members, "T1", str(tmp_path), chunk_size=2, order=ORDER)
    assert gu.read_manifest(str(tmp_path / "chunks.manifest")) == {"chunk_00000": 2, "chunk_00001": 1}
    for chunk in ("chunk_00000", "chunk_00001"):
        _, recs = _records(str(tmp_path / "chunks" / f"{chunk}.vcf"))
        ids = {r[2] for r in recs}
        assert all((i[:-2] + "_1") in ids and (i[:-2] + "_2") in ids for i in ids)
    assert os.path.exists(tmp_path / "T1_union_svs.vcf.gz.tbi")
    assert len(pd.read_csv(tmp_path / "T1_union_svs.bedpe.gz", sep="\t")) == 3
    assert len(pd.read_csv(tmp_path / "T1_union_svs.members.tsv.gz", sep="\t")) == 3


def test_empty_union_writes_one_empty_chunk(tmp_path):
    union, members = _union([])
    assert len(union) == 0
    sv_union.write_outputs(union, members, "T1", str(tmp_path), chunk_size=10, order=ORDER)
    assert gu.read_manifest(str(tmp_path / "chunks.manifest")) == {"chunk_00000": 0}
    head, recs = _records(str(tmp_path / "chunks" / "chunk_00000.vcf"))
    assert recs == [] and head


def test_synthetic_fixture_union(synth, tmp_path):
    tsvs = []
    for (who, caller), vcf in sorted(synth["sv_vcfs"].items()):
        tumor = synth["sv_target"] if who == "T" else synth["other_tumor"]
        out = tmp_path / f"{who}.{caller}.tsv"
        sv_normalize.main(["--vcf", vcf, "--caller", caller, "--tumor-id", tumor, "--normal-id", "N1",
                           "--contigs", ",".join(str(i) for i in range(1, 23)), "--out-tsv", str(out),
                           "--out-stats", str(tmp_path / f"{who}.{caller}.json")])
        tsvs.append(str(out))
    od = tmp_path / "u"
    od.mkdir()
    sv_union.main(["--inputs", *tsvs, "--sample-id", synth["sv_target"],
                   "--contigs", ",".join(str(i) for i in range(1, 23)), "--outdir", str(od)])
    union = pd.read_csv(od / f"{synth['sv_target']}_union_svs.bedpe.gz", sep="\t", dtype={"chromosome_1": str})
    exp = synth["sv_expected"]
    assert len(union) == exp["n_union"]
    assert int(union["called_in_target"].sum()) == len(exp["in_target"])
    shared = union[union["n_source_samples"] == 2]
    assert len(shared) == 1 and shared.iloc[0]["rep_caller"] == "gridss"
