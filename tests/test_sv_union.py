import gzip
import os
import random

import pandas as pd
import pytest

import geno_utils as gu
import sv_normalize
import sv_union

ORDER = {"1": 0, "2": 1, "3": 2, "4": 3}
PRIORITY = ["gridss", "brass", "svaba"]


def _row(c1, p1, s1, c2, p2, s2, caller="brass", tumor="T1", rid="r", imprecise=False, ci=0):
    return {"chrom1": c1, "pos1": p1, "strand1": s1, "chrom2": c2, "pos2": p2, "strand2": s2,
            "ci1_lo": -ci, "ci1_hi": ci, "ci2_lo": -ci, "ci2_hi": ci, "imprecise": imprecise,
            "sv_class": sv_normalize.sv_class(c1, s1, c2, s2), "caller": caller, "tumor_id": tumor,
            "normal_id": "N1", "record_id": rid, "filter": "PASS", "qual": ".", "ins_len": 0,
            "ins_seq": "", "ins_seq_available": True}


def _ins(chrom, pos, ins_len, seq="", caller="generic", tumor="T1", rid="i", imprecise=False, del_len=0):
    r = _row(chrom, pos, "+", chrom, pos + del_len + 1, "-", caller, tumor, rid, imprecise)
    r.update(sv_class="INS", ins_len=ins_len, ins_seq=seq, ins_seq_available=bool(seq))
    return r


def _union(rows, slop=100, sample="T1"):
    return sv_union.cluster_union(rows, sample, ORDER, slop, PRIORITY)


def test_calls_within_slop_on_same_strands_merge_with_provenance():
    rows = [_row("1", 1000, "+", "1", 5000, "-", "brass", "T1", "b1"),
            _row("1", 1015, "+", "1", 5015, "-", "svaba", "T1", "s1", imprecise=True, ci=150),
            _row("1", 1030, "+", "1", 4990, "-", "brass", "T2", "b9")]
    union, adjs, members = _union(rows)
    assert len(union) == 1
    u = union.iloc[0]
    assert u["sv_id"] == "SV000001"
    assert (u["position_1"], u["position_2"], u["rep_caller"]) == (1000, 5000, "brass")
    assert u["callers"] == "brass;svaba" and u["n_callers"] == 2
    assert u["source_samples"] == "T1;T2" and bool(u["called_in_target"]) is True
    assert (u["spread_1"], u["spread_2"]) == (30, 25)
    assert u["n_members"] == 3
    assert set(members["sv_id"]) == {"SV000001"} and len(members) == 3
    # every distinct breakpoint is its own adjacency (genotyped separately), mapped to its source
    assert (u["n_adjacencies"], u["rep_adj_id"], u["n_linked_clusters"]) == (3, "SV000001.1", 0)
    assert adjs[["adj_id", "position_1", "position_2"]].values.tolist() == [
        ["SV000001.1", 1000, 5000], ["SV000001.2", 1015, 5015], ["SV000001.3", 1030, 4990]]
    assert adjs["is_cluster_rep"].tolist() == [True, False, False]
    assert adjs["called_in_target"].tolist() == [True, True, False]
    assert adjs["sources"].tolist() == ["T1_vs_N1:brass:b1", "T1_vs_N1:svaba:s1", "T2_vs_N1:brass:b9"]
    assert dict(zip(members["record_id"], members["adj_id"])) == {"b1": "SV000001.1", "s1": "SV000001.2",
                                                                  "b9": "SV000001.3"}


def test_identical_breakpoints_share_one_adjacency_with_the_best_members_ci():
    rows = [_row("1", 1000, "+", "1", 5000, "-", "svaba", "T1", "s1", imprecise=True, ci=150),
            _row("1", 1000, "+", "1", 5000, "-", "brass", "T2", "b1", ci=5),
            _row("1", 1040, "+", "1", 5000, "-", "gridss", "T1", "g1")]
    union, adjs, members = _union(rows)
    assert len(union) == 1 and len(adjs) == 2 and len(members) == 3
    a = adjs.set_index("position_1")
    assert a.loc[1000, "n_members"] == 2 and a.loc[1000, "callers"] == "brass;svaba"
    assert (a.loc[1000, "imprecise"], a.loc[1000, "ci1_lo"], a.loc[1000, "ci1_hi"]) == (False, -5, 5)
    assert union.iloc[0]["rep_adj_id"] == a.loc[1040, "adj_id"]         # gridss first by priority


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


def test_chains_do_not_merge_transitively():
    # 8 junctions 60 bp apart: single linkage would make one cluster spanning 420 bp
    rows = [_row("1", 1000 + 60 * i, "+", "2", 7000 + 60 * i, "-", rid=f"r{i}") for i in range(8)]
    union, adjs, _ = _union(rows)
    assert len(union) == 4 and len(adjs) == 8
    assert (union["spread_1"] <= 100).all() and (union["spread_2"] <= 100).all()
    assert union["n_members"].tolist() == [2, 2, 2, 2]
    # neighboring clusters within slop of each other are flagged for review
    assert union["n_linked_clusters"].tolist() == [1, 2, 2, 1]


def test_every_member_is_within_slop_of_every_other():
    rng = random.Random(7)
    rows = [_row("1", rng.randint(1000, 1600), "+", "1", rng.randint(9000, 9600), "-", rid=f"r{i}")
            for i in range(60)]
    union, adjs, members = _union(rows)
    for _, m in members.groupby("sv_id"):
        assert m["pos1"].max() - m["pos1"].min() <= 100 and m["pos2"].max() - m["pos2"].min() <= 100


def test_clustering_is_independent_of_input_order():
    rng = random.Random(1)
    rows = [_row("1", 1000 + 60 * i, "+", "2", 7000 + 60 * i, "-", rid=f"r{i}") for i in range(8)]
    rows += [_row("3", 50 + 500 * i, "-", "3", 9000, "+", rid=f"q{i}") for i in range(4)]
    rows += [_row("2", 300 + rng.randint(0, 150), "+", "3", 800 + rng.randint(0, 150), "+", rid=f"z{i}")
             for i in range(20)]
    base = _union(rows)
    for _ in range(5):
        rng.shuffle(rows)
        for got, exp in zip(_union(rows), base):
            pd.testing.assert_frame_equal(got, exp)


def test_ids_follow_genomic_order():
    rows = [_row("2", 10, "+", "2", 900, "-", rid="late"), _row("1", 500, "+", "3", 10, "-", rid="mid"),
            _row("1", 100, "+", "1", 900, "-", rid="early")]
    u = _union(rows)[0]
    assert u[["sv_id", "chromosome_1", "position_1"]].values.tolist() == [
        ["SV000001", "1", 100], ["SV000002", "1", 500], ["SV000003", "2", 10]]
    assert u["svtyper_type"].tolist() == ["DEL", "BND", "DEL"]
    assert u["sv_length"].tolist()[0] == 800 and pd.isna(u["sv_length"].tolist()[1])


def test_insertions_cluster_apart_and_go_to_insertion_chunks(tmp_path):
    seq = "ACGT" * 15
    rows = [_ins("1", 1000, 60, seq, caller="gridss", rid="g"),
            _ins("1", 1000, 60, tumor="T2", rid="len"),          # length-only call of the same event
            _ins("1", 1030, 300, imprecise=True, rid="imp"),
            _row("1", 1000, "+", "1", 5000, "-", rid="del")]
    union, adjs, members = _union(rows)
    assert len(union) == 2
    ins = union[union["sv_class"] == "INS"].iloc[0]
    assert (ins["genotyper"], ins["n_adjacencies"], ins["sv_length"]) == ("realign", 3, 60)
    assert pd.isna(ins["svtyper_type"]) or ins["svtyper_type"] == ""
    a = adjs[adjs["sv_id"] == ins["sv_id"]].set_index("ins_seq_available")
    rep = adjs.set_index("adj_id").loc[ins["rep_adj_id"]]
    assert rep["ins_seq"] == seq                                # sequence-resolved preferred
    sv_union.write_outputs(union, adjs, members, "T1", str(tmp_path), chunk_size=10, order=ORDER,
                           ins_chunk_size=2)
    assert gu.read_manifest(str(tmp_path / "chunks.manifest")) == {"chunk_00000": 1}
    assert gu.read_manifest(str(tmp_path / "ins_chunks.manifest")) == {"ins_chunk_00000": 2, "ins_chunk_00001": 1}
    sites = pd.concat([pd.read_csv(tmp_path / "ins_chunks" / f, sep="\t", keep_default_na=False)
                       for f in ("ins_chunk_00000.tsv", "ins_chunk_00001.tsv")])
    assert list(sites.columns) == ["adj_id", "chromosome", "position", "del_len", "ins_len", "ins_seq",
                                   "ins_seq_available", "imprecise"]
    assert sorted(sites["ins_len"].tolist()) == [60, 60, 300] and set(sites["del_len"]) == {0}
    _, recs = _records(str(tmp_path / "T1_union_svs.vcf.gz"))
    assert [r[2] for r in recs] == [union[union["sv_class"] == "DEL"].iloc[0]["rep_adj_id"]]


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
    union, adjs, members = _union(rows)
    path = str(tmp_path / "u.vcf")
    sv_union.write_svtyper_vcf(adjs, path, ORDER)
    head, recs = _records(path)
    by_id = {r[2]: r for r in recs}
    d = by_id["SV000001.1"]
    assert (d[1], d[4]) == ("100", "<DEL>") and "END=899" in d[7] and "CIPOS=-5,5" in d[7] and "CIEND=-5,5" in d[7]
    assert "EVENT=SV000001.1" in d[7] and "CLUSTER=SV000001" in d[7]
    u = by_id["SV000002.1"]
    assert (u[1], u[4]) == ("1999", "<DUP>") and "END=2600" in u[7]
    i1, i2 = by_id["SV000003.1_1"], by_id["SV000003.1_2"]
    assert (i1[0], i1[1], i1[4]) == ("2", "100", "N]2:700]")
    assert (i2[0], i2[1], i2[4]) == ("2", "700", "N]2:100]")
    t1, t2 = by_id["SV000004.1_1"], by_id["SV000004.1_2"]
    # '-' breakends are written one base left of the junction (svtyper adds 1 back)
    assert (t1[0], t1[1], t1[4]) == ("2", "999", "[3:399[N")
    assert (t2[0], t2[1], t2[4]) == ("3", "399", "[2:999[N")
    for r in recs:
        assert "CIPOS=" in r[7]
        if "SVTYPE=BND" in r[7]:
            assert "MATEID=" in r[7]
    assert any(h.startswith("##INFO=<ID=CIPOS") for h in head)
    assert any(h.startswith("##INFO=<ID=CLUSTER") for h in head)


def test_chunks_keep_bnd_mates_together(tmp_path):
    rows = [_row("1", 100 + 1000 * i, "-", "2", 400, "-", rid=f"t{i}") for i in range(3)]
    # one cluster with five distinct translocation breakpoints: five BND pairs
    rows += [_row("1", 5000 + 10 * i, "+", "3", 800 + 5 * i, "-", rid=f"c{i}") for i in range(5)]
    union, adjs, members = _union(rows)
    assert len(union) == 4 and len(adjs) == 8
    sv_union.write_outputs(union, adjs, members, "T1", str(tmp_path), chunk_size=3, order=ORDER)
    assert gu.read_manifest(str(tmp_path / "chunks.manifest")) == {"chunk_00000": 3, "chunk_00001": 3,
                                                                   "chunk_00002": 2}
    seen = []
    for chunk in ("chunk_00000", "chunk_00001", "chunk_00002"):
        _, recs = _records(str(tmp_path / "chunks" / f"{chunk}.vcf"))
        ids = {r[2] for r in recs}
        assert all((i[:-2] + "_1") in ids and (i[:-2] + "_2") in ids for i in ids)
        seen += sorted({i[:-2] for i in ids})
    assert sorted(seen) == sorted(adjs["adj_id"])
    assert os.path.exists(tmp_path / "T1_union_svs.vcf.gz.tbi")
    assert len(pd.read_csv(tmp_path / "T1_union_svs.bedpe.gz", sep="\t")) == 4
    assert len(pd.read_csv(tmp_path / "T1_union_sv_adjacencies.tsv.gz", sep="\t")) == 8
    m = pd.read_csv(tmp_path / "T1_union_svs.members.tsv.gz", sep="\t")
    assert len(m) == 8 and list(m.columns[:2]) == ["sv_id", "adj_id"]


def test_empty_union_writes_one_empty_chunk(tmp_path):
    union, adjs, members = _union([])
    assert len(union) == 0 and len(adjs) == 0
    sv_union.write_outputs(union, adjs, members, "T1", str(tmp_path), chunk_size=10, order=ORDER)
    assert gu.read_manifest(str(tmp_path / "chunks.manifest")) == {"chunk_00000": 0}
    assert gu.read_manifest(str(tmp_path / "ins_chunks.manifest")) == {"ins_chunk_00000": 0}
    head, recs = _records(str(tmp_path / "chunks" / "chunk_00000.vcf"))
    assert recs == [] and head
    assert len(pd.read_csv(tmp_path / "ins_chunks" / "ins_chunk_00000.tsv", sep="\t")) == 0


def _normalize_all(vcfs, tmp_path, contigs, fasta=None):
    tsvs = []
    for (tumor, caller, name), vcf in vcfs:
        out = tmp_path / f"{name}.tsv"
        sv_normalize.main(["--vcf", vcf, "--caller", caller, "--tumor-id", tumor, "--normal-id", "N1",
                           "--contigs", contigs, "--out-tsv", str(out), "--out-stats", str(tmp_path / f"{name}.json"),
                           *(["--reference", fasta] if fasta else [])])
        tsvs.append(str(out))
    return tsvs


def test_synthetic_fixture_union(synth, tmp_path):
    vcfs = [((synth["sv_target"] if who == "T" else synth["other_tumor"], caller, f"{who}.{caller}"), vcf)
            for (who, caller), vcf in sorted(synth["sv_vcfs"].items())]
    tsvs = _normalize_all(vcfs, tmp_path, ",".join(str(i) for i in range(1, 23)))
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
    # two deletions were called at shifted breakpoints by a second caller: 2 adjacencies each
    adjs = pd.read_csv(od / f"{synth['sv_target']}_union_sv_adjacencies.tsv.gz", sep="\t")
    assert len(adjs) == exp["n_union"] + 2 and sorted(union["n_adjacencies"].tolist()) == [1, 1, 1, 1, 2, 2]


def test_synthetic_insertion_union(synth, tmp_path):
    v = synth["sv_ins_vcfs"]
    vcfs = [(("T1", "gridss", "g"), v["gridss"]), (("T1", "generic", "t1"), v["T1_generic"]),
            (("T2", "generic", "t2"), v["T2_generic"])]
    tsvs = _normalize_all(vcfs, tmp_path, "1,2,3,4", fasta=synth["fasta"])
    od = tmp_path / "u"
    sv_union.main(["--inputs", *tsvs, "--sample-id", "T1", "--reference", synth["fasta"], "--outdir", str(od)])
    union = pd.read_csv(od / "T1_union_svs.bedpe.gz", sep="\t")
    adjs = pd.read_csv(od / "T1_union_sv_adjacencies.tsv.gz", sep="\t", keep_default_na=False)
    exp = synth["sv_ins_expected"]
    assert (len(union), len(adjs)) == (exp["n_clusters"], exp["n_adjacencies"])
    assert set(adjs["sv_class"]) == {"INS"} and set(adjs["genotyper"]) == {"realign"}
    ins1 = adjs[adjs["position_1"] == synth["ins_sites"]["INS1"][1]]
    assert sorted(ins1["ins_seq_available"].tolist()) == [False, True]
    assert sorted(ins1["called_in_target"].tolist()) == [False, True]
    ins2 = adjs[adjs["position_1"] == synth["ins_sites"]["INS2"][1]].iloc[0]
    assert ins2["callers"] == "generic" and ins2["source_samples"] == "T1;T2" and ins2["n_members"] == 2
    assert gu.read_manifest(str(od / "ins_chunks.manifest")) == {"ins_chunk_00000": 5}
    assert gu.read_manifest(str(od / "chunks.manifest")) == {"chunk_00000": 0}
