import gzip
import json

import pandas as pd
import pytest
import yaml

import sv_finalize
import sv_ins_count
import sv_normalize
import sv_union

ORDER = {"1": 0, "2": 1}
FMT = "GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:FS:FP:FC:AOF:ROF:FX"
SEQ = "ACGT" * 15


def _row(c1, p1, s1, c2, p2, s2, rid):
    return {"chrom1": c1, "pos1": p1, "strand1": s1, "chrom2": c2, "pos2": p2, "strand2": s2,
            "ci1_lo": 0, "ci1_hi": 0, "ci2_lo": 0, "ci2_hi": 0, "imprecise": False,
            "sv_class": sv_normalize.sv_class(c1, s1, c2, s2), "caller": "brass", "tumor_id": "T1",
            "normal_id": "N1", "record_id": rid, "filter": "PASS", "qual": ".", "ins_len": 0,
            "ins_seq": "", "ins_seq_available": True}


def _ins(c, p, rid, seq="", ins_len=60):
    r = _row(c, p, "+", c, p + 1, "-", rid)
    r.update(sv_class="INS", caller="generic", ins_len=ins_len, ins_seq=seq, ins_seq_available=bool(seq))
    return r


@pytest.fixture
def union_dir(tmp_path):
    rows = [_row("1", 100, "+", "1", 900, "-", "a"),       # SV000001.1 DEL (representative)
            _row("1", 130, "+", "1", 920, "-", "e"),       # SV000001.2 DEL: another breakpoint, same cluster
            _row("1", 2000, "+", "2", 50, "-", "b"),       # SV000002.1 BND
            _row("2", 500, "+", "2", 5000, "-", "c"),      # SV000003.1 DEL
            _ins("2", 3000, "f", SEQ),                     # SV000004.1 INS, sequence-resolved
            _ins("2", 3000, "g"),                          # SV000004.2 INS, length only
            _row("2", 8000, "+", "2", 9000, "-", "d")]     # SV000005.1 DEL
    union, adjs, members = sv_union.cluster_union(rows, "T1", ORDER, 100, ["brass"])
    sv_union.write_outputs(union, adjs, members, "T1", str(tmp_path), chunk_size=2, order=ORDER)
    return tmp_path


def _gt(path, records):
    head = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT1\n")
    path.write_text(head + "".join("\t".join(r) + "\n" for r in records))
    return str(path)


def _rec(chrom, pos, rid, sample):
    return [chrom, str(pos), rid, "N", "<DEL>", "10", ".", "SVTYPE=DEL", FMT, sample]


def _sample(gt="0/1", dp=12, ro=8, ao=4, aof=3, rof=7, fs=1, fp=2, fc=0, fx=1):
    return f"{gt}:20:30.5:-3,-1,-2:{dp}:{ro}:{ao}:{ro}:{ao}:4:1:0:4:3:0.33:{fs}:{fp}:{fc}:{aof}:{rof}:{fx}"


def _gt_files(union_dir, tmp_path):
    c0 = _gt(tmp_path / "chunk_00000.gt.vcf", [
        _rec("1", 100, "SV000001.1", _sample()),
        _rec("1", 130, "SV000001.2", _sample(aof=6, rof=2)),
    ])
    c1 = _gt(tmp_path / "chunk_00001.gt.vcf", [
        _rec("1", 2000, "SV000002.1_1", _sample(aof=5, rof=5)),
        _rec("2", 49, "SV000002.1_2", _sample(aof=9, rof=1)),
        _rec("2", 500, "SV000003.1", _sample(gt="0/0", dp=0, ro=0, ao=0, aof=0, rof=0, fs=0, fp=0, fx=0)),
    ])
    c2 = tmp_path / "chunk_00002.gt.vcf"
    c2.write_text("")                                      # SV000005.1 missing from svtyper's output
    return [c0, c1, str(c2)]


def _ins_files(tmp_path, rows=None):
    rows = rows if rows is not None else [
        ("SV000004.1", 6, 7, 1, 0, "genotyped"), ("SV000004.2", 6, 6, 2, 0, "genotyped")]
    p = tmp_path / "ins_counts_ins_chunk_00000.tsv"
    pd.DataFrame(rows, columns=sv_ins_count.COLUMNS).to_csv(p, sep="\t", index=False)
    return [str(p)]


def _finalize(union_dir, gt_files, ins_files, **kw):
    return sv_finalize.finalize(str(union_dir / "T1_union_svs.bedpe.gz"),
                                str(union_dir / "T1_union_sv_adjacencies.tsv.gz"), gt_files,
                                str(union_dir / "chunks.manifest"), "T1", kw.get("stats", []), ins_files,
                                str(union_dir / "ins_chunks.manifest"))


V1_2_COLUMNS = ["sample_id", "sv_id", "chromosome_1", "position_1", "strand_1", "chromosome_2", "position_2",
                "strand_2", "sv_class", "svtyper_type", "sv_length",
                "alt_fragments", "ref_fragments", "alt_split_fragments", "alt_span_fragments",
                "alt_clip_fragments", "conflict_fragments", "frag_vaf",
                "svt_gt", "svt_gq", "svt_sq", "svt_dp", "svt_ro", "svt_ao", "svt_rs", "svt_as", "svt_asc",
                "svt_rp", "svt_ap", "svt_ab", "genotype_status",
                "n_members", "n_callers", "callers", "n_source_samples", "source_samples", "sources",
                "called_in_target", "spread_1", "spread_2", "rep_caller", "rep_source"]


def test_cluster_table_keeps_the_v1_2_columns_first():
    assert sv_finalize.FINAL_COLUMNS[:len(V1_2_COLUMNS)] == V1_2_COLUMNS


def test_cluster_and_adjacency_tables(union_dir, tmp_path):
    clusters, adjs, qc = _finalize(union_dir, _gt_files(union_dir, tmp_path), _ins_files(tmp_path))
    assert list(clusters.columns) == sv_finalize.FINAL_COLUMNS
    assert list(adjs.columns) == sv_finalize.ADJ_FINAL_COLUMNS
    c = clusters.set_index("sv_id")
    # the cluster row is the representative adjacency's genotype (columns as in v1.2) ...
    assert c.loc["SV000001", "alt_fragments"] == 3 and c.loc["SV000001", "ref_fragments"] == 7
    assert c.loc["SV000001", "frag_vaf"] == pytest.approx(0.3)
    assert c.loc["SV000001", "genotype_status"] == "genotyped" and c.loc["SV000001", "rep_adj_id"] == "SV000001.1"
    # ... plus a summary of every adjacency of the cluster
    assert (c.loc["SV000001", "n_adjacencies"], c.loc["SV000001", "n_adjacencies_genotyped"],
            c.loc["SV000001", "n_adjacencies_with_alt"], c.loc["SV000001", "max_alt_fragments"]) == (2, 2, 2, 6)
    assert c.loc["SV000001", "max_frag_vaf"] == pytest.approx(0.75)
    assert (c.loc["SV000002", "alt_fragments"], c.loc["SV000002", "ref_fragments"]) == (5, 5)   # _1 record only
    assert c.loc["SV000003", "genotype_status"] == "no_evidence"
    assert c.loc["SV000005", "genotype_status"] == "missing_from_output"
    assert pd.isna(c.loc["SV000005", "alt_fragments"]) and c.loc["SV000005", "n_adjacencies_genotyped"] == 0
    ins = c.loc["SV000004"]
    assert (ins["genotyper"], ins["sv_class"], ins["alt_fragments"], ins["ref_fragments"]) == ("realign", "INS", 7, 6)
    assert ins["frag_vaf"] == pytest.approx(7 / 13) and pd.isna(ins["svt_gt"]) and pd.isna(ins["alt_split_fragments"])
    a = adjs.set_index("adj_id")
    assert len(a) == 7 and a.loc["SV000001.2", "alt_fragments"] == 6 and not a.loc["SV000001.2", "is_cluster_rep"]
    assert (a.loc["SV000004.2", "alt_fragments"], a.loc["SV000004.2", "uninformative_fragments"]) == (6, 2)
    assert bool(a.loc["SV000004.1", "ins_seq_available"]) and not bool(a.loc["SV000004.2", "ins_seq_available"])
    assert pd.isna(a.loc["SV000001.1", "uninformative_fragments"])
    assert a.loc["SV000005.1", "genotype_status"] == "missing_from_output"
    assert qc["status_counts"] == {"genotyped": 3, "missing_from_output": 1, "no_evidence": 1}
    assert qc["positive_control_fraction"] == pytest.approx(3 / 5)
    assert qc["positive_control_fraction_any_adjacency"] == pytest.approx(3 / 5)
    assert qc["n_adjacencies"] == 7 and qc["clusters_with_multiple_adjacencies"] == 2
    assert qc["adjacency_status_counts"] == {"genotyped": 5, "missing_from_output": 1, "no_evidence": 1}
    assert qc["insertions"] == {"adjacencies": 2, "sequence_resolved": 1, "length_only": 1,
                                "status_counts": {"genotyped": 2}}
    assert qc["clusters_linked_to_other_clusters"] == 0


def test_cluster_alt_support_from_a_non_representative_adjacency(union_dir, tmp_path):
    files = _gt_files(union_dir, tmp_path)
    _gt(tmp_path / "chunk_00000.gt.vcf", [_rec("1", 100, "SV000001.1", _sample(aof=0, rof=9, fs=0, fp=0, fx=0)),
                                          _rec("1", 130, "SV000001.2", _sample(aof=4, rof=5))])
    clusters, _, qc = _finalize(union_dir, files, _ins_files(tmp_path))
    c = clusters.set_index("sv_id").loc["SV000001"]
    assert (c["alt_fragments"], c["n_adjacencies_with_alt"], c["max_alt_fragments"]) == (0, 1, 4)
    assert qc["positive_control_fraction"] == pytest.approx(2 / 5)
    assert qc["positive_control_fraction_any_adjacency"] == pytest.approx(3 / 5)


def test_skipped_insertions_have_no_counts(union_dir, tmp_path):
    ins = _ins_files(tmp_path, [("SV000004.1", None, None, None, None, "skipped_no_reference"),
                                ("SV000004.2", None, None, None, None, "skipped_no_reference")])
    clusters, adjs, qc = _finalize(union_dir, _gt_files(union_dir, tmp_path), ins)
    c = clusters.set_index("sv_id").loc["SV000004"]
    assert c["genotype_status"] == "skipped_no_reference" and pd.isna(c["alt_fragments"])
    assert qc["insertions"]["status_counts"] == {"skipped_no_reference": 2}


def test_not_genotyped_when_svtyper_gives_no_counts(union_dir, tmp_path):
    files = _gt_files(union_dir, tmp_path)
    bad = "./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:."
    _gt(tmp_path / "chunk_00001.gt.vcf", [_rec("2", 500, "SV000003.1", bad)])
    clusters, _, _ = _finalize(union_dir, files, _ins_files(tmp_path))
    c = clusters.set_index("sv_id")
    assert c.loc["SV000003", "genotype_status"] == "not_genotyped"
    assert c.loc["SV000002", "genotype_status"] == "missing_from_output"


def test_missing_chunks_fail(union_dir, tmp_path):
    files = _gt_files(union_dir, tmp_path)
    with pytest.raises(ValueError, match="missing"):
        _finalize(union_dir, files[:1], _ins_files(tmp_path))
    with pytest.raises(ValueError, match="missing"):
        _finalize(union_dir, files, [])


def test_insertion_counts_must_cover_every_insertion(union_dir, tmp_path):
    files = _gt_files(union_dir, tmp_path)
    with pytest.raises(ValueError, match="insertion"):
        _finalize(union_dir, files, _ins_files(tmp_path, [("SV000004.1", 6, 7, 1, 0, "genotyped")]))
    with pytest.raises(ValueError, match="insertion"):
        _finalize(union_dir, files, _ins_files(tmp_path, [("SV000004.1", 6, 7, 1, 0, "genotyped"),
                                                          ("SV000004.1", 6, 7, 1, 0, "genotyped"),
                                                          ("SV000004.2", 6, 7, 1, 0, "genotyped")]))


def _cli(union_dir, gt_files, ins_files, out):
    sv_finalize.main(["--union", str(union_dir / "T1_union_svs.bedpe.gz"),
                      "--adjacencies", str(union_dir / "T1_union_sv_adjacencies.tsv.gz"),
                      "--genotyped", *gt_files, "--manifest", str(union_dir / "chunks.manifest"),
                      "--ins-counts", *ins_files, "--ins-manifest", str(union_dir / "ins_chunks.manifest"),
                      "--sample-id", "T1", "--contigs", "1,2", "--outdir", str(out)])


def test_cli_writes_tables_yaml_qc_and_sorted_vcf(union_dir, tmp_path):
    out = tmp_path / "final"
    out.mkdir()
    _cli(union_dir, _gt_files(union_dir, tmp_path), _ins_files(tmp_path), out)
    for name, cols in (("T1_sv_genotyping.csv.gz", sv_finalize.FINAL_COLUMNS),
                       ("T1_sv_adjacency_genotyping.csv.gz", sv_finalize.ADJ_FINAL_COLUMNS)):
        table = pd.read_csv(out / name)
        assert list(table.columns) == cols
        meta = yaml.safe_load(open(out / f"{name}.yaml"))
        assert [c["name"] for c in meta["columns"]] == cols
    assert len(pd.read_csv(out / "T1_sv_genotyping.csv.gz")) == 5
    assert len(pd.read_csv(out / "T1_sv_adjacency_genotyping.csv.gz")) == 7
    assert json.load(open(out / "T1_genotyping_qc.json"))["type"] == "sv_genotyping"
    with gzip.open(out / "T1_svtyper.vcf.gz", "rt") as fh:
        recs = [l.split("\t")[:3] for l in fh if not l.startswith("#")]
    assert recs == [["1", "100", "SV000001.1"], ["1", "130", "SV000001.2"], ["1", "2000", "SV000002.1_1"],
                    ["2", "49", "SV000002.1_2"], ["2", "500", "SV000003.1"]]
    assert (out / "T1_svtyper.vcf.gz.tbi").exists()


def test_cli_with_all_empty_outputs_writes_header_only_tables(tmp_path):
    union, adjs, members = sv_union.cluster_union([], "T1", ORDER, 100, ["brass"])
    sv_union.write_outputs(union, adjs, members, "T1", str(tmp_path), chunk_size=5, order=ORDER)
    (tmp_path / "chunk_00000.gt.vcf").write_text("")
    ins = tmp_path / "ins_counts_ins_chunk_00000.tsv"
    ins.write_text("\t".join(sv_ins_count.COLUMNS) + "\n")          # sv_ins_count on the empty chunk
    out = tmp_path / "final"
    out.mkdir()
    _cli(tmp_path, [str(tmp_path / "chunk_00000.gt.vcf")], [str(ins)], out)
    for name, cols in (("T1_sv_genotyping.csv.gz", sv_finalize.FINAL_COLUMNS),
                       ("T1_sv_adjacency_genotyping.csv.gz", sv_finalize.ADJ_FINAL_COLUMNS)):
        table = pd.read_csv(out / name)
        assert len(table) == 0 and list(table.columns) == cols
    qc = json.load(open(out / "T1_genotyping_qc.json"))
    assert qc["positive_control_fraction"] is None and qc["n_adjacencies"] == 0


# --- BEDPE (one row per adjacency; 0-based half-open junction bases) --------------------------

BEDPE_HEAD = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"]


def test_adjacency_table_carries_the_union_ci(tmp_path):
    imprecise = _row("1", 100, "+", "1", 900, "-", "a")
    imprecise.update(ci1_lo=-20, ci1_hi=30, ci2_lo=-5, ci2_hi=5, imprecise=True)
    union, adjs, members = sv_union.cluster_union([imprecise], "T1", ORDER, 100, ["brass"])
    sv_union.write_outputs(union, adjs, members, "T1", str(tmp_path), chunk_size=5, order=ORDER)
    gt = _gt(tmp_path / "chunk_00000.gt.vcf", [_rec("1", 100, "SV000001.1", _sample())])
    ins = tmp_path / "ins_counts_ins_chunk_00000.tsv"
    ins.write_text("\t".join(sv_ins_count.COLUMNS) + "\n")
    _, a, _ = _finalize(tmp_path, [gt], [str(ins)])
    assert ["ci1_lo", "ci1_hi", "ci2_lo", "ci2_hi"] == [c for c in sv_finalize.ADJ_FINAL_COLUMNS if c.startswith("ci")]
    assert a.iloc[0][["ci1_lo", "ci1_hi", "ci2_lo", "ci2_hi"]].tolist() == [-20, 30, -5, 5]
    b = sv_finalize.bedpe_table(a, ORDER).iloc[0]
    # the interval is the junction base only; the CI is reported in its own columns
    assert (b["start1"], b["end1"], b["start2"], b["end2"]) == (99, 100, 899, 900)
    assert (b["ci1_lo"], b["ci1_hi"], b["ci2_lo"], b["ci2_hi"]) == (-20, 30, -5, 5)


def test_bedpe_table_is_the_adjacency_table_in_bedpe_layout(union_dir, tmp_path):
    _, adjs, _ = _finalize(union_dir, _gt_files(union_dir, tmp_path), _ins_files(tmp_path))
    b = sv_finalize.bedpe_table(adjs, ORDER)
    assert list(b.columns) == sv_finalize.BEDPE_COLUMNS
    assert sv_finalize.BEDPE_COLUMNS[:10] == BEDPE_HEAD
    extra = sv_finalize.BEDPE_COLUMNS[10:]
    coords = {"adj_id", "chromosome_1", "position_1", "strand_1", "chromosome_2", "position_2", "strand_2"}
    assert extra == [c for c in sv_finalize.ADJ_FINAL_COLUMNS if c not in coords]
    assert len(b) == len(adjs) == 7
    a = adjs.set_index("adj_id")
    for r in b.to_dict("records"):
        s = a.loc[r["name"]]
        assert (r["chrom1"], r["start1"], r["end1"], r["strand1"]) == (
            s["chromosome_1"], s["position_1"] - 1, s["position_1"], s["strand_1"])
        assert (r["chrom2"], r["start2"], r["end2"], r["strand2"]) == (
            s["chromosome_2"], s["position_2"] - 1, s["position_2"], s["strand_2"])
        for c in extra:
            assert (pd.isna(r[c]) and pd.isna(s[c])) or r[c] == s[c], (r["name"], c)
    by = b.set_index("name")
    # translocation: one base on each contig
    assert by.loc["SV000002.1", ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]].tolist() == \
        ["1", 1999, 2000, "2", 49, 50]
    # insertion: the bases left and right of the inserted sequence
    assert by.loc["SV000004.1", ["start1", "end1", "start2", "end2", "strand1", "strand2"]].tolist() == \
        [2999, 3000, 3000, 3001, "+", "-"]
    # score = alt_fragments (missing when not counted)
    assert by.loc["SV000001.2", "score"] == 6 and pd.isna(by.loc["SV000005.1", "score"])


def test_bedpe_table_sorts_in_reference_order(union_dir, tmp_path):
    _, adjs, _ = _finalize(union_dir, _gt_files(union_dir, tmp_path), _ins_files(tmp_path))
    names = sv_finalize.bedpe_table(adjs, {"2": 0, "1": 1})["name"].tolist()
    assert names == ["SV000003.1", "SV000004.1", "SV000004.2", "SV000005.1",      # contig 2 first, by start
                     "SV000001.1", "SV000001.2", "SV000002.1"]
    shuffled = adjs.sample(frac=1, random_state=3)
    assert sv_finalize.bedpe_table(shuffled, {"2": 0, "1": 1})["name"].tolist() == names


def test_write_bedpe_header_missing_values_and_round_trip(union_dir, tmp_path):
    _, adjs, _ = _finalize(union_dir, _gt_files(union_dir, tmp_path), _ins_files(tmp_path))
    b = sv_finalize.bedpe_table(adjs, ORDER)
    path = tmp_path / "T1.bedpe"
    sv_finalize.write_bedpe(b, str(path))
    lines = path.read_text().splitlines()
    assert lines[0] == "#" + "\t".join(sv_finalize.BEDPE_COLUMNS)
    assert len(lines) == 8 and all(len(l.split("\t")) == len(sv_finalize.BEDPE_COLUMNS) for l in lines)
    rec = dict(zip(sv_finalize.BEDPE_COLUMNS, next(l for l in lines if "\tSV000005.1\t" in l).split("\t")))
    assert (rec["score"], rec["alt_fragments"], rec["svt_gt"], rec["frag_vaf"]) == (".", ".", ".", ".")
    assert not any(tok in l.split("\t") for l in lines for tok in ("nan", "NaN", "<NA>", "None", ""))
    back = pd.read_csv(path, sep="\t", na_values=".", keep_default_na=False,
                       dtype={"#chrom1": str, "chrom2": str})
    assert back["name"].tolist() == b["name"].tolist()
    assert back[["start1", "end1", "start2", "end2"]].values.tolist() == \
        b[["start1", "end1", "start2", "end2"]].values.tolist()


def test_cli_writes_bedpe(union_dir, tmp_path):
    out = tmp_path / "final"
    out.mkdir()
    _cli(union_dir, _gt_files(union_dir, tmp_path), _ins_files(tmp_path), out)
    lines = (out / "T1_sv_adjacency_genotyping.bedpe").read_text().splitlines()
    assert lines[0].startswith("#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\t")
    assert [l.split("\t")[6] for l in lines[1:]] == ["SV000001.1", "SV000001.2", "SV000002.1", "SV000003.1",
                                                     "SV000004.1", "SV000004.2", "SV000005.1"]


def test_cli_with_empty_union_writes_header_only_bedpe(tmp_path):
    union, adjs, members = sv_union.cluster_union([], "T1", ORDER, 100, ["brass"])
    sv_union.write_outputs(union, adjs, members, "T1", str(tmp_path), chunk_size=5, order=ORDER)
    (tmp_path / "chunk_00000.gt.vcf").write_text("")
    ins = tmp_path / "ins_counts_ins_chunk_00000.tsv"
    ins.write_text("\t".join(sv_ins_count.COLUMNS) + "\n")
    out = tmp_path / "final"
    out.mkdir()
    _cli(tmp_path, [str(tmp_path / "chunk_00000.gt.vcf")], [str(ins)], out)
    assert (out / "T1_sv_adjacency_genotyping.bedpe").read_text() == \
        "#" + "\t".join(sv_finalize.BEDPE_COLUMNS) + "\n"
