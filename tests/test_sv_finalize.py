import gzip
import json

import pandas as pd
import pytest
import yaml

import sv_finalize
import sv_normalize
import sv_union

ORDER = {"1": 0, "2": 1}
FMT = "GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:FS:FP:FC:AOF:ROF:FX"


def _row(c1, p1, s1, c2, p2, s2, rid):
    return {"chrom1": c1, "pos1": p1, "strand1": s1, "chrom2": c2, "pos2": p2, "strand2": s2,
            "ci1_lo": 0, "ci1_hi": 0, "ci2_lo": 0, "ci2_hi": 0, "imprecise": False,
            "sv_class": sv_normalize.sv_class(c1, s1, c2, s2), "caller": "brass", "tumor_id": "T1",
            "normal_id": "N1", "record_id": rid, "filter": "PASS", "qual": ".", "ins_len": 0}


@pytest.fixture
def union_dir(tmp_path):
    rows = [_row("1", 100, "+", "1", 900, "-", "a"),       # SV000001 DEL
            _row("1", 2000, "+", "2", 50, "-", "b"),       # SV000002 BND
            _row("2", 500, "+", "2", 5000, "-", "c"),      # SV000003 DEL
            _row("2", 8000, "+", "2", 9000, "-", "d")]     # SV000004 DEL
    union, members = sv_union.cluster_union(rows, "T1", ORDER, 100, ["brass"])
    sv_union.write_outputs(union, members, "T1", str(tmp_path), chunk_size=2, order=ORDER)
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
        _rec("1", 100, "SV000001", _sample()),
        _rec("1", 2000, "SV000002_1", _sample(aof=5, rof=5)),
        _rec("2", 49, "SV000002_2", _sample(aof=9, rof=1)),
    ])
    c1 = _gt(tmp_path / "chunk_00001.gt.vcf", [
        _rec("2", 500, "SV000003", _sample(gt="0/0", dp=0, ro=0, ao=0, aof=0, rof=0, fs=0, fp=0, fx=0)),
    ])
    return [c0, c1]


def test_parse_uses_first_bnd_record_and_statuses(union_dir, tmp_path):
    files = _gt_files(union_dir, tmp_path)
    df, qc = sv_finalize.finalize(str(union_dir / "T1_union_svs.bedpe.gz"), files,
                                  str(union_dir / "chunks.manifest"), "T1", [])
    df = df.set_index("sv_id")
    assert df.loc["SV000001", "alt_fragments"] == 3 and df.loc["SV000001", "ref_fragments"] == 7
    assert df.loc["SV000001", "frag_vaf"] == pytest.approx(0.3)
    assert df.loc["SV000001", "genotype_status"] == "genotyped"
    assert (df.loc["SV000002", "alt_fragments"], df.loc["SV000002", "ref_fragments"]) == (5, 5)
    assert df.loc["SV000003", "genotype_status"] == "no_evidence"
    assert df.loc["SV000004", "genotype_status"] == "missing_from_output"
    assert pd.isna(df.loc["SV000004", "alt_fragments"])
    assert qc["status_counts"] == {"genotyped": 2, "missing_from_output": 1, "no_evidence": 1}
    assert qc["positive_control_fraction"] == pytest.approx(2 / 4)


def test_not_genotyped_when_svtyper_gives_no_counts(union_dir, tmp_path):
    files = _gt_files(union_dir, tmp_path)
    bad = "./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:."
    _gt(tmp_path / "chunk_00001.gt.vcf", [_rec("2", 500, "SV000003", bad)])
    df, _ = sv_finalize.finalize(str(union_dir / "T1_union_svs.bedpe.gz"), files,
                                 str(union_dir / "chunks.manifest"), "T1", [])
    assert df.set_index("sv_id").loc["SV000003", "genotype_status"] == "not_genotyped"


def test_empty_svtyper_output_file_is_allowed(union_dir, tmp_path):
    files = _gt_files(union_dir, tmp_path)
    (tmp_path / "chunk_00001.gt.vcf").write_text("")
    df, _ = sv_finalize.finalize(str(union_dir / "T1_union_svs.bedpe.gz"), files,
                                 str(union_dir / "chunks.manifest"), "T1", [])
    assert df.set_index("sv_id").loc["SV000003", "genotype_status"] == "missing_from_output"


def test_missing_chunk_fails(union_dir, tmp_path):
    files = _gt_files(union_dir, tmp_path)
    with pytest.raises(ValueError, match="missing"):
        sv_finalize.finalize(str(union_dir / "T1_union_svs.bedpe.gz"), files[:1],
                             str(union_dir / "chunks.manifest"), "T1", [])


def test_cli_writes_table_yaml_qc_and_sorted_vcf(union_dir, tmp_path):
    files = _gt_files(union_dir, tmp_path)
    out = tmp_path / "final"
    out.mkdir()
    sv_finalize.main(["--union", str(union_dir / "T1_union_svs.bedpe.gz"), "--genotyped", *files,
                      "--manifest", str(union_dir / "chunks.manifest"), "--sample-id", "T1",
                      "--contigs", "1,2", "--outdir", str(out)])
    table = pd.read_csv(out / "T1_sv_genotyping.csv.gz")
    assert list(table.columns) == sv_finalize.FINAL_COLUMNS
    meta = yaml.safe_load(open(out / "T1_sv_genotyping.csv.gz.yaml"))
    assert [c["name"] for c in meta["columns"]] == sv_finalize.FINAL_COLUMNS
    assert json.load(open(out / "T1_genotyping_qc.json"))["type"] == "sv_genotyping"
    with gzip.open(out / "T1_svtyper.vcf.gz", "rt") as fh:
        recs = [l.split("\t")[:2] for l in fh if not l.startswith("#")]
    assert recs == [["1", "100"], ["1", "2000"], ["2", "49"], ["2", "500"]]
    assert (out / "T1_svtyper.vcf.gz.tbi").exists()


def test_cli_with_all_empty_outputs_writes_header_only_table(tmp_path):
    union, members = sv_union.cluster_union([], "T1", ORDER, 100, ["brass"])
    sv_union.write_outputs(union, members, "T1", str(tmp_path), chunk_size=5, order=ORDER)
    (tmp_path / "chunk_00000.gt.vcf").write_text("")
    out = tmp_path / "final"
    out.mkdir()
    sv_finalize.main(["--union", str(tmp_path / "T1_union_svs.bedpe.gz"), "--genotyped",
                      str(tmp_path / "chunk_00000.gt.vcf"), "--manifest", str(tmp_path / "chunks.manifest"),
                      "--sample-id", "T1", "--outdir", str(out)])
    table = pd.read_csv(out / "T1_sv_genotyping.csv.gz")
    assert len(table) == 0 and list(table.columns) == sv_finalize.FINAL_COLUMNS
