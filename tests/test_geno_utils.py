import gzip
import json

import pandas as pd
import pytest
import yaml

import geno_utils as gu


def test_contig_order_from_fai_keeps_file_order(tmp_path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("2\t100\t3\t60\t61\n1\t100\t110\t60\t61\nX\t50\t220\t60\t61\n")
    assert gu.contig_order(fai=str(fai)) == {"2": 0, "1": 1, "X": 2}


def test_contig_order_from_list_when_no_fai():
    assert gu.contig_order(contigs=["1", "2", "X"]) == {"1": 0, "2": 1, "X": 2}


def test_parse_contigs_handles_empty_and_spaces():
    assert gu.parse_contigs("1, 2,X") == ["1", "2", "X"]
    assert gu.parse_contigs("") is None
    assert gu.parse_contigs(None) is None


def test_site_sort_key_puts_unknown_contigs_last():
    key = gu.site_sort_key({"1": 0, "2": 1})
    sites = [("2", 5), ("GL1", 1), ("1", 10), ("1", 2)]
    assert sorted(sites, key=lambda s: key(*s)) == [("1", 2), ("1", 10), ("2", 5), ("GL1", 1)]


@pytest.mark.parametrize(
    "alt, expected",
    [
        ("G[22:200[", ("+", "22", 200, "-", "")),
        ("G]22:200]", ("+", "22", 200, "+", "")),
        ("]22:200]C", ("-", "22", 200, "+", "")),
        ("[22:200[C", ("-", "22", 200, "-", "")),
        ("CCCG[7:55[", ("+", "7", 55, "-", "CCG")),
        ("]X:9]AGGT", ("-", "X", 9, "+", "AGG")),
    ],
)
def test_parse_bnd_alt(alt, expected):
    p = gu.parse_bnd_alt(alt)
    assert (p["strand1"], p["chrom2"], p["pos2"], p["strand2"], p["ins_seq"]) == expected


@pytest.mark.parametrize("alt", ["<DEL>", "A", ".ACGT", "ACGT.", "G[22:200]", "N"])
def test_parse_bnd_alt_rejects_non_bnd(alt):
    assert gu.parse_bnd_alt(alt) is None


def test_single_breakend_detection():
    assert gu.is_single_breakend(".AACG")
    assert gu.is_single_breakend("AACG.")
    assert not gu.is_single_breakend("A[1:5[")
    assert not gu.is_single_breakend("<DEL>")


@pytest.mark.parametrize(
    "strand1, strand2, expected",
    [("+", "-", "N[3:10["), ("+", "+", "N]3:10]"), ("-", "+", "]3:10]N"), ("-", "-", "[3:10[N")],
)
def test_format_bnd_alt_round_trips(strand1, strand2, expected):
    alt = gu.format_bnd_alt("N", strand1, "3", 10, strand2)
    assert alt == expected
    p = gu.parse_bnd_alt(alt)
    assert (p["strand1"], p["strand2"], p["chrom2"], p["pos2"]) == (strand1, strand2, "3", 10)


def test_write_table_writes_gz_csv_and_csverve_yaml(tmp_path):
    df = pd.DataFrame({"chromosome": ["1"], "position": [5], "vaf": [0.5]})
    out = tmp_path / "x.csv.gz"
    gu.write_table(df, str(out), {"chromosome": "category", "position": "int", "vaf": "float"})
    with gzip.open(out, "rt") as fh:
        assert fh.readline().strip() == "chromosome,position,vaf"
    meta = yaml.safe_load(open(str(out) + ".yaml"))
    assert meta == {
        "header": True,
        "sep": ",",
        "columns": [
            {"name": "chromosome", "dtype": "category"},
            {"name": "position", "dtype": "int"},
            {"name": "vaf", "dtype": "float"},
        ],
    }


def test_write_table_rejects_missing_dtype(tmp_path):
    df = pd.DataFrame({"a": [1], "b": [2]})
    with pytest.raises(ValueError, match="b"):
        gu.write_table(df, str(tmp_path / "x.csv.gz"), {"a": "int"})


def test_read_samplesheet_validates_columns(tmp_path):
    ss = tmp_path / "ss.csv"
    ss.write_text("tumor_id,normal_id,caller,vcf\nT1,N1,mutect,/a.vcf\n")
    rows = gu.read_samplesheet(str(ss))
    assert rows == [{"tumor_id": "T1", "normal_id": "N1", "caller": "mutect", "vcf": "/a.vcf"}]
    bad = tmp_path / "bad.csv"
    bad.write_text("tumor_id,caller,vcf\nT1,mutect,/a.vcf\n")
    with pytest.raises(ValueError, match="normal_id"):
        gu.read_samplesheet(str(bad))


def test_manifest_round_trip_and_gather_check(tmp_path):
    m = tmp_path / "chunks.manifest"
    gu.write_manifest(str(m), [("chunk_00000", 3), ("chunk_00001", 1)])
    assert gu.read_manifest(str(m)) == {"chunk_00000": 3, "chunk_00001": 1}
    gu.check_gather(str(m), ["chunk_00001", "chunk_00000"])
    with pytest.raises(ValueError, match="missing.*chunk_00001"):
        gu.check_gather(str(m), ["chunk_00000"])
    with pytest.raises(ValueError, match="unexpected.*chunk_00009"):
        gu.check_gather(str(m), ["chunk_00000", "chunk_00001", "chunk_00009"])


def test_open_text_reads_plain_and_gzip(tmp_path):
    p = tmp_path / "a.txt"
    p.write_text("hello\n")
    g = tmp_path / "a.txt.gz"
    with gzip.open(g, "wt") as fh:
        fh.write("hello\n")
    for path in (p, g):
        with gu.open_text(str(path)) as fh:
            assert fh.read() == "hello\n"


def test_write_json_is_sorted_and_stable(tmp_path):
    p = tmp_path / "a.json"
    gu.write_json({"b": 1, "a": {"d": 2, "c": 3}}, str(p))
    assert json.loads(p.read_text()) == {"a": {"c": 3, "d": 2}, "b": 1}
    assert p.read_text().index('"a"') < p.read_text().index('"b"')


def test_iter_pass_alleles_filters_records_and_splits_alts(tmp_path):
    from collections import Counter

    vcf = tmp_path / "x.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "1\t10\t.\tac\tA,acg\t.\tPASS\t.\n"
                   "1\t20\t.\tA\tG\t.\tLowQual\t.\n"
                   "1\t30\t.\tA\tG\t.\t.\t.\n"
                   "MT\t5\t.\tA\tG\t.\tPASS\t.\n\n")
    tally = Counter()
    got = list(gu.iter_pass_alleles(str(vcf), tally, contigs=["1"]))
    assert got == [("1", 10, "AC", "A"), ("1", 10, "AC", "ACG")]
    assert (tally["records"], tally["not_pass"], tally["off_contig"]) == (4, 2, 1)
    tally = Counter()
    got = list(gu.iter_pass_alleles(str(vcf), tally, missing_filter_is_pass=True))
    assert [g[1] for g in got] == [10, 10, 30, 5]
    assert tally["not_pass"] == 1 and tally["off_contig"] == 0
