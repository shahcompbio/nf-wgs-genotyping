import json

import pytest

import sv_normalize

ORDER = {"1": 0, "2": 1, "3": 2}


def _vcf(tmp_path, lines, samples=("NORMAL", "TUMOUR"), name="in.vcf"):
    p = tmp_path / name
    head = "##fileformat=VCFv4.2\n" + "\t".join(
        ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]) + "\n"
    p.write_text(head + "".join("\t".join(map(str, l)) + "\n" for l in lines))
    return str(p)


def _norm(path, caller, **kw):
    kw.setdefault("order", ORDER)
    return sv_normalize.normalize(path, caller, "T1", "N1", **kw)


def _ends(rows):
    return [(r["chrom1"], r["pos1"], r["strand1"], r["chrom2"], r["pos2"], r["strand2"]) for r in rows]


def test_brass_missing_filter_is_pass_and_mates_collapse(tmp_path):
    vcf = _vcf(tmp_path, [
        ("1", 1000, "b_1", "N", "N[1:5000[", "8", ".", "SVTYPE=BND;MATEID=b_2", "GT", "0/0", "0/1"),
        ("1", 5000, "b_2", "N", "]1:1000]N", "8", ".", "SVTYPE=BND;MATEID=b_1", "GT", "0/0", "0/1"),
    ])
    rows, stats = _norm(vcf, "brass")
    assert _ends(rows) == [("1", 1000, "+", "1", 5000, "-")]
    assert rows[0]["sv_class"] == "DEL" and rows[0]["caller"] == "brass"
    assert stats["dropped"]["mate_duplicate"] == 1
    rows, stats = _norm(vcf, "brass", brass_missing_filter_is_pass=False)
    assert rows == [] and stats["dropped"]["not_pass"] == 2


def test_svaba_requires_pass_and_recovers_missing_mate_from_alt(tmp_path):
    vcf = _vcf(tmp_path, [
        ("2", 300, "s:1", "C", "]1:900]C", "55", "PASS", "SVTYPE=BND;MATEID=s:2;IMPRECISE", "GT", "0/0", "0/1"),
        ("1", 50, "x:1", "C", "C[1:4000[", "55", "LowQual", "SVTYPE=BND;MATEID=x:2", "GT", "0/0", "0/1"),
    ], samples=("/a/n.bam", "/a/t.bam"))
    rows, stats = _norm(vcf, "svaba", imprecise_ci=150)
    # ends are put in canonical order (chrom 1 before chrom 2); CI travels with its end
    assert _ends(rows) == [("1", 900, "+", "2", 300, "-")]
    r = rows[0]
    assert r["sv_class"] == "TRA" and r["imprecise"] is True
    assert (r["ci1_lo"], r["ci1_hi"], r["ci2_lo"], r["ci2_hi"]) == (-150, 150, -150, 150)
    assert stats["dropped"]["not_pass"] == 1


def test_gridss_somatic_filter_on_normal_vf_and_single_breakends(tmp_path):
    lines = [
        ("1", 100, "g1o", "A", "A[1:2000[", "500", "PASS", "PARID=g1h;CIPOS=-3,4;CIRPOS=-5,6", "GT:VF", ".:0", ".:20"),
        ("1", 2000, "g1h", "A", "]1:100]A", "500", "PASS", "PARID=g1o;CIPOS=-5,6;CIRPOS=-3,4", "GT:VF", ".:0", ".:20"),
        ("1", 3000, "g2o", "A", "A]2:500]", "900", "PASS", "PARID=g2h", "GT:VF", ".:25", ".:20"),
        ("1", 3500, "g3b", "A", ".AACCGGTT", "900", "PASS", "", "GT:VF", ".:0", ".:20"),
        ("1", 3600, "g4o", "A", "A[1:9000[", "50", "LOW_QUAL", "PARID=g4h", "GT:VF", ".:0", ".:3"),
    ]
    vcf = _vcf(tmp_path, lines, samples=("N1", "T1"))
    rows, stats = _norm(vcf, "gridss")
    assert _ends(rows) == [("1", 100, "+", "1", 2000, "-")]
    assert (rows[0]["ci1_lo"], rows[0]["ci1_hi"], rows[0]["ci2_lo"], rows[0]["ci2_hi"]) == (-3, 4, -5, 6)
    assert stats["dropped"]["germline_normal_support"] == 1
    assert stats["dropped"]["single_breakend"] == 1
    assert stats["dropped"]["not_pass"] == 1
    rows, _ = _norm(vcf, "gridss", gridss_max_normal_vf=100)
    assert len(rows) == 2


def test_gridss_without_normal_column_is_an_error(tmp_path):
    vcf = _vcf(tmp_path, [("1", 100, "g", "A", "A[1:2000[", "5", "PASS", "PARID=h", "GT:VF", ".:0", ".:2")],
               samples=("OTHER", "T1"))
    with pytest.raises(ValueError, match="N1"):
        _norm(vcf, "gridss")


def test_symbolic_alleles_expand_to_junctions(tmp_path):
    vcf = _vcf(tmp_path, [
        ("1", 100, "d", "N", "<DEL>", ".", "PASS", "SVTYPE=DEL;END=500;CIPOS=-1,1;CIEND=-2,2", "GT", "0/0", "0/1"),
        ("2", 100, "u", "N", "<DUP:TANDEM>", ".", "PASS", "SVTYPE=DUP;END=500", "GT", "0/0", "0/1"),
        ("3", 100, "v", "N", "<INV>", ".", "PASS", "SVTYPE=INV;END=500", "GT", "0/0", "0/1"),
        ("3", 900, "i", "N", "<INS>", ".", "PASS", "SVTYPE=INS;END=900", "GT", "0/0", "0/1"),
    ])
    rows, stats = _norm(vcf, "generic")
    assert _ends(rows) == [
        ("1", 100, "+", "1", 501, "-"),
        ("2", 101, "-", "2", 500, "+"),
        ("3", 100, "+", "3", 500, "+"),
        ("3", 101, "-", "3", 501, "-"),
    ]
    assert [r["sv_class"] for r in rows] == ["DEL", "DUP", "INV", "INV"]
    assert (rows[0]["ci2_lo"], rows[0]["ci2_hi"]) == (-2, 2)
    assert stats["dropped"]["unsupported_allele"] == 1


def test_small_and_off_contig_events_are_dropped(tmp_path):
    vcf = _vcf(tmp_path, [
        ("1", 100, "a", "N", "N[1:130[", ".", "PASS", "SVTYPE=BND", "GT", "0/0", "0/1"),
        ("1", 100, "b", "N", "N[GL1:130[", ".", "PASS", "SVTYPE=BND", "GT", "0/0", "0/1"),
        ("1", 100, "c", "N", "N[1:1000[", ".", "PASS", "SVTYPE=BND", "GT", "0/0", "0/1"),
    ])
    rows, stats = _norm(vcf, "generic", contigs=["1", "2", "3"], min_size=50)
    assert _ends(rows) == [("1", 100, "+", "1", 1000, "-")]
    assert stats["dropped"]["small"] == 1 and stats["dropped"]["off_contig"] == 1


def test_synthetic_fixtures(synth):
    files = {k: v for k, v in synth["sv_vcfs"].items()}
    rows, stats = sv_normalize.normalize(files[("T", "gridss")], "gridss", "NA12878", "N1",
                                         order={"2": 0})
    assert len(rows) == 1 and stats["dropped"]["germline_normal_support"] == 2   # both mate records
    rows, _ = sv_normalize.normalize(files[("T", "brass")], "brass", "NA12878", "N1", order={"2": 0})
    assert len(rows) == 2


def test_cli_writes_bedpe_and_stats(synth, tmp_path):
    out, st = tmp_path / "o.tsv", tmp_path / "o.json"
    sv_normalize.main(["--vcf", synth["sv_vcfs"][("T", "svaba")], "--caller", "svaba", "--tumor-id", "NA12878",
                       "--normal-id", "N1", "--contigs", "1,2,3", "--out-tsv", str(out), "--out-stats", str(st)])
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == sv_normalize.COLUMNS
    assert len(lines) == 3
    assert json.loads(st.read_text())["kept"] == 2
