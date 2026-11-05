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
        ("3", 950, "x", "N", "<CNV>", ".", "PASS", "SVTYPE=CNV;END=2000", "GT", "0/0", "0/1"),
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
    # an <INS> with neither a length nor a sequence cannot be genotyped or size-filtered
    assert stats["dropped"]["insertion_length_unknown"] == 1


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


# ----------------------------------------------------------------------------- insertions
SEQ60 = "GATTACAGGCTTCAAGCTAGCATCGGATCCTTAGGCTAACGTTGCAATGCGTACGTCAGT"
assert len(SEQ60) == 60


def _ins(rows):
    return [(r["chrom1"], r["pos1"], r["pos2"], r["sv_class"], r["ins_len"], r["ins_seq"], r["ins_seq_available"])
            for r in rows]


def test_symbolic_insertions_keep_length_and_sequence(tmp_path):
    vcf = _vcf(tmp_path, [
        ("1", 900, "a", "N", "<INS>", ".", "PASS", "SVTYPE=INS;END=900;SVLEN=300", "GT", "0/0", "0/1"),
        ("1", 1900, "b", "N", "<INS:ME:ALU>", ".", "PASS", "SVTYPE=INS;SVLEN=-310", "GT", "0/0", "0/1"),
        ("2", 500, "c", "N", "<INS>", ".", "PASS", f"SVTYPE=INS;END=500;SVINSSEQ={SEQ60.lower()}", "GT", "0/0", "0/1"),
        ("2", 700, "d", "N", "<INS>", ".", "PASS", "SVTYPE=INS;INSLEN=80", "GT", "0/0", "0/1"),
        ("2", 800, "e", "N", "<INS>", ".", "PASS", "SVTYPE=INS;SVLEN=49", "GT", "0/0", "0/1"),
        ("2", 850, "f", "N", "<INS>", ".", "PASS", "SVTYPE=INS;SVLEN=50", "GT", "0/0", "0/1"),
    ])
    rows, stats = _norm(vcf, "generic", min_size=50)
    assert _ins(rows) == [
        ("1", 900, 901, "INS", 300, "", False),
        ("1", 1900, 1901, "INS", 310, "", False),
        ("2", 500, 501, "INS", 60, SEQ60, True),
        ("2", 700, 701, "INS", 80, "", False),
        ("2", 850, 851, "INS", 50, "", False),
    ]
    assert all((r["strand1"], r["strand2"]) == ("+", "-") for r in rows)
    assert stats["dropped"]["small"] == 1          # the 49-bp insertion: indel workflow territory
    assert stats["insertions"] == {"kept": 5, "sequence_resolved": 1, "length_only": 4, "left_aligned": 0}


def test_imprecise_insertion_ci_applies_to_both_ends(tmp_path):
    vcf = _vcf(tmp_path, [
        ("1", 900, "a", "N", "<INS>", ".", "PASS", "SVTYPE=INS;SVLEN=300;IMPRECISE", "GT", "0/0", "0/1"),
        ("1", 1900, "b", "N", "<INS>", ".", "PASS", "SVTYPE=INS;SVLEN=300;IMPRECISE;CIPOS=-20,30", "GT", "0/0", "0/1"),
    ])
    rows, _ = _norm(vcf, "generic", imprecise_ci=150)
    assert [(r["imprecise"], r["ci1_lo"], r["ci1_hi"], r["ci2_lo"], r["ci2_hi"]) for r in rows] == [
        (True, -150, 150, -150, 150), (True, -20, 30, -20, 30)]


def test_insertion_breakend_mates_collapse_to_one_insertion(tmp_path):
    vcf = _vcf(tmp_path, [
        ("1", 1000, "g1o", "A", f"A{SEQ60}[1:1001[", "500", "PASS", "SVTYPE=BND;PARID=g1h", "GT:VF", ".:0", ".:20"),
        ("1", 1001, "g1h", "C", f"]1:1000]{SEQ60}C", "500", "PASS", "SVTYPE=BND;PARID=g1o", "GT:VF", ".:0", ".:20"),
        # 5 deleted bases replaced by 80 inserted bases: an insertion (net +75) with del_len 5
        ("1", 2000, "g2o", "A", "A" + "T" * 80 + "[1:2006[", "500", "PASS", "SVTYPE=BND", "GT:VF", ".:0", ".:20"),
        # 60 deleted bases: a deletion junction that also carries inserted sequence
        ("1", 3000, "g3o", "A", "A" + "G" * 100 + "[1:3061[", "500", "PASS", "SVTYPE=BND", "GT:VF", ".:0", ".:20"),
        # 10 inserted bases: too small for the SV pipeline
        ("1", 4000, "g4o", "A", "A" + "G" * 10 + "[1:4001[", "500", "PASS", "SVTYPE=BND", "GT:VF", ".:0", ".:20"),
    ], samples=("N1", "T1"))
    rows, stats = _norm(vcf, "gridss")
    assert _ins(rows) == [
        ("1", 1000, 1001, "INS", 60, SEQ60, True),
        ("1", 2000, 2006, "INS", 80, "T" * 80, True),
        ("1", 3000, 3061, "DEL", 100, "G" * 100, True),
    ]
    assert stats["dropped"]["mate_duplicate"] == 1 and stats["dropped"]["small"] == 1


def test_sequence_resolved_alleles(tmp_path):
    vcf = _vcf(tmp_path, [
        ("1", 100, "a", "A", "A" + SEQ60, ".", "PASS", "SVTYPE=INS", "GT", "0/0", "0/1"),
        ("1", 200, "b", "A", "A" + SEQ60[:49], ".", "PASS", "SVTYPE=INS", "GT", "0/0", "0/1"),
        ("1", 300, "c", "AT", "A" + SEQ60 + "T", ".", "PASS", ".", "GT", "0/0", "0/1"),
        ("1", 400, "d", "A" + SEQ60, "A", ".", "PASS", "SVTYPE=DEL", "GT", "0/0", "0/1"),
    ])
    rows, stats = _norm(vcf, "generic")
    assert _ins(rows) == [("1", 100, 101, "INS", 60, SEQ60, True), ("1", 300, 301, "INS", 60, SEQ60, True)]
    assert stats["dropped"]["small"] == 1 and stats["dropped"]["unsupported_allele"] == 1


@pytest.mark.parametrize("alt,info,kept", [
    ("<DEL>", "SVTYPE=DEL;END=149", False),   # 49 deleted bases: indel workflow (|len change| < 50)
    ("<DEL>", "SVTYPE=DEL;END=150", True),    # 50 deleted bases
    ("<DUP>", "SVTYPE=DUP;END=149", False),   # 49 duplicated bases
    ("<DUP>", "SVTYPE=DUP;END=150", True),    # 50 duplicated bases (previously dropped)
    ("N[1:150[", "SVTYPE=BND", False),        # 49 deleted bases between the retained ends
    ("N[1:151[", "SVTYPE=BND", True),         # 50
])
def test_size_boundary_is_the_net_length_change(tmp_path, alt, info, kept):
    vcf = _vcf(tmp_path, [("1", 100, "x", "N", alt, ".", "PASS", info, "GT", "0/0", "0/1")])
    rows, stats = _norm(vcf, "generic", min_size=50)
    assert (len(rows) == 1) is kept and stats["dropped"]["small"] == (0 if kept else 1)


def _fasta(tmp_path, seq, name="ref.fa"):
    import pysam
    p = tmp_path / name
    p.write_text(">1\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    pysam.faidx(str(p))
    return str(p)


def test_insertions_in_a_repeat_are_left_aligned_with_a_fasta(tmp_path):
    left = "ACGTTGCAGTCCATGGACTTAGCGATCGTACGGATTCAGCTAGGCATTCAGGACTTCGAT"   # 60 bp, no CA run
    seq = left * 3 + "CA" * 20 + left * 3                                 # (CA)x20 at 181..220
    fa = _fasta(tmp_path, seq)
    ins = "CA" * 30
    vcf = _vcf(tmp_path, [
        # the same +(CA)x30 written at the start, the middle and the end of the repeat
        ("1", 180, "a", seq[179], seq[179] + ins, ".", "PASS", ".", "GT", "0/0", "0/1"),
        ("1", 200, "b", seq[199], seq[199] + ins, ".", "PASS", ".", "GT", "0/0", "0/1"),
        ("1", 220, "c", seq[219], f"{seq[219]}{ins}[1:221[", ".", "PASS", "SVTYPE=BND", "GT", "0/0", "0/1"),
    ])
    rows, stats = sv_normalize.normalize(vcf, "generic", "T1", "N1", order={"1": 0}, fasta=fa)
    assert _ins(rows) == [("1", 180, 181, "INS", 60, ins, True)]
    assert stats["dropped"]["mate_duplicate"] == 2 and stats["insertions"]["left_aligned"] == 2
    rows, _ = sv_normalize.normalize(vcf, "generic", "T1", "N1", order={"1": 0})       # no FASTA
    assert [r["pos1"] for r in rows] == [180, 200, 220]


def test_cli_columns_include_insertion_fields(tmp_path):
    vcf = _vcf(tmp_path, [("1", 900, "a", "N", "<INS>", ".", "PASS", "SVTYPE=INS;SVLEN=300", "GT", "0/0", "0/1")])
    out, st = tmp_path / "o.tsv", tmp_path / "o.json"
    sv_normalize.main(["--vcf", vcf, "--caller", "generic", "--tumor-id", "T1", "--normal-id", "N1",
                       "--out-tsv", str(out), "--out-stats", str(st)])
    header, row = out.read_text().splitlines()
    assert header.split("\t")[-2:] == ["ins_seq", "ins_seq_available"]
    assert row.split("\t")[-3:] == ["300", "", "False"]
