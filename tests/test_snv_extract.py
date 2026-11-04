import gzip
import json

import pytest

import snv_extract

HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _vcf(tmp_path, lines, name="in.vcf", gz=False):
    text = HEAD + "".join("\t".join(map(str, l)) + "\n" for l in lines)
    p = tmp_path / name
    if gz:
        with gzip.open(p, "wt") as fh:
            fh.write(text)
    else:
        p.write_text(text)
    return str(p)


def _ref(synth, chrom, pos):
    return synth["seqs"][chrom][pos - 1]


def _alt(base, k=0):
    return [b for b in "ACGT" if b != base][k]


def test_keeps_only_pass_biallelic_snvs(synth, tmp_path):
    r = lambda p: _ref(synth, "1", p)  # noqa: E731
    vcf = _vcf(tmp_path, [
        ("1", 100, ".", r(100), _alt(r(100)), ".", "PASS", "."),
        ("1", 110, ".", r(110), _alt(r(110)), ".", "LowQual", "."),
        ("1", 120, ".", r(120), r(120) + "T", ".", "PASS", "."),
        ("1", 130, ".", r(130) + r(131), "AA", ".", "PASS", "."),
        ("1", 140, ".", r(140), f"{_alt(r(140), 0)},{_alt(r(140), 1)}", ".", "PASS", "."),
        ("1", 150, ".", r(150), "<DEL>", ".", "PASS", "."),
        ("1", 160, ".", r(160), _alt(r(160)), ".", ".", "."),
        ("1", 170, ".", r(170), _alt(r(170)), ".", "PASS;LowQual", "."),
    ])
    rows, stats = snv_extract.extract(vcf, "mutect", "T1", "N1", fasta=synth["fasta"])
    assert [(x["chrom"], x["pos"], x["alt"]) for x in rows] == [
        ("1", 100, _alt(r(100))), ("1", 140, _alt(r(140), 0)), ("1", 140, _alt(r(140), 1))]
    assert all(x["tumor_id"] == "T1" and x["normal_id"] == "N1" and x["caller"] == "mutect" for x in rows)
    assert stats["records"] == 8
    assert stats["kept"] == 3
    assert stats["dropped"]["not_pass"] == 3
    assert stats["dropped"]["not_snv"] == 3


def test_missing_filter_counts_as_pass_only_when_asked(synth, tmp_path):
    p = 160
    vcf = _vcf(tmp_path, [("1", p, ".", _ref(synth, "1", p), _alt(_ref(synth, "1", p)), ".", ".", ".")])
    assert snv_extract.extract(vcf, "generic", "T", "N")[1]["kept"] == 0
    assert snv_extract.extract(vcf, "generic", "T", "N", missing_filter_is_pass=True)[1]["kept"] == 1


def test_off_contig_records_are_dropped(synth, tmp_path):
    vcf = _vcf(tmp_path, [("MT", 10, ".", "A", "G", ".", "PASS", "."),
                          ("1", 5, ".", _ref(synth, "1", 5), _alt(_ref(synth, "1", 5)), ".", "PASS", ".")])
    rows, stats = snv_extract.extract(vcf, "strelka", "T", "N", contigs=["1", "2"])
    assert [x["chrom"] for x in rows] == ["1"]
    assert stats["dropped"]["off_contig"] == 1


def test_ref_mismatch_is_dropped_and_fails_above_threshold(synth, tmp_path):
    good = [("1", p, ".", _ref(synth, "1", p), _alt(_ref(synth, "1", p)), ".", "PASS", ".")
            for p in range(200, 220)]
    wrong = _alt(_ref(synth, "1", 300))
    bad = [("1", 300, ".", wrong, _alt(wrong), ".", "PASS", ".")]
    vcf = _vcf(tmp_path, good + bad)
    with pytest.raises(ValueError, match="REF mismatch"):
        snv_extract.extract(vcf, "caveman", "T", "N", fasta=synth["fasta"], max_ref_mismatch_frac=0.01)
    rows, stats = snv_extract.extract(vcf, "caveman", "T", "N", fasta=synth["fasta"],
                                      max_ref_mismatch_frac=0.1)
    assert len(rows) == 20
    assert stats["dropped"]["ref_mismatch"] == 1


def test_duplicate_records_collapse_and_gz_input(synth, tmp_path):
    p = 400
    rec = ("1", p, ".", _ref(synth, "1", p), _alt(_ref(synth, "1", p)), ".", "PASS", ".")
    vcf = _vcf(tmp_path, [rec, rec], name="in.vcf.gz", gz=True)
    rows, stats = snv_extract.extract(vcf, "mutect", "T", "N")
    assert len(rows) == 1
    assert stats["dropped"]["duplicate"] == 1


def test_cli_writes_tsv_and_stats(synth, tmp_path):
    out, st = tmp_path / "o.tsv", tmp_path / "o.json"
    snv_extract.main(["--vcf", synth["snv_vcfs"]["caveman"], "--caller", "caveman", "--tumor-id", "T1",
                      "--normal-id", "N1", "--reference", synth["fasta"], "--contigs", "1,2",
                      "--out-tsv", str(out), "--out-stats", str(st)])
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == ["chrom", "pos", "ref", "alt", "tumor_id", "normal_id", "caller"]
    assert len(lines) == 3
    assert json.loads(st.read_text())["kept"] == 2
