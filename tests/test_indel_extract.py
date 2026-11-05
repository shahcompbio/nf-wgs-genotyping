import json

import pytest

import indel_extract

HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _vcf(tmp_path, lines, name="in.vcf"):
    p = tmp_path / name
    p.write_text(HEAD + "".join("\t".join(map(str, l)) + "\n" for l in lines))
    return str(p)


def _keys(rows):
    return [(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in rows]


def test_pindel_keeps_pass_indels_normalized(synth):
    rows, stats = indel_extract.extract(synth["indel_vcfs"]["pindel"], "pindel", "T1", "N1",
                                        fasta=synth["fasta"], contigs=["1", "2", "3"])
    s = synth["indel_sites"]
    assert _keys(rows) == [tuple(s[k]) for k in ("D1", "HP1", "I3", "C4")]
    assert all(r["tumor_id"] == "T1" and r["normal_id"] == "N1" and r["caller"] == "pindel" for r in rows)
    assert stats["records"] == 6 and stats["kept"] == 4
    assert (stats["dropped"]["not_pass"], stats["dropped"]["not_indel"]) == (1, 1)
    assert stats["normalized"] == 1                           # C4 lost its shared first base


def test_strelka_right_shifted_insertion_is_left_aligned(synth):
    rows, stats = indel_extract.extract(synth["indel_vcfs"]["strelka"], "strelka", "T1", "N1",
                                        fasta=synth["fasta"], contigs=["1", "2", "3"])
    s = synth["indel_sites"]
    assert _keys(rows) == [tuple(s["D1"]), tuple(s["I3"])]
    assert stats["dropped"]["off_contig"] == 1 and stats["normalized"] == 1


def test_mutect_multiallelic_split_mnv_and_long_events(synth):
    rows, stats = indel_extract.extract(synth["indel_vcfs"]["mutect"], "mutect", "T2", "N1",
                                        fasta=synth["fasta"], max_length=50)
    s = synth["indel_sites"]
    assert _keys(rows) == [tuple(s["HP1"]), tuple(s["HP2"])]
    assert (stats["dropped"]["mnv"], stats["dropped"]["too_long"]) == (1, 1)
    rows, _ = indel_extract.extract(synth["indel_vcfs"]["mutect"], "mutect", "T2", "N1",
                                    fasta=synth["fasta"], max_length=100)
    assert len(rows) == 3                                     # the 60-bp deletion is kept


def test_equivalent_representations_collapse_to_one_row(synth, tmp_path):
    seq = synth["seqs"]["3"]
    vcf = _vcf(tmp_path, [("3", 1000, ".", "CA", "C", ".", "PASS", "."),
                          ("3", 1007, ".", "AA", "A", ".", "PASS", "."),        # same deletion, run end
                          ("3", 1000, ".", "CAAAAAAAAG", "CAAAAAAAG", ".", "PASS", ".")])
    assert seq[999:1009] == "CAAAAAAAAG"
    rows, stats = indel_extract.extract(vcf, "generic", "T", "N", fasta=synth["fasta"])
    assert _keys(rows) == [("3", 1000, "CA", "C")]
    assert stats["dropped"]["duplicate"] == 2


def test_ref_mismatch_is_dropped_and_fails_above_threshold(synth, tmp_path):
    seq = synth["seqs"]["3"]
    good = [("3", p, ".", seq[p - 1:p + 1], seq[p - 1], ".", "PASS", ".") for p in range(100, 300, 10)]
    wrong = [b for b in "ACGT" if b != seq[399]][0]           # REF base that is not the FASTA's
    bad = [("3", 400, ".", wrong + seq[400], wrong, ".", "PASS", ".")]
    vcf = _vcf(tmp_path, good + bad)
    with pytest.raises(ValueError, match="REF mismatch"):
        indel_extract.extract(vcf, "pindel", "T", "N", fasta=synth["fasta"], max_ref_mismatch_frac=0.01)
    rows, stats = indel_extract.extract(vcf, "pindel", "T", "N", fasta=synth["fasta"],
                                        max_ref_mismatch_frac=0.1)
    assert len(rows) == 20 and stats["dropped"]["ref_mismatch"] == 1


def test_missing_filter_counts_as_pass_only_when_asked(synth, tmp_path):
    seq = synth["seqs"]["3"]
    vcf = _vcf(tmp_path, [("3", 100, ".", seq[99:101], seq[99], ".", ".", ".")])
    assert indel_extract.extract(vcf, "generic", "T", "N", fasta=synth["fasta"])[1]["kept"] == 0
    assert indel_extract.extract(vcf, "generic", "T", "N", fasta=synth["fasta"],
                                 missing_filter_is_pass=True)[1]["kept"] == 1


def test_cli_requires_reference_and_writes_tsv_and_stats(synth, tmp_path):
    out, st = tmp_path / "o.tsv", tmp_path / "o.json"
    args = ["--vcf", synth["indel_vcfs"]["pindel"], "--caller", "pindel", "--tumor-id", "T1",
            "--normal-id", "N1", "--contigs", "1,2,3", "--out-tsv", str(out), "--out-stats", str(st)]
    with pytest.raises(SystemExit):
        indel_extract.main(args)                              # --reference is required
    indel_extract.main(args + ["--reference", synth["fasta"], "--max-length", "50"])
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == ["chrom", "pos", "ref", "alt", "tumor_id", "normal_id", "caller"]
    assert len(lines) == 5
    stats = json.loads(st.read_text())
    assert stats["kept"] == 4 and stats["variant_type"] == "indel"
