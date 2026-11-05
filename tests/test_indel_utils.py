"""Indel normalization (left-align + trim), repeat-aware variable spans and haplotypes, on
hand-written sequences (positions in comments are 1-based VCF positions)."""
import pytest

import indel_utils as iu

#            1234567890123456789012345678901234567890
SEQ = ("GGATC"            # 1-5
       "CAAAAAAAAG"       # 6-15: pad C (6), A x8 (7-14), G (15)
       "TCTGA"            # 16-20
       "TCACACACACACAG"   # 21-34: pad T (21), (CA) x6 (22-33), G (34)
       "GTCAGTT"          # 35-41: G (35), TCA (36-38), G (39) -> unique 3-bp deletion
       "GTACGTTCCT")      # 42-51: G (42) TACG (43-46) T (47)


def fetch(start0, end0):
    return SEQ[max(start0, 0):end0]


@pytest.mark.parametrize("given,expected", [
    ((13, "AA", "A"), (6, "CA", "C")),            # homopolymer deletion written at the run's end
    ((10, "AAA", "A"), (6, "CAA", "C")),          # 2-bp deletion written mid-run
    ((33, "A", "ACA"), (21, "T", "TCA")),         # dinucleotide insertion written right-shifted
    ((6, "CA", "C"), (6, "CA", "C")),             # already normalized
    ((35, "GTCAG", "GG"), (35, "GTCA", "G")),     # shared trailing base trimmed
    ((42, "GTACG", "GC"), (43, "TACG", "C")),     # complex: shared leading base trimmed
    ((7, "AAAAAAAAG", "AAAAAAAG"), (6, "CA", "C")),  # over-long alleles: trimmed and left-aligned
])
def test_normalize_left_aligns_and_trims(given, expected):
    assert iu.normalize(fetch, *given) == expected


def test_normalize_stops_at_contig_start():
    # deleting the first G of "GG..." cannot be shifted left of position 1
    assert iu.normalize(fetch, 1, "GG", "G") == (1, "GG", "G")
    assert iu.normalize(fetch, 2, "GA", "A") == (1, "GG", "G")


@pytest.mark.parametrize("allele,span", [
    ((6, "CA", "C"), (6, 14)),          # 0-based half-open: the whole A run (7-14)
    ((6, "CAA", "C"), (6, 14)),
    ((6, "C", "CA"), (6, 14)),          # insertion into the run: anywhere in it
    ((21, "T", "TCA"), (21, 33)),       # the whole (CA) x6 repeat (22-33)
    ((35, "GTCA", "G"), (35, 38)),      # unique deletion: exactly the deleted bases (36-38)
    ((43, "TACG", "C"), (42, 46)),      # complex: the changed bases (43-46)
    ((19, "G", "GC"), (19, 19)),        # unique insertion between 19 (G) and 20 (A): empty span
])
def test_variable_span_covers_every_equivalent_placement(allele, span):
    assert iu.variable_span(fetch, len(SEQ), *allele) == span


def test_variable_span_grows_its_window_for_long_repeats():
    seq = "C" + "A" * 300 + "G" + "TTCCGGAT"
    f = lambda a, b: seq[max(a, 0):b]           # noqa: E731
    assert iu.variable_span(f, len(seq), 1, "CA", "C", window=8) == (1, 301)


def test_footprint_adds_one_flanking_base_each_side():
    assert iu.footprint((6, 14), 51) == (5, 15)
    assert iu.footprint((19, 19), 51) == (18, 20)
    assert iu.footprint((0, 3), 51) == (0, 4)


def test_apply_allele_builds_the_alternate_haplotype_window():
    ws, we = 3, 20
    assert iu.apply_allele(SEQ[ws:we], ws, 6, "CA", "C") == SEQ[3:6] + SEQ[7:20]
    assert iu.apply_allele(SEQ[ws:we], ws, 6, "C", "CAA") == SEQ[3:6] + "AA" + SEQ[6:20]
    with pytest.raises(ValueError, match="outside"):
        iu.apply_allele(SEQ[10:20], 10, 6, "CA", "C")


@pytest.mark.parametrize("ref,alt,kind", [
    ("A", "G", "snv"), ("CA", "C", "indel"), ("C", "CAT", "indel"), ("GTACG", "GC", "indel"),
    ("AC", "GT", "mnv"), ("A", "<DEL>", "other"), ("A", "*", "other"), ("N[1:5[", "A", "other"),
    ("A", "A", "other"),
])
def test_allele_kind(ref, alt, kind):
    assert iu.allele_kind(ref, alt) == kind
