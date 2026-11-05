"""Indel representation helpers shared by indel_extract, indel_union and indel_count.

Coordinates: VCF positions are 1-based; spans/footprints are 0-based half-open reference
intervals. ``fetch(start0, end0)`` returns the reference bases of one contig (upper case);
out-of-range ends are clipped by the caller (pysam.FastaFile.fetch and str slicing both do).
"""
import re

_PLAIN = re.compile(r"^[ACGT]+$")


def allele_kind(ref, alt):
    """'snv', 'indel' (lengths differ, plain bases; complex del-ins included), 'mnv' (equal
    lengths > 1) or 'other' (symbolic, breakend, '*', non-ACGT bases, REF == ALT)."""
    if not (_PLAIN.match(ref) and _PLAIN.match(alt)) or ref == alt:
        return "other"
    if len(ref) == len(alt):
        return "snv" if len(ref) == 1 else "mnv"
    return "indel"


def normalize(fetch, pos, ref, alt):
    """Left-align and trim a biallelic indel (Tan, Abecasis & Kang, Bioinformatics 2015; the rule
    of `vt normalize` / `bcftools norm`). Returns (pos, ref, alt) with one shared leading base
    when the event is a pure insertion/deletion. At position 1 the allele cannot be extended to
    the left, so shifting stops there."""
    while True:
        changed = False
        if ref and alt and ref[-1] == alt[-1]:
            if (len(ref) == 1 or len(alt) == 1) and pos == 1:
                break
            ref, alt = ref[:-1], alt[:-1]
            changed = True
        if not ref or not alt:
            base = fetch(pos - 2, pos - 1).upper()
            ref, alt, pos = base + ref, base + alt, pos - 1
            changed = True
        if not changed:
            break
    while len(ref) >= 2 and len(alt) >= 2 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def _common_prefix(a, b):
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


def _common_suffix(a, b):
    n = min(len(a), len(b))
    i = 0
    while i < n and a[-1 - i] == b[-1 - i]:
        i += 1
    return i


def variable_span(fetch, contig_len, pos, ref, alt, window=50):
    """0-based half-open reference interval covering every placement of the event that yields the
    same haplotype (e.g. the whole homopolymer for a 1-bp deletion in it; empty for an insertion
    with a unique placement). The window around the event doubles until the repeat no longer
    reaches its edge (or it covers the contig)."""
    p0 = pos - 1
    while True:
        ws, we = max(0, p0 - window), min(contig_len, p0 + len(ref) + window)
        r = fetch(ws, we).upper()
        a = r[:p0 - ws] + alt + r[p0 - ws + len(ref):]
        cp, cs = _common_prefix(r, a), _common_suffix(r, a)
        n = min(len(r), len(a))
        # a repeat that runs into the window edge makes the prefix (or suffix) cover all of the
        # shorter sequence; grow the window unless that edge is the contig's end
        if not ((cp >= n and we < contig_len) or (cs >= n and ws > 0)):
            break
        window *= 2
    d = max(0, len(r) - len(a))
    start = min(cp, len(r) - cs - d)
    end = max(cp + d, len(r) - cs)
    return ws + start, ws + end


def footprint(span, contig_len):
    """Span plus one flanking reference base on each side: the bases a read must reach to be
    examined at the locus (for an insertion, the two bases around it)."""
    return max(0, span[0] - 1), min(contig_len, span[1] + 1)


def apply_allele(window_seq, ws, pos, ref, alt):
    """Haplotype of reference window [ws, ws + len(window_seq)) carrying the allele."""
    p0 = pos - 1
    if p0 < ws or p0 + len(ref) > ws + len(window_seq):
        raise ValueError(f"allele {pos}:{ref}>{alt} is outside the window starting at {ws}")
    return window_seq[:p0 - ws] + alt + window_seq[p0 - ws + len(ref):]
