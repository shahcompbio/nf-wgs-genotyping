#!/usr/bin/env python3
"""Deterministic synthetic fixtures for nf-wgs-genotyping tests.

    python tests/make_fixtures.py tests/data/synthetic

SNV fixtures are fully synthetic (seeded random reference, hand-built reads with known
composition). SV fixtures are caller-style VCFs (BRASS/GRIDSS/SvABA conventions, VCF-spec BND
positions) describing deletions from svtyper's MIT-licensed NA12878 fixture
(tests/data/svtyper), so they can be genotyped against NA12878.target_loci.sorted.bam.
Nothing here is derived from patient data.
"""
import argparse
import array
import csv
import json
import os
import random

import pysam

SEED = 20261104
# contig 3 is generated after 1 and 2 by the same RNG (so 1 and 2 are unchanged), then the indel
# motifs below are written into it (see _design_indel_contig)
CONTIGS = [("1", 3000), ("2", 3000), ("3", 3000)]
READ_LEN = 50
TARGET = "T1"          # target sample id == BAM @RG SM
NORMAL = "N1"
OTHER_TUMOR = "T2"

S1 = ("1", 1000)       # complex composition site
S2 = ("1", 1500)       # 4 reference fragments only
S3 = ("2", 500)        # no coverage; multiallelic call

# Expected counts at S1 (see _snv_reads): fragments counted once, explicit filters.
TRUTH_S1 = {"ref_count": 10, "alt_count": 9, "other_count": 1, "conflict_count": 1}
TRUTH_S2 = {"ref_count": 4, "alt_count": 0, "other_count": 0, "conflict_count": 0}

# Indel loci on contig 3 (1-based VCF coordinates, normalized = left-aligned, minimal):
#   D1   unique 3-bp deletion            pad G at 500, TCA deleted, G at 504
#   HP   A x8 at 1001-1008 (pad C 1000, G 1009): two alleles, -A (HP1) and -AA (HP2)
#   I3   (CA) x6 at 1501-1512 (pad T 1500, G 1513): +CA
#   C4   complex del-ins GTACG -> GC at 2000 (normalized: 2001 TACG -> C), T at 2005
#   LONG 60-bp deletion at 2500 (>= indel_max_length 50: dropped by extraction)
INDEL_SITES = {
    "D1": ("3", 500, "GTCA", "G"),
    "HP1": ("3", 1000, "CA", "C"),
    "HP2": ("3", 1000, "CAA", "C"),
    "I3": ("3", 1500, "T", "TCA"),
    "C4": ("3", 2001, "TACG", "C"),
}
# Expected fragment counts (see _indel_reads). other = fragments matching another allele of the
# locus; uninformative/conflict are per locus (HP1 and HP2 share them).
TRUTH_INDELS = {
    "D1": {"ref_count": 8, "alt_count": 7, "other_count": 0, "uninformative_count": 1, "conflict_count": 1},
    "HP1": {"ref_count": 6, "alt_count": 4, "other_count": 3, "uninformative_count": 2, "conflict_count": 0},
    "HP2": {"ref_count": 6, "alt_count": 3, "other_count": 4, "uninformative_count": 2, "conflict_count": 0},
    "I3": {"ref_count": 5, "alt_count": 5, "other_count": 0, "uninformative_count": 1, "conflict_count": 0},
    "C4": {"ref_count": 3, "alt_count": 3, "other_count": 0, "uninformative_count": 0, "conflict_count": 0},
}

P, PROPER, UNMAP, MUNMAP, REV, MREV, R1, R2 = 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80
SECONDARY, QCFAIL, DUP, SUPP = 0x100, 0x200, 0x400, 0x800

# NA12878 deletions from svtyper tests/data/example.gt.vcf: (id, chrom, POS, END, NA12878 GT)
NA12878_DELS = {
    "72354": ("2", 2798774, 2799241, "1/1"),
    "72623": ("2", 3732805, 3738027, "0/1"),
    "73459": ("2", 6400508, 6402192, "0/0"),
    "76210": ("2", 14704171, 14710115, "0/1"),
    "76667": ("2", 18906149, 18907347, "0/1"),
    "81625": ("2", 30934974, 30935094, "0/0"),
    "82677": ("2", 33224462, 33227291, "0/1"),
}
SV_TARGET = "NA12878"   # SM of the svtyper fixture BAM


# ----------------------------------------------------------------------------- reference
def _reference(outdir):
    rng = random.Random(SEED)
    seqs = {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in CONTIGS}
    _design_indel_contig(seqs)
    path = os.path.join(outdir, "ref.fa")
    with open(path, "w") as fh:
        for c, s in seqs.items():
            fh.write(f">{c}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    pysam.faidx(path)
    return path, seqs


def _bases(ref_base):
    others = [b for b in "ACGT" if b != ref_base]
    return others[0], others[1]          # alt, other


def _put(seq, pos1, motif):
    i = pos1 - 1
    return seq[:i] + motif + seq[i + len(motif):]


def _design_indel_contig(seqs):
    """Write the indel motifs into contig 3 (see INDEL_SITES)."""
    s = seqs["3"]
    s = _put(s, 500, "GTCAG")                    # D1: pad G, TCA deleted, G (G != A, T != G: unique)
    s = _put(s, 1000, "C" + "A" * 8 + "G")       # HP: pad C, A x8, G
    s = _put(s, 1500, "T" + "CA" * 6 + "G")      # I3: pad T, (CA) x6, G
    s = _put(s, 2000, "GTACGT")                  # C4: GTACG -> GC, then T
    seqs["3"] = s


def _hap_read(seq, start0, var=None, n=READ_LEN, aligned=None):
    """(query, cigar) of an n-base read whose first base is at reference 0-based start0, taken
    from the haplotype that carries var = (pos1, REF, ALT) (VCF style, shared first base), or
    from the reference. CIGAR: M through the shared base, I/D for the rest of ALT/REF, then M.
    aligned=k instead writes M k + S (n - k): an aligner that soft-clipped the indel."""
    if var is None:
        return seq[start0:start0 + n], [(0, n)]
    pos1, ref, alt = var
    p0 = pos1 - 1
    assert ref[0] == alt[0] and start0 <= p0, (var, start0)
    hap = seq[:p0] + alt + seq[p0 + len(ref):]
    query = hap[start0:start0 + n]
    if aligned is not None:
        return query, [(0, aligned), (4, n - aligned)]
    m1 = p0 + 1 - start0
    ins, dele = len(alt) - 1, len(ref) - 1
    rest = n - m1 - ins
    assert rest > 0, (var, start0)
    cigar = [(0, m1)] + ([(1, ins)] if ins else []) + ([(2, dele)] if dele else []) + [(0, rest)]
    return query, cigar


def _indel_reads(header, seqs):
    """Reads at the contig-3 indel loci; the comments give each fragment's expected class."""
    chrom, s = "3", seqs["3"]
    reads, n = [], [0]

    def frag(start0, var=None, *, aligned=None, flag_extra=0, mapq=60, mask_from=None,
             mate_start0=None, mate_var=None):
        """read1 forward at start0 from haplotype `var`; read2 reverse at mate_start0 (default
        start0 + 250, outside every locus) from haplotype `mate_var` (default: reference)."""
        n[0] += 1
        name = f"indel{n[0]:03d}"
        m0 = mate_start0 if mate_start0 is not None else start0 + 250
        f1 = P | R1 | MREV | PROPER | flag_extra
        f2 = P | R2 | REV | PROPER | flag_extra
        q1, c1 = _hap_read(s, start0, var, aligned=aligned)
        quals = [35] * READ_LEN
        if mask_from is not None:                      # BQ 2 from this query index on
            quals[mask_from:] = [2] * (READ_LEN - mask_from)
        q2, c2 = _hap_read(s, m0, mate_var)
        reads.append(_seg(header, name, chrom, start0, q1, f1, mapq, quals, (chrom, m0), cigar=c1))
        reads.append(_seg(header, name, chrom, m0, q2, f2, mapq, None, (chrom, start0), cigar=c2))

    d1 = INDEL_SITES["D1"][1:]
    for st in (462, 464, 468, 472, 474, 476, 478):
        frag(st)                                       # D1 ref: 7
    frag(466, mate_start0=476)                         # D1 ref, overlapping mates: +1 (once)
    for st in (461, 465, 469, 473):
        frag(st, d1)                                   # D1 alt (CIGAR D): 4
    frag(463, d1, mate_start0=471, mate_var=d1)        # D1 alt, overlapping mates: +1 (once)
    frag(455, d1, aligned=45)                          # D1 alt, clipped right after the pad base: +1
    frag(455, d1, aligned=40)                          # D1 alt, clipped 5 bases before the pad: +1
    frag(467, d1, mask_from=33)                        # D1 alt with BQ 2 after the pad: uninformative
    frag(470, d1, mate_start0=470)                     # mates alt vs ref: conflict
    for st in (464, 466):
        frag(st, d1, flag_extra=DUP)                   # duplicates: excluded
    frag(468, d1, mapq=10)                             # MAPQ 10: excluded
    frag(472, d1, flag_extra=QCFAIL)                   # QC fail: excluded
    for kind, flag in (("isec", SECONDARY), ("isup", SUPP)):   # primaries elsewhere: excluded
        name = f"{kind}001"
        reads.append(_seg(header, name, chrom, 2700, s[2700:2750], P | R1 | MREV | PROPER, mate=(chrom, 2850)))
        reads.append(_seg(header, name, chrom, 2850, s[2850:2900], P | R2 | REV | PROPER, mate=(chrom, 2700)))
        q, c = _hap_read(s, 471, d1)
        reads.append(_seg(header, name, chrom, 471, q, P | R1 | MREV | PROPER | flag, mate=(chrom, 2850), cigar=c))

    hp1, hp2 = INDEL_SITES["HP1"][1:], INDEL_SITES["HP2"][1:]
    hp1_right = (1007, "AA", "A")                      # same haplotype, deletion at the run's end
    for st in (966, 968, 970, 972, 974, 976):
        frag(st)                                       # HP ref: 6
    frag(967, hp1)
    frag(971, hp1)
    frag(969, hp1_right)
    frag(973, hp1_right)                               # HP1 alt: 4 (HP2 other: 4)
    for st in (965, 975, 977):
        frag(st, hp2)                                  # HP2 alt: 3 (HP1 other: 3)
    frag(955)                                          # ends inside the A run: uninformative
    frag(1003)                                         # starts inside the A run: uninformative

    i3 = INDEL_SITES["I3"][1:]
    i3_right = (1512, "A", "ACA")                      # same haplotype, insertion at the repeat's end
    for st in (1468, 1470, 1472, 1474, 1476):
        frag(st)                                       # I3 ref: 5
    frag(1467, i3)
    frag(1471, i3)
    frag(1469, i3_right)
    frag(1473, i3_right)                               # I3 alt (CIGAR I): 4
    frag(1472, i3, aligned=28)                         # I3 alt, insertion soft-clipped: +1
    frag(1458)                                         # ends inside the repeat: uninformative

    c4 = (2000, "GTACG", "GC")                         # caller representation (shared pad base)
    for st in (1965, 1968, 1971):
        frag(st)                                       # C4 ref: 3
    for st in (1966, 1969, 1972):
        frag(st, c4)                                   # C4 alt: 3
    return reads


# ----------------------------------------------------------------------------- reads
def _seg(header, name, chrom, start0, seq, flag, mapq=60, quals=None, mate=None, cigar=None):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.query_sequence = seq
    a.flag = flag
    a.reference_id = header.get_tid(chrom)
    a.reference_start = start0
    a.mapping_quality = mapq
    a.cigartuples = cigar or [(0, len(seq))]
    a.query_qualities = array.array("B", quals or [35] * len(seq))
    if mate is not None:
        mchrom, mstart0 = mate
        a.next_reference_id = header.get_tid(mchrom)
        a.next_reference_start = mstart0
        if mchrom == chrom:
            lo = min(start0, mstart0)
            hi = max(start0, mstart0) + READ_LEN
            a.template_length = (hi - lo) if start0 <= mstart0 else -(hi - lo)
    a.set_tag("RG", TARGET)
    return a


def _read_seq(seqs, chrom, start0, site=None, base=None):
    s = seqs[chrom][start0:start0 + READ_LEN]
    if site is not None:
        off = site - 1 - start0
        s = s[:off] + base + s[off + 1:]
    return s


def _quals(start0, site, bq):
    q = [35] * READ_LEN
    q[site - 1 - start0] = bq
    return q


def _snv_reads(header, seqs):
    chrom, pos = S1
    ref = seqs[chrom][pos - 1]
    alt, other = _bases(ref)
    reads = []
    n = [0]

    def pair(base, start0, *, flag_extra=0, mapq=60, bq=35, mate_start0=None, mate_chrom=None,
             proper=True, mate_base=None, mate_bq=35):
        """read1 forward covering the site; read2 reverse at mate_start0 (covers the site only
        if mate_base is given)."""
        n[0] += 1
        name = f"frag{n[0]:03d}"
        mchrom = mate_chrom or chrom
        m0 = mate_start0 if mate_start0 is not None else start0 + 250
        f1 = P | R1 | MREV | (PROPER if proper else 0) | flag_extra
        f2 = P | R2 | REV | (PROPER if proper else 0) | flag_extra
        r1 = _seg(header, name, chrom, start0, _read_seq(seqs, chrom, start0, pos, base), f1, mapq,
                  _quals(start0, pos, bq), (mchrom, m0))
        if mate_base is not None:
            s2, q2 = _read_seq(seqs, chrom, m0, pos, mate_base), _quals(m0, pos, mate_bq)
        else:
            s2, q2 = _read_seq(seqs, mchrom, m0), None
        r2 = _seg(header, name, mchrom, m0, s2, f2, mapq, q2, (chrom, start0))
        reads.extend([r1, r2])

    for i in range(10):
        pair(ref, 975 + i)                                  # 10 ref
    for i in range(5):
        pair(alt, 960 + i)                                  # 5 alt
    for i in range(2):
        pair(alt, 980 + i, mate_start0=990 + i, mate_base=alt)   # overlapping mates agree: +2 alt
    pair(alt, 978, mate_start0=992, mate_base=ref, bq=30, mate_bq=30)   # tie -> conflict
    pair(alt, 977, mate_start0=993, mate_base=ref, bq=35, mate_bq=25)   # higher BQ alt -> +1 alt
    pair(other, 970)                                        # other: 1
    for i in range(3):
        pair(alt, 965 + i, flag_extra=DUP)                  # duplicates: excluded
    for i in range(2):
        pair(alt, 968 + i, mapq=10)                         # MAPQ 10: excluded
    for i in range(2):
        pair(alt, 971 + i, bq=10)                           # BQ 10: excluded
    pair(alt, 974, flag_extra=QCFAIL)                       # qcfail: excluded
    pair(alt, 976, proper=False, mate_chrom="2", mate_start0=100)   # improper pair: +1 alt

    # secondary and supplementary records at the site whose primaries lie elsewhere: excluded
    for kind, flag in (("sec", SECONDARY), ("sup", SUPP)):
        name = f"{kind}001"
        reads.append(_seg(header, name, chrom, 2400, _read_seq(seqs, chrom, 2400), P | R1 | MREV | PROPER,
                          mate=(chrom, 2600)))
        reads.append(_seg(header, name, chrom, 2600, _read_seq(seqs, chrom, 2600), P | R2 | REV | PROPER,
                          mate=(chrom, 2400)))
        reads.append(_seg(header, name, chrom, 979, _read_seq(seqs, chrom, 979, pos, alt),
                          P | R1 | MREV | PROPER | flag, mate=(chrom, 2600)))

    # S2: 4 reference fragments
    c2, p2 = S2
    for i in range(4):
        pair_start = 1480 + i
        n[0] += 1
        name = f"frag{n[0]:03d}"
        reads.append(_seg(header, name, c2, pair_start, _read_seq(seqs, c2, pair_start),
                          P | R1 | MREV | PROPER, mate=(c2, pair_start + 250)))
        reads.append(_seg(header, name, c2, pair_start + 250, _read_seq(seqs, c2, pair_start + 250),
                          P | R2 | REV | PROPER, mate=(c2, pair_start)))
    return reads, ref, alt


def _bam(outdir, seqs):
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": c, "LN": n} for c, n in CONTIGS],
        "RG": [{"ID": TARGET, "SM": TARGET}],
    })
    reads, ref, alt = _snv_reads(header, seqs)
    reads += _indel_reads(header, seqs)
    tmp = os.path.join(outdir, "target.unsorted.bam")
    with pysam.AlignmentFile(tmp, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    bam = os.path.join(outdir, "target.bam")
    pysam.sort("--no-PG", "-o", bam, tmp)
    os.remove(tmp)
    pysam.index(bam)
    return bam, ref, alt


# ----------------------------------------------------------------------------- SNV VCFs
_VCF_HEAD = ("##fileformat=VCFv4.2\n"
             "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
             "##contig=<ID=1,length=3000>\n##contig=<ID=2,length=3000>\n##contig=<ID=3,length=3000>\n")


def _vcf(path, records, samples=("NORMAL", "TUMOR")):
    with open(path, "w") as fh:
        fh.write(_VCF_HEAD)
        fh.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                            *samples]) + "\n")
        for chrom, pos, ref, alt, filt in records:
            fh.write("\t".join([chrom, str(pos), ".", ref, alt, ".", filt, ".", "GT", "0/0", "0/1"])
                     + "\n")
    return path


def _snv_vcfs(outdir, seqs, alt1):
    r = lambda c, p: seqs[c][p - 1]                        # noqa: E731
    a2, _ = _bases(r(*S2))
    s3a, s3b = _bases(r(*S3))
    caveman = _vcf(os.path.join(outdir, "T1_vs_N1.caveman.vcf"), [
        (*S1, r(*S1), alt1, "PASS"),
        (*S2, r(*S2), a2, "PASS"),
        ("1", 1200, r("1", 1200), _bases(r("1", 1200))[0], "MQ"),          # not PASS
    ])
    mutect = _vcf(os.path.join(outdir, "T1_vs_N1.mutect.vcf"), [
        (*S1, r(*S1), alt1, "PASS"),
        ("1", 1300, r("1", 1300), r("1", 1300) + "T", "PASS"),              # indel: dropped
        ("1", 1400, r("1", 1400) + r("1", 1401), "AC" if r("1", 1400) != "A" else "GC", "PASS"),
    ])
    strelka = _vcf(os.path.join(outdir, "T2_vs_N1.strelka.vcf"), [
        ("1", 1700, r("1", 1700), _bases(r("1", 1700))[0], "."),            # missing FILTER
        (*S3, r(*S3), f"{s3a},{s3b}", "PASS"),                              # multiallelic
    ])
    with open(os.path.join(outdir, "T2_vs_N1.strelka.vcf"), "a") as fh:
        fh.write("\t".join(["MT", "10", ".", "A", "G", ".", "PASS", ".", "GT", "0/0", "0/1"]) + "\n")
    rows = [(TARGET, NORMAL, "caveman", caveman), (TARGET, NORMAL, "mutect", mutect),
            (OTHER_TUMOR, NORMAL, "strelka", strelka)]
    ss = os.path.join(outdir, "snv_samplesheet.csv")
    with open(ss, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tumor_id", "normal_id", "caller", "vcf"])
        for t, nrm, caller, path in rows:
            w.writerow([t, nrm, caller, os.path.basename(path)])
    union = {
        (S1[0], S1[1], r(*S1), alt1): {"callers": "caveman;mutect", "sources": [TARGET]},
        (S2[0], S2[1], r(*S2), a2): {"callers": "caveman", "sources": [TARGET]},
        (S3[0], S3[1], r(*S3), s3a): {"callers": "strelka", "sources": [OTHER_TUMOR]},
        (S3[0], S3[1], r(*S3), s3b): {"callers": "strelka", "sources": [OTHER_TUMOR]},
    }
    return ss, {"caveman": caveman, "mutect": mutect, "strelka": strelka}, union


# ----------------------------------------------------------------------------- indel VCFs
def _indel_vcfs(outdir, seqs, snv_ss):
    """Pindel / Strelka / Mutect-style indel VCFs on contig 3 and the combined samplesheet
    (SNV rows of snv_ss + indel rows, with a variant_type column)."""
    s = seqs["3"]
    r = lambda p, k=1: s[p - 1:p - 1 + k]               # noqa: E731
    mnv_alt = "".join(_bases(b)[0] for b in r(2400, 2))
    pindel = _vcf(os.path.join(outdir, "T1_vs_N1.pindel.vcf"), [
        ("3", 500, "GTCA", "G", "PASS"),                  # D1
        ("3", 1000, "CA", "C", "PASS"),                   # HP1
        ("3", 1500, "T", "TCA", "PASS"),                  # I3
        ("3", 2000, "GTACG", "GC", "PASS"),               # C4 (normalizes to 2001 TACG>C)
        ("3", 2200, r(2200, 3), r(2200), "FF010"),        # not PASS
        ("3", 2300, r(2300), _bases(r(2300))[0], "PASS"),  # SNV: not an indel
    ])
    strelka = _vcf(os.path.join(outdir, "T1_vs_N1.strelka.indels.vcf"), [
        ("3", 500, "GTCA", "G", "PASS"),                  # D1 (second caller)
        ("3", 1512, "A", "ACA", "PASS"),                  # I3 right-shifted: normalizes to 1500 T>TCA
        ("MT", 10, "AT", "A", "PASS"),                    # off contig
    ])
    mutect = _vcf(os.path.join(outdir, "T2_vs_N1.mutect.indels.vcf"), [
        ("3", 1000, "CAA", "CA,C", "PASS"),               # multiallelic: HP1 (after trimming) + HP2
        ("3", 2400, r(2400, 2), mnv_alt, "PASS"),         # MNV: dropped
        ("3", 2500, r(2500, 61), r(2500), "PASS"),        # 60-bp deletion: too long
    ])
    rows = []
    with open(snv_ss, newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append([row["tumor_id"], row["normal_id"], row["caller"], row["vcf"], "snv"])
    rows += [[TARGET, NORMAL, "pindel", os.path.basename(pindel), "indel"],
             [TARGET, NORMAL, "strelka", os.path.basename(strelka), "indel"],
             [OTHER_TUMOR, NORMAL, "mutect", os.path.basename(mutect), "indel"]]
    ss = os.path.join(outdir, "samplesheet.csv")
    with open(ss, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tumor_id", "normal_id", "caller", "vcf", "variant_type"])
        w.writerows(rows)
    union = {
        "D1": {"callers": "pindel;strelka", "source_samples": TARGET},
        "HP1": {"callers": "mutect;pindel", "source_samples": f"{TARGET};{OTHER_TUMOR}"},
        "HP2": {"callers": "mutect", "source_samples": OTHER_TUMOR},
        "I3": {"callers": "pindel;strelka", "source_samples": TARGET},
        "C4": {"callers": "pindel", "source_samples": TARGET},
    }
    return ss, {"pindel": pindel, "strelka": strelka, "mutect": mutect}, union


# ----------------------------------------------------------------------------- SV VCFs
def _bnd_pair(chrom, pos, end, rid, shift=0):
    """DEL (POS=pad base, END=last deleted base) as a VCF-spec BND pair: (pos,+)--(end+1,-)."""
    a, b = pos + shift, end + 1 + shift
    return [(chrom, a, f"{rid}_1", f"N[{chrom}:{b}[", f"{rid}_2"),
            (chrom, b, f"{rid}_2", f"]{chrom}:{a}]N", f"{rid}_1")]


def _sv_vcfs(outdir):
    d = NA12878_DELS

    def write(path, samples, recs, fmt=None):
        with open(path, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                                "FORMAT", *samples]) + "\n")
            for r in recs:
                fh.write("\t".join(r) + "\n")
        return path

    def brass(recs):
        out = []
        for chrom, pos, rid, alt, mate in recs:
            out.append([chrom, str(pos), rid, "N", alt, "8", ".",
                        f"SVTYPE=BND;MATEID={mate}", "GT", "0/0", "0/1"])
        return out

    def svaba(recs, imprecise=False):
        out = []
        for chrom, pos, rid, alt, mate in recs:
            info = f"SVTYPE=BND;MATEID={mate}" + (";IMPRECISE" if imprecise else "")
            out.append([chrom, str(pos), rid, "N", alt, "55", "PASS", info, "GT", "0/0", "0/1"])
        return out

    def gridss(recs, normal_vf, tumor_vf):
        out = []
        for chrom, pos, rid, alt, mate in recs:
            out.append([chrom, str(pos), rid, "N", alt, "500", "PASS",
                        f"PARID={mate};CIPOS=0,0;CIRPOS=0,0;EVENT={rid[:-2]}", "GT:VF",
                        f".:{normal_vf}", f".:{tumor_vf}"])
        return out

    files = {}
    files[("T", "brass")] = write(os.path.join(outdir, "NA12878_vs_N1.brass.vcf"), ["NORMAL", "TUMOUR"],
                                  brass(_bnd_pair(*d["72354"][:3], "b72354") + _bnd_pair(*d["76667"][:3], "b76667")))
    files[("T", "svaba")] = write(os.path.join(outdir, "NA12878_vs_N1.svaba.vcf"), ["n.bam", "t.bam"],
                                  svaba(_bnd_pair(*d["72354"][:3], "s72354", shift=15)
                                        + _bnd_pair(*d["76210"][:3], "s76210")))
    files[("T", "gridss")] = write(os.path.join(outdir, "NA12878_vs_N1.gridss.vcf"), [NORMAL, SV_TARGET],
                                   gridss(_bnd_pair(*d["72623"][:3], "g72623"), 0, 12)
                                   + gridss(_bnd_pair(*d["82677"][:3], "g82677"), 25, 20))
    files[("T2", "brass")] = write(os.path.join(outdir, "T2_vs_N1.brass.vcf"), ["NORMAL", "TUMOUR"],
                                   brass(_bnd_pair(*d["73459"][:3], "b73459")))
    files[("T2", "svaba")] = write(os.path.join(outdir, "T2_vs_N1.svaba.vcf"), ["n.bam", "t.bam"],
                                   svaba(_bnd_pair(*d["81625"][:3], "s81625"), imprecise=True)
                                   + svaba(_bnd_pair(*d["72623"][:3], "s72623", shift=-30)))
    files[("T2", "gridss")] = write(os.path.join(outdir, "T2_vs_N1.gridss.vcf"), [NORMAL, OTHER_TUMOR], [])
    ss = os.path.join(outdir, "sv_samplesheet.csv")
    with open(ss, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tumor_id", "normal_id", "caller", "vcf"])
        for (t, caller), path in sorted(files.items()):
            w.writerow([SV_TARGET if t == "T" else OTHER_TUMOR, NORMAL, caller, os.path.basename(path)])
    expected = {
        "n_union": 6,
        "in_target": {"72354", "72623", "76210", "76667"},
        "not_in_target": {"73459", "81625"},
        "filtered_germline": "82677",
    }
    return ss, files, expected


def build_all(outdir):
    os.makedirs(outdir, exist_ok=True)
    fasta, seqs = _reference(outdir)
    bam, ref1, alt1 = _bam(outdir, seqs)
    snv_ss, snv_vcfs, snv_union = _snv_vcfs(outdir, seqs, alt1)
    ss, indel_vcfs, indel_union = _indel_vcfs(outdir, seqs, snv_ss)
    sv_ss, sv_vcfs, sv_expected = _sv_vcfs(outdir)
    truth = {
        "fasta": fasta, "bam": bam, "target": TARGET, "normal": NORMAL, "other_tumor": OTHER_TUMOR,
        "S1": S1, "S2": S2, "S3": S3, "S1_ref": ref1, "S1_alt": alt1,
        "truth_S1": TRUTH_S1, "truth_S2": TRUTH_S2,
        "snv_samplesheet": snv_ss, "snv_vcfs": snv_vcfs, "snv_union": snv_union,
        "samplesheet": ss, "indel_vcfs": indel_vcfs, "indel_sites": INDEL_SITES,
        "truth_indels": TRUTH_INDELS, "indel_union": indel_union,
        "sv_samplesheet": sv_ss, "sv_vcfs": sv_vcfs, "sv_expected": sv_expected,
        "sv_target": SV_TARGET, "seqs": seqs,
    }
    with open(os.path.join(outdir, "truth.json"), "w") as fh:
        portable = {k: (os.path.basename(v) if isinstance(v, str) and os.path.isabs(v) else v)
                    for k, v in truth.items()
                    if k not in ("seqs", "snv_union", "sv_vcfs", "snv_vcfs", "indel_vcfs")}
        json.dump(portable,
                  fh, indent=2, sort_keys=True, default=lambda o: sorted(o) if isinstance(o, set) else str(o))
    return truth


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("outdir")
    build_all(ap.parse_args().outdir)
