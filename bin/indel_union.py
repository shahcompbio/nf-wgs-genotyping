#!/usr/bin/env python3
"""Union extracted (normalized) indels across callers and source samples, group alleles into
loci, and split the union into locus-whole chunks for scatter/gather.

A locus is a set of alleles whose footprints (the reference span over which the event can slide,
plus one flanking base each side; see indel_utils.variable_span) share at least two bases, e.g.
a 1-bp and a 2-bp deletion in the same homopolymer. indel_count genotypes a locus's alleles
jointly, so a chunk never splits a locus.

Outputs (in --outdir):
  <sample>_union_indels.tsv.gz       one row per (chromosome, position, ref, alt): locus, provenance
  <sample>_union_indels.vcf.gz(.tbi) sites-only VCF (INFO: CALLERS, SAMPLES, NCALLERS, NSAMPLES, INTARGET)
  chunks/chunk_NNNNN.tsv             chromosome, position, ref, alt, locus_id (>= 1 chunk)
  chunks.manifest                    chunk_id, n_items (used by indel_finalize to verify the gather)
"""
import argparse
import os
import sys

import pandas as pd
import pysam

import geno_utils as gu
import indel_utils as iu
import snv_union

KEYS = ["chromosome", "position", "ref", "alt"]
LOCUS_COLUMNS = ["indel_length", "locus_id", "var_start", "var_end"]
UNION_COLUMNS = KEYS + LOCUS_COLUMNS + ["n_callers", "callers", "n_source_samples", "source_samples",
                                        "sources", "called_in_target"]
CHUNK_COLUMNS = KEYS + ["locus_id"]


def build_union(tsvs, sample_id, order):
    """Group extracted alleles by (chromosome, position, ref, alt) with provenance (as SNVs)."""
    return snv_union.build_union(tsvs, sample_id, order)


def _merge(intervals):
    """intervals: [(start0, end0)] -> list of merged locus intervals, one per input (same order).
    Intervals join a locus when they share >= 2 bases with it."""
    order = sorted(range(len(intervals)), key=lambda i: intervals[i])
    out = [None] * len(intervals)
    cur = None
    members = []
    for i in order:
        s, e = intervals[i]
        if cur is not None and s <= cur[1] - 2:
            cur = (cur[0], max(cur[1], e))
            members.append(i)
            continue
        for m in members:
            out[m] = cur
        cur, members = (s, e), [i]
    for m in members:
        out[m] = cur
    return out


def annotate_loci(df, fasta):
    """Add indel_length (len(alt) - len(ref)), var_start/var_end (1-based, inclusive footprint)
    and locus_id ('<chrom>:<start>-<end>' of the merged footprints)."""
    df = df.copy()
    if df.empty:
        for c in LOCUS_COLUMNS:
            df[c] = pd.Series(dtype=object)
        return df[UNION_COLUMNS]
    fa = pysam.FastaFile(fasta)
    lengths = dict(zip(fa.references, fa.lengths))
    fps = []
    for chrom, pos, ref, alt in df[KEYS].itertuples(index=False):
        chrom = str(chrom)
        fetch = lambda a, b, c=chrom: fa.fetch(c, max(a, 0), b)       # noqa: E731
        span = iu.variable_span(fetch, lengths[chrom], int(pos), ref, alt)
        fps.append(iu.footprint(span, lengths[chrom]))
    loci = [None] * len(df)
    chroms = df["chromosome"].astype(str).tolist()
    for chrom in dict.fromkeys(chroms):
        idx = [i for i, c in enumerate(chroms) if c == chrom]
        for i, locus in zip(idx, _merge([fps[i] for i in idx])):
            loci[i] = f"{chrom}:{locus[0] + 1}-{locus[1]}"
    df["indel_length"] = df["alt"].str.len() - df["ref"].str.len()
    df["var_start"] = [f[0] + 1 for f in fps]
    df["var_end"] = [f[1] for f in fps]
    df["locus_id"] = loci
    return df[UNION_COLUMNS]


def write_outputs(df, sample_id, outdir, chunk_size, contigs_order, fai=None):
    os.makedirs(os.path.join(outdir, "chunks"), exist_ok=True)
    df.to_csv(os.path.join(outdir, f"{sample_id}_union_indels.tsv.gz"), sep="\t", index=False,
              compression="gzip")
    snv_union._write_vcf(df, os.path.join(outdir, f"{sample_id}_union_indels.vcf"), contigs_order, fai,
                         source="indel_union")
    groups = [g for _, g in df.groupby("locus_id", sort=False)] if len(df) else []
    parts, current = [], []
    for g in groups:
        if current and sum(len(x) for x in current) + len(g) > chunk_size:
            parts.append(current)
            current = []
        current.append(g)
    if current or not parts:
        parts.append(current)
    chunks = []
    for i, part in enumerate(parts):
        cid = f"chunk_{i:05d}"
        frame = pd.concat(part)[CHUNK_COLUMNS] if part else pd.DataFrame(columns=CHUNK_COLUMNS)
        frame.to_csv(os.path.join(outdir, "chunks", f"{cid}.tsv"), sep="\t", index=False)
        chunks.append((cid, len(frame)))
    gu.write_manifest(os.path.join(outdir, "chunks.manifest"), chunks)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", nargs="+", required=True, help="indel_extract TSVs")
    ap.add_argument("--sample-id", required=True, help="target sample id (sets called_in_target)")
    ap.add_argument("--reference", required=True, help="indexed FASTA (repeat spans; .fai sets contig order)")
    ap.add_argument("--chunk-size", type=int, default=1000, help="max alleles per chunk (a locus is never split)")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)
    if a.chunk_size < 1:
        ap.error("--chunk-size must be >= 1")
    fai = f"{a.reference}.fai"
    order = gu.contig_order(fai=fai)
    df = annotate_loci(build_union(a.inputs, a.sample_id, order), a.reference)
    write_outputs(df, a.sample_id, a.outdir, a.chunk_size, order, fai)
    return 0


if __name__ == "__main__":
    sys.exit(main())
