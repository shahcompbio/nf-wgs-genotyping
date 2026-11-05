#!/usr/bin/env python3
"""Union extracted SNVs across callers and source samples, keep provenance, and split the union
into chunks for scatter/gather.

Outputs (in --outdir):
  <sample>_union_snvs.tsv.gz       one row per (chromosome, position, ref, alt) with provenance
  <sample>_union_snvs.vcf.gz(.tbi) sites-only VCF (INFO: CALLERS, SAMPLES, NCALLERS, NSAMPLES, INTARGET)
  chunks/chunk_NNNNN.tsv           chromosome, position, ref, alt (always >= 1 chunk, may be header-only)
  chunks.manifest                  chunk_id, n_items (used by snv_finalize to verify the gather)
"""
import argparse
import os
import sys

import pandas as pd
import pysam

import geno_utils as gu

KEYS = ["chromosome", "position", "ref", "alt"]
UNION_COLUMNS = KEYS + ["n_callers", "callers", "n_source_samples", "source_samples", "sources",
                        "called_in_target"]


def build_union(tsvs, sample_id, order):
    frames = [pd.read_csv(t, sep="\t", dtype={"chrom": str, "tumor_id": str, "normal_id": str})
              for t in tsvs]
    calls = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if calls.empty:
        return pd.DataFrame({c: pd.Series(dtype=object) for c in UNION_COLUMNS})
    calls = calls.rename(columns={"chrom": "chromosome", "pos": "position"})
    calls["source"] = calls["tumor_id"] + "_vs_" + calls["normal_id"] + ":" + calls["caller"]
    rows = []
    for key, g in calls.groupby(KEYS, sort=False):
        callers = sorted(set(g["caller"]))
        samples = sorted(set(g["tumor_id"]))
        rows.append(dict(zip(KEYS, key), n_callers=len(callers), callers=";".join(callers),
                         n_source_samples=len(samples), source_samples=";".join(samples),
                         sources=";".join(sorted(set(g["source"]))), called_in_target=sample_id in samples))
    key = gu.site_sort_key(order)
    rows.sort(key=lambda r: (key(r["chromosome"], r["position"]), r["ref"], r["alt"]))
    return pd.DataFrame(rows, columns=UNION_COLUMNS)


def _write_vcf(df, path_vcf, order, fai=None, source="snv_union"):
    lengths = {}
    if fai:
        with open(fai) as fh:
            lengths = {l.split("\t")[0]: l.split("\t")[1] for l in fh if l.strip()}
    contigs = list(order) + sorted(set(df["chromosome"]) - set(order))
    with open(path_vcf, "w") as fh:
        fh.write(f"##fileformat=VCFv4.2\n##source=nf-wgs-genotyping {source}\n")
        for c in contigs:
            fh.write(f"##contig=<ID={c}" + (f",length={lengths[c]}" if c in lengths else "") + ">\n")
        fh.write('##INFO=<ID=CALLERS,Number=.,Type=String,Description="Callers with a PASS call">\n'
                 '##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Tumor samples with a PASS call">\n'
                 '##INFO=<ID=NCALLERS,Number=1,Type=Integer,Description="Number of callers">\n'
                 '##INFO=<ID=NSAMPLES,Number=1,Type=Integer,Description="Number of source samples">\n'
                 '##INFO=<ID=INTARGET,Number=0,Type=Flag,Description="Called in the target sample">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in df.itertuples(index=False):
            info = (f"CALLERS={r.callers.replace(';', ',')};SAMPLES={r.source_samples.replace(';', ',')};"
                    f"NCALLERS={r.n_callers};NSAMPLES={r.n_source_samples}"
                    + (";INTARGET" if r.called_in_target else ""))
            fh.write(f"{r.chromosome}\t{r.position}\t.\t{r.ref}\t{r.alt}\t.\tPASS\t{info}\n")
    pysam.tabix_index(path_vcf, preset="vcf", force=True)


def write_outputs(df, sample_id, outdir, chunk_size, contigs_order, fai=None):
    os.makedirs(os.path.join(outdir, "chunks"), exist_ok=True)
    df.to_csv(os.path.join(outdir, f"{sample_id}_union_snvs.tsv.gz"), sep="\t", index=False,
              compression="gzip")
    _write_vcf(df, os.path.join(outdir, f"{sample_id}_union_snvs.vcf"), contigs_order, fai)
    chunks = []
    n_chunks = max(1, -(-len(df) // chunk_size))
    for i in range(n_chunks):
        part = df.iloc[i * chunk_size:(i + 1) * chunk_size][KEYS]
        cid = f"chunk_{i:05d}"
        part.to_csv(os.path.join(outdir, "chunks", f"{cid}.tsv"), sep="\t", index=False)
        chunks.append((cid, len(part)))
    gu.write_manifest(os.path.join(outdir, "chunks.manifest"), chunks)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", nargs="+", required=True, help="snv_extract TSVs")
    ap.add_argument("--sample-id", required=True, help="target sample id (sets called_in_target)")
    ap.add_argument("--reference", default=None, help="FASTA whose .fai sets contig order")
    ap.add_argument("--contigs", default=None, help="contig order when no reference is given")
    ap.add_argument("--chunk-size", type=int, default=5000)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)
    if a.chunk_size < 1:
        ap.error("--chunk-size must be >= 1")
    fai = f"{a.reference}.fai" if a.reference else None
    order = gu.contig_order(fai=fai, contigs=gu.parse_contigs(a.contigs))
    df = build_union(a.inputs, a.sample_id, order)
    write_outputs(df, a.sample_id, a.outdir, a.chunk_size, order, fai)
    return 0


if __name__ == "__main__":
    sys.exit(main())
