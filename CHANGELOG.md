# Changelog

## v1.2.0 (2026-11-04)

- Indels in `snv_genotyping.nf`: samplesheet rows with the new optional column `variant_type`
  `indel` (callers `pindel`, `mutect`, `strelka`, `generic`) are extracted (PASS; multiallelic
  split; plain-base alleles of different lengths, complex del-ins included; MNVs and changes
  ≥ `--indel_max_length` 50 bp dropped), left-aligned and trimmed against the FASTA, unioned with
  provenance, grouped into loci (alleles whose repeat-aware footprints overlap) and genotyped by
  semi-global realignment of every fragment to the reference and allele haplotypes (Biopython
  `PairwiseAligner`, BWA-MEM scores; soft clips included; BQ < `--min_baseq` masked;
  `--indel_min_score_diff` 5). New outputs `<S>_indel_genotyping.csv.gz(.yaml)`,
  `<S>_indel_genotyping.vcf.gz(.tbi)`, `<S>_indel_qc.json`, `union/<S>_union_indels.*`,
  `inputs/*.indels.stats.json`. New params `indel_max_length`, `indel_min_score_diff`,
  `indel_chunk_size`. `--reference_fasta` is required when there are indel rows.
- Without `variant_type` (v1.1 samplesheets) every row is `snv`: SNV outputs are unchanged and no
  indel step runs. `sv_genotyping.nf` accepts `variant_type` `sv` only (the default).
- Image `mattmyers/nf-wgs-genotyping:1.2.0` (adds `biopython=1.85`); `params.container` updated.
- Internal: `geno_utils.iter_pass_alleles` (record filtering shared by SNV and indel extraction),
  `snv_finalize.join_counts` / `write_vcf(header_lines, extra)` reused by `indel_finalize.py`;
  test fixtures gain contig 3 with designed indel loci (contigs 1–2 and SNV truth unchanged).

## v1.1.0 (2026-11-04)

- SNV: `<S>_snv_genotyping.vcf.gz` (+ `.tbi`) next to the CSV — one record per CSV row, same order.
  Header from the union VCF (contigs with lengths; INFO `CALLERS`, `SAMPLES`, `NCALLERS`,
  `NSAMPLES`, `INTARGET`) plus INFO `SOURCES`; FORMAT `AD` (ref,alt fragments), `DP`
  (ref+alt+other), `AF` (alt/(ref+alt), `.` when 0), `OC` (other). No `GT`: the pipeline counts
  fragments and makes no genotype call. `snv_finalize.py` takes `--union-vcf` (required).
- SV outputs unchanged. Image unchanged: `params.container` stays
  `mattmyers/nf-wgs-genotyping:1.0.1` (`bin/` is staged by Nextflow; pysam is already in the image).

## v1.0.1 (2026-11-04)

- `params.container = 'mattmyers/nf-wgs-genotyping:1.0.1'` (Docker Hub; Nextflow adds `docker://`
  for singularity). Record the pushed digest here: `mattmyers/nf-wgs-genotyping:1.0.1@sha256:…`.
- `--python`: interpreter for `bin/` scripts (default `python3`); `procps` comes from apt in the
  Dockerfile so `environment.yml` also solves on macOS.
- Test: `params.container` must be a pullable `[registry/]<namespace>/<image>:<tag>` reference.

## v1.0.0 (2026-11-04)

- Two entry points: `snv_genotyping.nf` (PASS SNV union, fragment-level pysam counting) and
  `sv_genotyping.nf` (BRASS/GRIDSS/SvABA PASS SV union, svtyper fork @27a0427).
- GRIDSS somatic filter on matched-normal VF; BRASS missing FILTER treated as PASS.
- Explicit, manifest-checked scatter/gather; csverve-style output YAML; QC JSON with a
  positive control (target's own calls).
- Profiles: iris, slurm, singularity, docker, conda, test, test_snv, test_sv, stub.
- `params.container` was a placeholder; this tag cannot pull an image on iris.
