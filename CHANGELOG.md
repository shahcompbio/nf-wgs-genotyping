# Changelog

## v1.3.0 (2026-11-05)

- **SNV and indel genotypes.** The tables gain `gt`, `gq`, `pl` (appended; existing columns
  unchanged), and the VCFs gain FORMAT `GT` (first), `GQ`, `PL` plus a
  `##nf-wgs-genotyping_gt_model` line. GT is **allele presence, not dosage**
  (`bin/presence_gt.py`):
  - `0/0`: alt only from error (beta-binomial, mean ε, overdispersion ρ);
  - `1/1`: ref only from error;
  - `0/1`: both alleles present at any fraction (f uniform on [`gt_min_vaf`, 1 − `gt_min_vaf`]);
  - `./.`: GQ < `gt_min_gq`, so low depth and ambiguous counts are no-calls.

  New params:
  - `gt_error_rate` (SNV; default 10^(−min_baseq/10)/3);
  - `indel_gt_error_rate` (0.01);
  - `gt_overdispersion` (0.01);
  - `gt_min_vaf` (0.05);
  - `gt_min_gq` (20).

  The QC JSONs add `gt_model`, `gt_counts`, `called_in_target_gt_counts`, and (SNV) `empirical_error_rate`.
- **SV insertions.** Previously dropped, they are now kept:
  - Sources: `<INS>` and its subtypes (length from SVLEN/SVINSLEN/INSLEN, sequence from
    SVINSSEQ), sequence-resolved REF/ALT alleles, and breakends carrying inserted sequence.
  - Representation: junctions `(c, pos, +)–(c, pos + del_len + 1, −)` with new member columns
    `ins_seq` and `ins_seq_available`. With a FASTA they are left-aligned.
  - Genotyping: realignment to ref/insertion haplotypes (`bin/sv_ins_count.py`, the indel
    machinery; length-only calls use an N block). New module `SV_INS_COUNT` and param
    `sv_ins_chunk_size` (200).
  - Needs `--reference_fasta`; without it the status is `skipped_no_reference`.
  - IMPRECISE insertions are reported, not genotyped (`skipped_imprecise`).
  - `<INS>` without length or sequence is dropped (`insertion_length_unknown`).
- **Every distinct SV adjacency is genotyped.** Previously only each cluster's representative was.
  - New `<S>_sv_adjacency_genotyping.csv.gz` (+ `.yaml`) and `union/<S>_union_sv_adjacencies.tsv.gz`.
  - `adj_id = <sv_id>.<k>` maps each genotype to its source calls (members table gains `adj_id`).
  - svtyper record ids are now `adj_id` (INFO `EVENT=<adj_id>`, new `CLUSTER=<sv_id>`).
  - `<S>_sv_genotyping.csv.gz` keeps one row per cluster with the v1.2 columns, from the
    representative adjacency (`rep_adj_id`). New columns are appended: `rep_adj_id`, `genotyper`,
    `n_adjacencies`, `n_adjacencies_genotyped`, `n_adjacencies_with_alt`, `max_alt_fragments`,
    `max_frag_vaf`, `n_linked_clusters`.
  - QC adds `n_adjacencies`, `adjacency_status_counts`, `positive_control_fraction_any_adjacency`,
    `clusters_with_multiple_adjacencies`, `clusters_linked_to_other_clusters` and `insertions`.
  - svtyper now runs on every adjacency, not one record per cluster: more records per run.
  - New `<S>_sv_adjacency_genotyping.bedpe`: the adjacency table as BEDPE (IGV, bedtools). Each end
    is the 0-based junction base `[position − 1, position)`; `name` = `adj_id`, `score` =
    `alt_fragments`; then the other adjacency columns; missing values `.`.
  - The adjacency table includes each end's confidence interval (`ci1_lo`, `ci1_hi`, `ci2_lo`, `ci2_hi`).
- **Behavior change — SV clusters have a bounded diameter.** Adjacencies are linked within
  `sv_slop` as before, but a cluster now keeps every pair of members within `sv_slop` at both ends.
  Before, union-find connected components could chain distant events into one cluster. The result
  is still independent of input order. `n_linked_clusters` flags clusters with another cluster
  within slop.
- **Behavior change — SV size boundary.** `sv_min_size` now applies to the net length change:
  - DEL: deleted bases (pos2 − pos1 − 1);
  - DUP: duplicated bases (pos2 − pos1 + 1);
  - INS: net inserted bases.

  Before, DEL/DUP used pos2 − pos1. So a 49-bp deletion was kept by both the SV and indel
  workflows, and a 50-bp tandem duplication by neither. With `sv_min_size` = `indel_max_length` =
  50, every insertion/deletion now belongs to exactly one workflow.
- New SV genotype statuses: `skipped_imprecise` and `skipped_no_reference`.
- Image `mattmyers/nf-wgs-genotyping:1.3.0`: same environment, `environment.yml` pins
  `scipy>=1.4`, and the Docker check imports scipy. `params.container` updated.
- Fixtures:
  - Contig 4 of the synthetic reference holds designed SV insertions, with their reads in
    `target.bam` (mates 230–270 bp away, so svtyper can estimate the insert-size spread) and
    `sv_ins_samplesheet.csv`.
  - Contigs 1–3, their reads and the SNV/indel truth are unchanged.
- Internal:
  - `indel_count` split into `locus_reads`, `haplotypes` and `assign_fragments`; behavior unchanged.
  - `snv_finalize.write_vcf` takes the genotype model.
  - `sv_union.cluster_union` returns `(clusters, adjacencies, members)`.
  - `sv_finalize.finalize` takes the adjacency table and the insertion counts and manifest.

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
