# Parameters

Both entry scripts share the input parameters; each has its own analysis parameters.

## Shared inputs

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | required | CSV with columns `tumor_id,normal_id,caller,vcf` and optional `variant_type` (`snv_genotyping.nf`: `snv` or `indel`, default `snv`; `sv_genotyping.nf`: `sv`, the default). Extra columns such as `analysis_pk,app_name,app_version` are ignored. Relative `vcf` paths resolve against the samplesheet's directory. One row per caller VCF per tumor–normal pair (and variant type). |
| `--target_bam` | required | BAM of the target sample to genotype. |
| `--target_bai` | `<bam>.bai`, then `<bam without .bam>.bai` | BAM index. |
| `--sample_id` | required | Target sample id (`[A-Za-z0-9._-]+`). Prefixes all outputs; a union variant is `called_in_target` when this id is one of its source tumors. |
| `--reference_fasta` | none | Indexed FASTA (`.fai` next to it). SNV: REF check. Indel: REF check, normalization and ref/alt haplotypes — **required** when the samplesheet has `indel` rows. SV: REF bases, svtyper `-T` (needed for clip matching), left-alignment of sequence-resolved insertions and the insertion haplotypes — without it insertions are reported with `genotype_status` `skipped_no_reference`. Strongly recommended. |
| `--contigs` | `1..22,X,Y` | Calls on other contigs are dropped. Empty string = keep all. |
| `--outdir` | `results` | Output directory. |
| `--container` | `mattmyers/nf-wgs-genotyping:1.3.0` | Image for docker/singularity profiles: `[registry/]<namespace>/<image>:<tag>`, no `docker:` prefix (Nextflow adds `docker://` for singularity). |
| `--conda_env` | `environment.yml` | Prebuilt env path for `-profile conda`. |
| `--python` | `python3` | Interpreter that runs the `bin/` scripts. Containers/`-profile conda` need no change. When running directly in an activated env (no container/conda profile), pass an absolute path (`--python $(which python)`): shell start-up files can change `PATH` inside tasks so that `python3` resolves to another interpreter (e.g. `/usr/bin/python3` without pyyaml). |

## SNV and indel (`snv_genotyping.nf`)

Samplesheet callers: `variant_type` `snv` (default): `caveman`, `mutect`, `strelka`, `generic`;
`variant_type` `indel`: `pindel`, `mutect`, `strelka`, `generic` (`generic` also accepts a missing
FILTER `.`). Without `indel` rows no indel step runs and no indel output is written.

| Parameter | Default | Description |
|---|---|---|
| `--snv_chunk_size` | 5000 | Union sites per counting task. |
| `--min_mapq` | 20 | Minimum read mapping quality (SNV and indel). |
| `--min_baseq` | 20 | Minimum base quality: SNV, at the site; indel, bases below it are masked (`N`, score 0) before realignment. |
| `--count_duplicates` | false | Count duplicate-flagged reads (SNV and indel). |
| `--max_ref_mismatch_frac` | 0.01 | Fail an input VCF whose REF disagrees with the FASTA for more than this fraction of SNV (or indel) alleles (wrong assembly). |
| `--indel_max_length` | 50 | Drop indel alleles whose length change is ≥ this (structural variants; the SV pipeline's `--sv_min_size` boundary). |
| `--indel_min_score_diff` | 5 | A fragment supports a haplotype when its summed alignment score beats the second-best haplotype by at least this. 5 = one high-quality discriminating base (match + 1 vs mismatch − 4). |
| `--indel_chunk_size` | 1000 | Union indel alleles per counting task; a locus (alleles genotyped jointly) is never split, so a chunk can exceed this only by one locus. |

### Genotypes (GT/GQ/PL of the SNV and indel tables and VCFs)

GT is **allele presence, not dosage** (`bin/presence_gt.py`): with k alt of n = ref + alt
fragments, `0/0` = k ~ BetaBinomial(n, mean ε, overdispersion ρ) (alt only from error), `1/1` =
n − k ~ BetaBinomial(n, ε, ρ) (ref only from error), `0/1` = k ~ Binomial(n, f) with f uniform on
[`gt_min_vaf`, 1 − `gt_min_vaf`] (both alleles present at any fraction: heterozygous, subclonal or
impure). PL are the Phred-scaled likelihoods; GQ = −10 log10(1 − P(best GT)) with a flat prior,
capped at 99; GT is `./.` below `--gt_min_gq` (low depth and ambiguous counts are not forced into
a genotype; n = 0 gives PL `0,0,0`). Conflict and other fragments are not in AD.

| Parameter | Default | Description |
|---|---|---|
| `--gt_error_rate` | null → 10^(−`min_baseq`/10)/3 | SNV ε: rate of one specific wrong base. The default is the worst case the base-quality filter allows (0.0033 at BQ 20), split over the three wrong bases: conservative. The QC JSON reports the sample's `empirical_error_rate` (other-base fragments / 2 × depth) for comparison. |
| `--indel_gt_error_rate` | 0.01 | Indel ε: rate of fragments supporting one specific wrong haplotype. It cannot be estimated from the counts (`other_count` is another allele of the locus, not an error); raise it for unions rich in short tandem repeats (polymerase slippage). |
| `--gt_overdispersion` | 0.01 | Beta-binomial ρ of the error rate between sites (0 = binomial). |
| `--gt_min_vaf` | 0.05 | Smallest allele fraction treated as present: `0/0` means the alt is not detected at ≥ this fraction. Must exceed the error rate. |
| `--gt_min_gq` | 20 | No-call (`./.`) below this GQ. With the defaults, `0/0` (0 alt) and `1/1` (0 ref) need ≥ 26 fragments (SNV) or ≥ 28 (indel); fewer give `./.`. |

## SV (`sv_genotyping.nf`)

Samplesheet callers: `brass`, `gridss`, `svaba`, `generic`.

| Parameter | Default | Description |
|---|---|---|
| `--sv_slop` | 100 | Clusters (summaries; every adjacency is genotyped on its own): two adjacencies with the same orientation link when both ends are within this distance (bp), and a cluster keeps every pair of members within it at both ends (bounded diameter: no transitive chaining). |
| `--sv_min_size` | 50 | Drop events whose net length change is below this: deleted bases for DEL, duplicated bases for DUP, net inserted bases for insertions (`pos2 − pos1` for INV). With `--indel_max_length` equal (both 50), every insertion/deletion belongs to exactly one workflow. |
| `--sv_imprecise_ci` | 150 | ± CI given to IMPRECISE calls without CIPOS. |
| `--sv_chunk_size` | 500 | Union adjacencies per svtyper task (BND mates never split). |
| `--sv_ins_chunk_size` | 200 | Union insertions per realignment task (`sv_ins_count.py`; needs `--reference_fasta`). Insertions use `--min_mapq`, `--min_baseq`, `--count_duplicates` and `--indel_min_score_diff`. |
| `--sv_rep_priority` | `gridss,brass,svaba` | Caller order for picking a cluster's representative coordinates (after precise > imprecise and sequence-resolved > length-only insertions). |
| `--gridss_max_normal_vf` | 0 | GRIDSS somatic filter: max supporting fragments (FORMAT/VF) in the matched normal column (named `normal_id`). |
| `--brass_missing_filter_is_pass` | true | BRASS writes FILTER `.`; treat it as PASS. |
| `--svtyper_cmd` | `svtyper` | Command used to launch svtyper (e.g. `python -m svtyper.classic`). |
| `--svtyper_args` | `--genotype_on_fragments --keep_all_ref --both_sides --clip_read_support --clip_partner_match --clip_min_base_quality 20` | Extra svtyper flags. Clip flags are dropped (with a warning) when no FASTA is given. |
| `--svtyper_max_reads` | null | svtyper `--max_reads`. Leave unset: when exceeded svtyper emits `./.` with no counts. |
