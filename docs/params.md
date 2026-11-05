# Parameters

Both entry scripts share the input parameters; each has its own analysis parameters.

## Shared inputs

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | required | CSV with columns `tumor_id,normal_id,caller,vcf` and optional `variant_type` (`snv_genotyping.nf`: `snv` or `indel`, default `snv`; `sv_genotyping.nf`: `sv`, the default). Extra columns such as `analysis_pk,app_name,app_version` are ignored. Relative `vcf` paths resolve against the samplesheet's directory. One row per caller VCF per tumor–normal pair (and variant type). |
| `--target_bam` | required | BAM of the target sample to genotype. |
| `--target_bai` | `<bam>.bai`, then `<bam without .bam>.bai` | BAM index. |
| `--sample_id` | required | Target sample id (`[A-Za-z0-9._-]+`). Prefixes all outputs; a union variant is `called_in_target` when this id is one of its source tumors. |
| `--reference_fasta` | none | Indexed FASTA (`.fai` next to it). SNV: REF check. Indel: REF check, normalization and ref/alt haplotypes — **required** when the samplesheet has `indel` rows. SV: REF bases and svtyper `-T` (needed for clip matching). Strongly recommended. |
| `--contigs` | `1..22,X,Y` | Calls on other contigs are dropped. Empty string = keep all. |
| `--outdir` | `results` | Output directory. |
| `--container` | `mattmyers/nf-wgs-genotyping:1.0.1` | Image for docker/singularity profiles: `[registry/]<namespace>/<image>:<tag>`, no `docker:` prefix (Nextflow adds `docker://` for singularity). |
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

## SV (`sv_genotyping.nf`)

Samplesheet callers: `brass`, `gridss`, `svaba`, `generic`.

| Parameter | Default | Description |
|---|---|---|
| `--sv_slop` | 100 | Max distance (bp) at both ends for two junctions with the same orientation to cluster. |
| `--sv_min_size` | 50 | Drop intrachromosomal junctions shorter than this. |
| `--sv_imprecise_ci` | 150 | ± CI given to IMPRECISE calls without CIPOS. |
| `--sv_chunk_size` | 500 | Union SVs per svtyper task (BND mates never split). |
| `--sv_rep_priority` | `gridss,brass,svaba` | Caller order for picking a cluster's representative coordinates (after precise > imprecise). |
| `--gridss_max_normal_vf` | 0 | GRIDSS somatic filter: max supporting fragments (FORMAT/VF) in the matched normal column (named `normal_id`). |
| `--brass_missing_filter_is_pass` | true | BRASS writes FILTER `.`; treat it as PASS. |
| `--svtyper_cmd` | `svtyper` | Command used to launch svtyper (e.g. `python -m svtyper.classic`). |
| `--svtyper_args` | `--genotype_on_fragments --keep_all_ref --both_sides --clip_read_support --clip_partner_match --clip_min_base_quality 20` | Extra svtyper flags. Clip flags are dropped (with a warning) when no FASTA is given. |
| `--svtyper_max_reads` | null | svtyper `--max_reads`. Leave unset: when exceeded svtyper emits `./.` with no counts. |
