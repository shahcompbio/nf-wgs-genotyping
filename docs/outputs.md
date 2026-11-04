# Outputs

All files are prefixed with `--sample_id` (`<S>`). Count tables are gzipped CSV with a csverve-style
`<file>.yaml` (`header`, `sep`, `columns: [{name, dtype}]`).

## SNV (`snv_genotyping.nf`)

| File | Content |
|---|---|
| `<S>_snv_genotyping.csv.gz` (+ `.yaml`) | One row per union SNV (columns below). |
| `<S>_genotyping_qc.json` | Union size, `called_in_target`, `positive_control_fraction` (share of the target's own calls with alt ≥ 1), conflict fragments, median depth, per-input extraction stats. |
| `union/<S>_union_snvs.tsv.gz` | Union with provenance. |
| `union/<S>_union_snvs.vcf.gz(.tbi)` | Sites-only union VCF (INFO `CALLERS`, `SAMPLES`, `NCALLERS`, `NSAMPLES`, `INTARGET`). |
| `inputs/*.snvs.stats.json` | Per caller VCF: records, kept, drops by reason. |
| `pipeline_info/` | Nextflow report/timeline/trace and `software_versions.yml`. |

`<S>_snv_genotyping.csv.gz` columns: `sample_id, chromosome, position (1-based), ref, alt,
ref_count, alt_count, other_count` (fragments passing filters, each fragment counted once),
`depth = ref+alt+other`, `vaf = alt/(ref+alt)` (NaN if 0), `n_callers, callers, n_source_samples,
source_samples, sources (T_vs_N:caller;…), called_in_target`.

## SV (`sv_genotyping.nf`)

| File | Content |
|---|---|
| `<S>_sv_genotyping.csv.gz` (+ `.yaml`) | One row per union SV (columns below). |
| `<S>_svtyper.vcf.gz(.tbi)` | svtyper output for all chunks, sorted. |
| `<S>.lib_info.json` | svtyper library (insert size) statistics of the target BAM. |
| `<S>_genotyping_qc.json` | Union size, positive control, `status_counts`, per-input stats. |
| `union/<S>_union_svs.bedpe.gz` | Union SVs (representative coordinates, provenance, cluster spread). |
| `union/<S>_union_svs.members.tsv.gz` | Every normalized input junction with its `sv_id`. |
| `union/<S>_union_svs.vcf.gz(.tbi)` | The svtyper-input VCF of the union. |
| `inputs/*.svs.stats.json` | Per caller VCF: records, kept, drops by reason. |

`<S>_sv_genotyping.csv.gz` columns:
- coordinates: `sample_id, sv_id, chromosome_1, position_1, strand_1, chromosome_2, position_2,
  strand_2, sv_class (DEL/DUP/INV/TRA), svtyper_type (DEL/DUP/BND), sv_length`;
  positions are the retained base next to the junction (VCF 4.2), strand `+` = sequence kept to
  the left of the position, `-` = to the right;
- fragment counts (svtyper `--fragment_counts`, each read pair once): `alt_fragments (AOF),
  ref_fragments (ROF), alt_split_fragments (FS), alt_span_fragments (FP), alt_clip_fragments (FC),
  conflict_fragments (FX), frag_vaf = AOF/(AOF+ROF)`;
- svtyper read-level fields: `svt_gt, svt_gq, svt_sq, svt_dp, svt_ro, svt_ao, svt_rs, svt_as,
  svt_asc, svt_rp, svt_ap, svt_ab`;
- `genotype_status`: `genotyped`, `no_evidence` (no informative fragments), `not_genotyped`
  (record without counts), `missing_from_output`;
- provenance: `n_members, n_callers, callers, n_source_samples, source_samples, sources
  (T_vs_N:caller:record_id;…), called_in_target, spread_1, spread_2, rep_caller, rep_source`.

For BND pairs svtyper genotypes both records from the same fragments; only the `_1` record is
reported (never the sum of both ends).
