# Outputs

All files are prefixed with `--sample_id` (`<S>`). Count tables are gzipped CSV with a csverve-style
`<file>.yaml` (`header`, `sep`, `columns: [{name, dtype}]`).

## SNV (`snv_genotyping.nf`)

| File | Content |
|---|---|
| `<S>_snv_genotyping.csv.gz` (+ `.yaml`) | One row per union SNV (columns below). |
| `<S>_snv_genotyping.vcf.gz(.tbi)` | The same rows as a single-sample VCF (below). |
| `<S>_genotyping_qc.json` | Union size, `called_in_target`, `positive_control_fraction` (share of the target's own calls with alt ≥ 1), conflict fragments, median depth, `empirical_error_rate` (other-base fragments / (2 × depth)), `gt_model` (parameters), `gt_counts` and `called_in_target_gt_counts` (per GT), per-input extraction stats. |
| `union/<S>_union_snvs.tsv.gz` | Union with provenance. |
| `union/<S>_union_snvs.vcf.gz(.tbi)` | Sites-only union VCF (INFO `CALLERS`, `SAMPLES`, `NCALLERS`, `NSAMPLES`, `INTARGET`). |
| `inputs/*.snvs.stats.json` | Per caller VCF: records, kept, drops by reason. |
| `pipeline_info/` | Nextflow report/timeline/trace and `software_versions.yml`. |

`<S>_snv_genotyping.csv.gz` columns: `sample_id, chromosome, position (1-based), ref, alt,
ref_count, alt_count, other_count` (fragments passing filters, each fragment counted once),
`depth = ref+alt+other`, `vaf = alt/(ref+alt)` (NaN if 0), `n_callers, callers, n_source_samples,
source_samples, sources (T_vs_N:caller;…), called_in_target`, then the genotype `gt` (`0/0`,
`0/1`, `1/1` or `./.`), `gq` and `pl` (`"PL00,PL01,PL11"`).

**GT is allele presence, not dosage** (model and parameters in [params.md](params.md#genotypes-gtgqpl-of-the-snv-and-indel-tables-and-vcfs)):
`0/0` = only the reference allele is detected (alt not above sequencing error at ≥ `gt_min_vaf`),
`0/1` = both alleles are present at any fraction (heterozygous, subclonal, or diluted by normal
cells), `1/1` = only the alt allele is detected (e.g. LOH), `./.` = not enough evidence
(GQ < `gt_min_gq`; always for sites with no ref/alt fragments). A tumor `0/1` at VAF 0.1 is a
present subclonal or impure variant, not a heterozygous germline genotype; use `vaf` and copy
number for dosage.

`<S>_snv_genotyping.vcf.gz`: one biallelic record per CSV row, same order (sorted, tabix-indexed),
FILTER `PASS`, sample column `<S>`. Header = union VCF header (contigs with lengths) plus:
- INFO `CALLERS`, `SAMPLES` (source tumors), `NCALLERS`, `NSAMPLES`, `SOURCES`
  (`T_vs_N:caller`), flag `INTARGET` (= `called_in_target`);
- FORMAT `GT` (first), `AD` = `ref_count,alt_count`, `DP` = `depth`, `AF` = `vaf` (`.` when
  ref+alt = 0), `OC` = `other_count`, `GQ`, `PL` (Number=G: 0/0, 0/1, 1/1);
- a `##nf-wgs-genotyping_gt_model=<model=allele_presence,error_rate=…,overdispersion=…,min_vaf=…,min_gq=…>`
  line recording the genotype model parameters.

## Indel (`snv_genotyping.nf`, samplesheet rows with `variant_type` `indel`)

Written only when the samplesheet has `indel` rows.

| File | Content |
|---|---|
| `<S>_indel_genotyping.csv.gz` (+ `.yaml`) | One row per union indel allele (columns below). |
| `<S>_indel_genotyping.vcf.gz(.tbi)` | The same rows as a single-sample VCF (below). |
| `<S>_indel_qc.json` | Union size, `loci`, `called_in_target`, `positive_control_fraction`, uninformative and conflict fragments (summed once per locus), median depth, `gt_model`, `gt_counts`, `called_in_target_gt_counts`, per-input stats. |
| `union/<S>_union_indels.tsv.gz` | Union with locus and provenance columns. |
| `union/<S>_union_indels.vcf.gz(.tbi)` | Sites-only union VCF (INFO as for SNVs). |
| `inputs/*.indels.stats.json` | Per caller VCF: records, kept, `normalized` (representation changed by left-alignment/trimming), drops by reason (`not_pass, not_indel, mnv, too_long, off_contig, ref_mismatch, duplicate`). |

Alleles are normalized (left-aligned, trimmed; a pure insertion/deletion keeps its shared first
base), so the same event from different callers is one row.

`<S>_indel_genotyping.csv.gz` columns: `sample_id, chromosome, position (1-based), ref, alt,
indel_length` (len(alt) − len(ref): + insertion, − deletion), `ref_count, alt_count, other_count`
(fragments, each counted once, supporting the reference, this allele, or another union allele of the
same locus), `depth = ref+alt+other`, `vaf = alt/(ref+alt)` (NaN if 0), `uninformative_count`
(fragments reaching the locus that support no haplotype by `--indel_min_score_diff`, e.g. ending
inside the repeat), `conflict_count` (mates supporting different haplotypes),
`locus_id` (`chrom:start-end`), `var_start, var_end` (1-based reference bases a read must reach:
every equivalent placement of the event plus one flanking base each side), then the provenance
columns of the SNV table, then `gt, gq, pl` (as for SNVs, from `ref_count`/`alt_count` with
`--indel_gt_error_rate`). `uninformative_count` and `conflict_count` are per locus: alleles of one
locus share them.

`<S>_indel_genotyping.vcf.gz`: as the SNV VCF (GT/GQ/PL included), plus INFO `LOCUS` (= `locus_id`)
and FORMAT `UF` (= `uninformative_count`); `OC` = `other_count` (other alleles of the locus).

## SV (`sv_genotyping.nf`)

| File | Content |
|---|---|
| `<S>_sv_adjacency_genotyping.csv.gz` (+ `.yaml`) | One row per distinct union **adjacency**: every adjacency is genotyped (columns below). |
| `<S>_sv_adjacency_genotyping.bedpe` | The adjacency table as BEDPE: one row per adjacency, 0-based junction bases (below). |
| `<S>_sv_genotyping.csv.gz` (+ `.yaml`) | One row per union SV **cluster**: the v1.2 columns from the cluster's representative adjacency, then a summary of all its adjacencies. |
| `<S>_svtyper.vcf.gz(.tbi)` | svtyper output for all chunks, sorted (record ids = `adj_id`). |
| `<S>.lib_info.json` | svtyper library (insert size) statistics of the target BAM. |
| `<S>_genotyping_qc.json` | Union size, positive control (representative and any adjacency), `status_counts` (clusters), `n_adjacencies`, `adjacency_status_counts`, `clusters_with_multiple_adjacencies`, `clusters_linked_to_other_clusters`, `insertions` (counts by kind and status), per-input stats. |
| `union/<S>_union_svs.bedpe.gz` | Union clusters (representative coordinates, provenance, spread, `rep_adj_id`, `n_adjacencies`, `n_linked_clusters`). |
| `union/<S>_union_sv_adjacencies.tsv.gz` | Union adjacencies with provenance and, for insertions, `ins_seq`. |
| `union/<S>_union_svs.members.tsv.gz` | Every normalized input junction with its `sv_id` and `adj_id` (`ins_len`, `ins_seq`, `ins_seq_available` for insertions). |
| `union/<S>_union_svs.vcf.gz(.tbi)` | The svtyper-input VCF of the union (non-insertion adjacencies). |
| `inputs/*.svs.stats.json` | Per caller VCF: records, kept, drops by reason (incl. `insertion_length_unknown`), `insertions` (kept, sequence-resolved, length-only, left-aligned). |

**Adjacencies and clusters.** An adjacency is a distinct normalized junction: identical ends
(and, for insertions, identical inserted sequence or length). Calls of the same adjacency by
several callers or samples are its members. Each adjacency is genotyped on its own (svtyper, or
realignment for insertions), so a cluster of nearby but different breakpoints reports every one
of them; `adj_id = <sv_id>.<k>` maps each genotype to its source calls (`sources`, and the
members table). Clusters group adjacencies whose members are all within `--sv_slop` of each
other (bounded diameter); `n_linked_clusters` > 0 flags clusters with another cluster within slop
(worth a look).

**Insertions** (`sv_class` `INS`, `genotyper` `realign`): written as junctions
`(c, pos, +)–(c, pos + del_len + 1, −)`, pos = the reference base left of the insertion (usually
`del_len` = 0, so `position_2 = position_1 + 1`); `sv_length` = net inserted bases; `ins_len` =
inserted bases; `ins_seq_available` false for length-only calls. svtyper skips SVTYPE=INS, so
insertions are genotyped by realigning reads to the reference and insertion haplotypes (as for
indels; an N block for length-only calls): `alt_fragments`, `ref_fragments`,
`uninformative_fragments`, `conflict_fragments`; the svtyper columns are empty.

`<S>_sv_adjacency_genotyping.csv.gz` columns: `sample_id, adj_id, sv_id`, the coordinate columns
below, `sv_class, genotyper (svtyper|realign), svtyper_type, sv_length, ins_len,
ins_seq_available, imprecise, ci1_lo, ci1_hi, ci2_lo, ci2_hi` (confidence interval of each end,
offsets from the position, from the adjacency's best member), the fragment counts (+ `uninformative_fragments`, insertions
only), `frag_vaf`, the svtyper read-level fields, `genotype_status`, `is_cluster_rep`, and
`n_members, n_callers, callers, n_source_samples, source_samples, sources
(T_vs_N:caller:record_id;…), called_in_target` of the adjacency.

`<S>_sv_adjacency_genotyping.bedpe` (plain text, tab-separated, for IGV, bedtools `pairtobed` /
`pairtopair` and similar): the rows of `<S>_sv_adjacency_genotyping.csv.gz` in BEDPE layout.
- Columns 1–10: `chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2`.
  - Each end is the junction base as a 0-based half-open interval: `start = position − 1`,
    `end = position`, so `end1`/`end2` equal `position_1`/`position_2`. An insertion is
    `[pos − 1, pos)` and `[pos + del_len, pos + del_len + 1)`.
  - The interval does not include the confidence interval; that is in `ci1_lo` … `ci2_hi`.
  - `name` = `adj_id`; `score` = `alt_fragments`.
  - Strands as in the tables: DEL `+ −`, DUP `− +` (the LUMPY/svtools convention).
- Then every other column of the adjacency table, in the same order.
- The header line starts with `#`. Missing values (e.g. `score` of an adjacency without counts)
  are `.`.
- Rows are sorted by reference order of (`chrom1`, `start1`), then (`chrom2`, `start2`), `name`.

`<S>_sv_genotyping.csv.gz` columns (cluster table):
- coordinates: `sample_id, sv_id, chromosome_1, position_1, strand_1, chromosome_2, position_2,
  strand_2, sv_class (DEL/DUP/INV/TRA/INS), svtyper_type (DEL/DUP/BND; empty for INS), sv_length`;
  positions are the retained base next to the junction (VCF 4.2), strand `+` = sequence kept to
  the left of the position, `-` = to the right;
- fragment counts (svtyper `--fragment_counts`, each read pair once): `alt_fragments (AOF),
  ref_fragments (ROF), alt_split_fragments (FS), alt_span_fragments (FP), alt_clip_fragments (FC),
  conflict_fragments (FX), frag_vaf = AOF/(AOF+ROF)`;
- svtyper read-level fields: `svt_gt, svt_gq, svt_sq, svt_dp, svt_ro, svt_ao, svt_rs, svt_as,
  svt_asc, svt_rp, svt_ap, svt_ab`;
- `genotype_status`: `genotyped`, `no_evidence` (no informative fragments), `not_genotyped`
  (record without counts), `missing_from_output`, `skipped_imprecise` (IMPRECISE insertion:
  position too uncertain to build a haplotype), `skipped_no_reference` (insertion without
  `--reference_fasta`);
- provenance: `n_members, n_callers, callers, n_source_samples, source_samples, sources
  (T_vs_N:caller:record_id;…), called_in_target, spread_1, spread_2, rep_caller, rep_source`;
- v1.3 (all adjacencies of the cluster): `rep_adj_id` (the adjacency the columns above come
  from), `genotyper`, `n_adjacencies`, `n_adjacencies_genotyped`, `n_adjacencies_with_alt`
  (alt ≥ 1), `max_alt_fragments`, `max_frag_vaf`, `n_linked_clusters`.

The counts, `svt_*` fields and `genotype_status` of a cluster row are those of its representative
adjacency, as in v1.2 (when a cluster has one adjacency, nothing changes); a cluster supported
only at another of its breakpoints shows `n_adjacencies_with_alt` ≥ 1 with `alt_fragments` 0.

For BND pairs svtyper genotypes both records from the same fragments; only the `_1` record is
reported (never the sum of both ends).
