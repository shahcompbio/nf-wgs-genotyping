# nf-wgs-genotyping

Genotype the **union of a patient's somatic calls** in one target bulk WGS sample. Two pipelines
share this repository, with distinct entry points:

| Entry point | Input calls | Genotyper | Output |
|---|---|---|---|
| `snv_genotyping.nf` | PASS SNVs from CaVEMan / Mutect / Strelka (or generic) VCFs, and PASS indels from Pindel / Strelka / Mutect (or generic) VCFs, of every tumor–normal pair | SNV: fragment-level pysam pileup; indel: fragment realignment to ref/alt haplotypes | ref/alt/other fragment counts, VAF, allele-presence GT/GQ/PL (separate SNV and indel tables + VCFs) |
| `sv_genotyping.nf` | PASS SVs (incl. insertions) from BRASS / GRIDSS / SvABA (or generic) VCFs of every tumor–normal pair | [svtyper fork](https://github.com/marcjwilliams1/svtyper) @ `27a0427`; insertions: fragment realignment | fragment + read counts and status of every distinct adjacency, summarized per cluster |

Built to be launched by the isabl apps `SNV_GENOTYPING` and `SV_GENOTYPING`
(papaemmelab `isabl_apps`), but runs standalone. The design follows the Shah-lab Mondrian
SNV-genotyping workflow (union of patient calls → scatter → gather) for bulk WGS, with explicit
read filters and an explicit, manifest-checked gather.

## Running

```bash
# from GitHub (pinned release); -main-script selects the pipeline
nextflow run shahcompbio/nf-wgs-genotyping -revision v1.3.0 -main-script snv_genotyping.nf -profile iris \
    --samplesheet samplesheet.csv --target_bam T1.bam --sample_id T1 \
    --reference_fasta /path/GRCh37.fasta --outdir results

# from a local checkout: run the entry script by path (-main-script only resolves inside remote repos)
nextflow run /path/to/nf-wgs-genotyping/sv_genotyping.nf -profile conda ...
```

`main.nf` deliberately refuses to run and prints these two commands.

**Samplesheet** (`tumor_id,normal_id,caller,vcf`, optional `variant_type`; extra columns are
ignored, relative VCF paths resolve against the samplesheet's directory). In `snv_genotyping.nf`,
`variant_type` is `snv` (default when the column is absent or blank) or `indel`; indel rows need
`--reference_fasta`:

```csv
tumor_id,normal_id,caller,vcf,variant_type
T1,N1,caveman,/data/T1_vs_N1/caveman.flagged.vcf.gz,snv
T1,N1,mutect,/data/T1_vs_N1/mutect.snvs.vcf.gz,snv
T1,N1,mutect,/data/T1_vs_N1/mutect.indels.vcf.gz,indel
T1,N1,pindel,/data/T1_vs_N1/pindel.flagged.vcf.gz,indel
T2,N1,strelka,/data/T2_vs_N1/strelka.snvs.vcf.gz,snv
T2,N1,strelka,/data/T2_vs_N1/strelka.indels.vcf.gz,indel
```

Profiles: `iris` (SLURM `componc_cpu` + Singularity), `slurm`, `singularity`, `docker`, `conda`,
`test` (local, small resources), `test_snv` / `test_sv` (committed test data), `stub`.
See [docs/params.md](docs/params.md) and [docs/outputs.md](docs/outputs.md).

## Method

### SNV
1. **Extract** (per VCF): FILTER exactly `PASS`; multiallelic records split; keep single-base
   ACGT REF≠ALT (indels/MNVs dropped); contigs restricted; with a FASTA, REF must match (an input
   with > `max_ref_mismatch_frac` mismatches fails: wrong assembly).
2. **Union**: group by (chrom, pos, ref, alt) keeping provenance (callers, source samples,
   `T_vs_N:caller`, `called_in_target`); sort by FASTA order; split into chunks + manifest.
3. **Count** (scatter): pysam pileup with *no implicit filters* (`stepper=nofilter`,
   `ignore_orphans=False`, `ignore_overlaps=False`, `min_base_quality=0`), then explicit filters:
   drop unmapped, secondary, supplementary, QC-fail, duplicate (unless `--count_duplicates`),
   MAPQ < 20, deletions/ref-skips, base quality < 20. **Each fragment is counted once**:
   overlapping mates that agree count once; disagreeing mates resolve to the higher-quality base;
   a quality tie is a conflict (not counted, reported in QC).
4. **Finalize** (gather): fails unless the gathered chunk ids equal the manifest and every union
   site has exactly one counts row.
5. **Genotype** (`bin/presence_gt.py`): GT is **allele presence, not dosage** — tumor purity, copy
   number and clonality move a somatic variant's VAF anywhere, so the diploid 0, ½, 1 model would
   call subclonal variants `0/0`. With k alt of n = ref + alt fragments: `0/0` = alt only from
   error (k ~ beta-binomial with mean ε = 10^(−`min_baseq`/10)/3, the worst case the BQ filter
   allows, and overdispersion 0.01), `1/1` = ref only from error, `0/1` = both alleles present at
   any fraction (k ~ Binomial(n, f), f uniform on [0.05, 0.95]). PL = Phred likelihoods, GQ =
   Phred error of the best GT (flat prior, cap 99), `./.` when GQ < 20: low depth and ambiguous
   counts are no-calls, never a forced genotype. Parameters in [docs/params.md](docs/params.md).

### Indel (same entry point, `variant_type` `indel` rows)
1. **Extract** (per VCF): FILTER exactly `PASS`; multiallelic records split; keep plain-base
   alleles whose REF and ALT lengths differ (insertions, deletions, complex del-ins; SNVs, MNVs,
   symbolic alleles dropped) and whose length change is < `indel_max_length` (50); REF must match
   the FASTA; **normalized** (left-aligned and trimmed, Tan et al. 2015 / `bcftools norm` rule) so
   the same event written differently by two callers is one allele.
2. **Union**: as for SNVs, plus each allele's *footprint* — every reference placement of the event
   that gives the same haplotype (e.g. the whole homopolymer), plus one flanking base each side.
   Alleles whose footprints share ≥ 2 bases form a **locus** and are genotyped jointly; chunks
   never split a locus.
3. **Count** (scatter): per locus, the reference window and one haplotype per allele (flanks of
   read length + allele length + 10). Every read whose **unclipped** span reaches the footprint
   (so indels the aligner soft-clipped are seen) passes the same read filters as SNVs, has bases
   with BQ < 20 masked to `N`, and is aligned **semi-globally** (whole read, free haplotype ends;
   BWA-MEM scores +1/−4/−6−k) to every haplotype with Biopython's `PairwiseAligner`. Mates' scores
   are summed, so **each fragment is counted once**; it supports the best haplotype when that beats
   the runner-up by ≥ `indel_min_score_diff` (5 = one high-quality discriminating base), else it is
   a conflict (mates support different haplotypes) or uninformative (e.g. it ends inside the
   repeat). `other` = fragments supporting another allele of the locus, never `alt`.
4. **Finalize**: manifest-checked gather as for SNVs; separate indel table, VCF and QC; GT/GQ/PL as
   for SNVs with a fixed indel error rate (`--indel_gt_error_rate` 0.01: wrong-haplotype
   assignment cannot be estimated from the counts).

### SV
1. **Normalize** (per VCF) to canonical junctions `(chrom1,pos1,strand1)–(chrom2,pos2,strand2)`:
   - pass rules: BRASS `PASS` or `.`; SvABA `PASS`; GRIDSS `PASS` **and** matched-normal
     `FORMAT/VF ≤ gridss_max_normal_vf` (GRIDSS 2.2.2 output is a joint normal+tumor call set);
   - breakends parsed from ALT (VCF 4.2: position = retained base next to the junction; `+` =
     sequence kept to the left), so a mate missing from the file is recovered; `<DEL>`, `<DUP>`,
     `<INV>` expanded; single breakends and off-contig junctions dropped;
   - **insertions** kept: `<INS>` (and subtypes; length from `SVLEN`/`SVINSLEN`/`INSLEN`, sequence
     from `SVINSSEQ`), sequence-resolved REF/ALT alleles and breakends carrying inserted sequence,
     as junctions `(c, pos, +)–(c, pos + del_len + 1, −)` with `ins_len`, `ins_seq` (or
     length-only) — left-aligned against the FASTA so callers placing one insertion at different
     points of a repeat agree; `<INS>` without length or sequence dropped;
   - size: an event is kept when its net length change is ≥ `sv_min_size` (deleted bases for
     DEL, duplicated bases for DUP, net inserted bases for INS) — the exact complement of the
     indel workflow's `|len(ALT) − len(REF)| < indel_max_length` when both are 50;
   - ends ordered canonically; mate records collapse to one junction.
2. **Union**: every distinct **adjacency** (identical ends; for insertions, identical sequence or
   length) is a genotyping candidate, with all callers/samples that reported it. Adjacencies with
   the same orientations are grouped into **clusters** for summaries: linked when both ends are
   within `sv_slop`, merged only while every member stays within `sv_slop` of every other at both
   ends (bounded diameter: no transitive chaining of distant events; independent of input order;
   `spread_1/2` and `n_linked_clusters` reported). Representative: precise > imprecise,
   sequence-resolved > length-only insertions, then `sv_rep_priority`, then narrowest CI.
   Non-insertion adjacencies are written for svtyper as `<DEL>` (intrachromosomal +/−), `<DUP>`
   (−/+) or BND pairs with `MATEID` (never split across chunks); every record has `CIPOS`
   (svtyper requires it). svtyper adds 1 to the POS of reverse breakends, so `−` ends are
   written one base left of the junction; this convention is pinned by a test against svtyper's
   own parser.
3. **svtyper**: library stats once (`SVTYPER_LIBINFO`), then per chunk with
   `--genotype_on_fragments --keep_all_ref --both_sides --clip_read_support --clip_partner_match
   --clip_min_base_quality 20` (+ `-T` when a FASTA is given):
   - fragment counts count each read pair once (unbiased fragment VAF);
   - `--keep_all_ref` stops svtyper zeroing one reference evidence type, which inflates low VAFs;
   - `--both_sides` stops split reads of a *different* junction that shares a breakend from
     counting (common in a multi-sample union);
   - partner-matched clips recover junction reads BWA could not split, rejecting clips that are
     just local reference.
   `--max_reads` stays off: when exceeded svtyper emits `./.` with no counts.
4. **Insertions** (svtyper skips SVTYPE=INS; needs `--reference_fasta`): each insertion adjacency
   is genotyped by realigning reads to the reference and insertion haplotypes with the indel
   machinery (soft clips included, fragments counted once, `indel_min_score_diff`); length-only
   calls use an N block (the aligner's wildcard), so only reads crossing the insertion point
   decide. IMPRECISE insertions are reported but not genotyped.
5. **Finalize**: manifest-checked gathers (svtyper and insertion chunks); BND pairs read from `_1`
   only; one row per adjacency (`<S>_sv_adjacency_genotyping.csv.gz`, `adj_id` → `sv_id` and
   source calls) and one per cluster (`<S>_sv_genotyping.csv.gz`: the representative adjacency's
   genotype, as in v1.2, plus a summary of all adjacencies), each with `genotype_status`; the
   adjacency table also as BEDPE (`<S>_sv_adjacency_genotyping.bedpe`, for IGV and bedtools).

### Fixes relative to Mondrian
- mondrian_nf `PYSAMGENOTYPER` scatter outputs were never gathered → explicit gather + check.
- mondrianutils `SnvGenotyper` relied on pysam defaults (orphans dropped, BQ ≥ 13, overlapping
  mates double-counted, supplementary kept) → all filters explicit, fragments counted once.
- mondrianutils `SvGenotyper` matched destruct read IDs (only works on the BAM destruct ran on)
  and copied breakend-1 fields into breakend-2 → replaced by breakpoint re-genotyping (svtyper).

## isabl integration

The isabl apps `SNV_GENOTYPING` / `SV_GENOTYPING` (papaemmelab `isabl_apps`, GRCh37) write
`samplesheet.csv` from the individual's SUCCEEDED caller analyses and call
`isabl_apps.nextflow.utils.create_nextflow_script` with `nf_tool="shahcompbio/nf-wgs-genotyping
-main-script <entry>.nf"`, `-revision v1.3.0`, `-profile iris` and a per-analysis work dir.
Result patterns match the file names in [docs/outputs.md](docs/outputs.md). `SNV_GENOTYPING`
(1.2.0+) also requires a SUCCEEDED `PINDEL` analysis per pair and writes `variant_type` rows for
the Pindel `flagged_vcf` and the Strelka/Mutect `indels` results.

```bash
# in the environment that runs `isabl` (with isabl_cli): install the branch with the apps
pip install -e /data1/shahs3/users/myersm2/repos/genotyping-dev/isabl_apps   # feature/wgs-genotyping
# admin, once: register SNVGenotypingGRCh37 / SVGenotypingGRCh37 in the client's
# INSTALLED_APPLICATIONS. No settings patch is needed: nf_tool, nf_version and the GRCh37
# reference have defaults. Do not run `isabl-apps patch-settings` from this branch -- it
# re-patches every app with this branch's settings.

# one target (tumor or normal WGS experiment); without --commit isabl only stages the analysis
# (creates the record, writes .script.sh and samplesheet.csv) -- inspect them, then submit
isabl apps-grch37 snv-genotyping-1.3.0 -fi system_id <TARGET_SYSTEM_ID>
isabl apps-grch37 snv-genotyping-1.3.0 -fi system_id <TARGET_SYSTEM_ID> --commit
isabl apps-grch37 sv-genotyping-1.1.0  -fi system_id <TARGET_SYSTEM_ID> --commit

# every WGS experiment of a patient as targets
isabl apps-grch37 snv-genotyping-1.3.0 -fi sample__individual__system_id <INDIVIDUAL_SYSTEM_ID> \
    -fi technique__method WG --commit

# genotype even if some of the patient's tumors have no caller analyses
isabl apps-grch37 sv-genotyping-1.1.0 -fi system_id <TARGET_SYSTEM_ID> --allow-unpaired-tumors --commit
```

Requires the `v1.3.0` tag of `shahcompbio/nf-wgs-genotyping` to be readable by Nextflow (a
public repo, or a GitHub token in `$NXF_HOME/scm` of the user running the head job; an anonymous
request to a private repo returns `Remote resource not found`) and the image
`mattmyers/nf-wgs-genotyping:1.3.0` on Docker Hub. App versions: `SNV_GENOTYPING` 1.3.0 (GT;
1.2.0 added indels), `SV_GENOTYPING` 1.1.0 (every adjacency, insertions).

## Software

- `environment.yml` (conda; also used by the Dockerfile). The svtyper fork is pip-installed from
  GitHub at commit `27a042710ddc6c2efa224b6b6cc50b8b45231c0c`.
- `Dockerfile`: `docker build -t mattmyers/nf-wgs-genotyping:1.3.0 . && docker push mattmyers/nf-wgs-genotyping:1.3.0`
  (`params.container`; another image: `--container [registry/]<namespace>/<image>:<tag>`, no `docker:` prefix).

## Testing

```bash
conda activate nf-wgs-genotyping                      # or any env with pysam/pandas/pyyaml/biopython/pytest
pytest -q tests                                       # unit + svtyper integration + nextflow runs
SVTYPER_SRC=/path/to/svtyper pytest -q tests          # if svtyper is not installed in the env
# in an activated env without containers, pin the interpreter used inside tasks:
nextflow run snv_genotyping.nf -profile test_snv --python $(which python)   # synthetic SNV + indel truth
nextflow run sv_genotyping.nf  -profile test_sv  --python $(which python)   # NA12878 deletions
nextflow run sv_genotyping.nf  -profile test_sv  --python $(which python) \  # synthetic insertions
    --samplesheet tests/data/synthetic/sv_ins_samplesheet.csv --target_bam tests/data/synthetic/target.bam \
    --sample_id T1 --reference_fasta tests/data/synthetic/ref.fa
nextflow run snv_genotyping.nf -profile stub -stub-run
```

`tests/make_fixtures.py` regenerates the synthetic data deterministically (contig 3 holds the
designed indel loci: a unique deletion, two deletions in one homopolymer, a dinucleotide-repeat
insertion written in two representations, a complex del-ins; contig 4 the SV insertions: a unique
60-bp insertion called with its sequence and as length-only, a 50-bp repeat insertion written at
both ends of the repeat, a 49-bp one, an imprecise one, one without length, one without reads). `tests/data/svtyper`
holds MIT-licensed fixtures from the svtyper fork (see its README). No patient data is included.

## Citation

SVTyper: Chiang et al., *Nature Methods* 2015 (doi:10.1038/nmeth.3505); fork by Marc J. Williams.
