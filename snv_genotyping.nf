#!/usr/bin/env nextflow
/*
 * nf-wgs-genotyping — SNV + indel entry point.
 *   nextflow run <org>/nf-wgs-genotyping -main-script snv_genotyping.nf \
 *     --samplesheet ss.csv --target_bam T.bam --sample_id T --reference_fasta ref.fa --outdir out
 * Union of PASS SNVs (CaVEMan/Mutect/Strelka/generic VCFs; samplesheet variant_type `snv`, the
 * default) genotyped with fragment-level pysam counting, and union of PASS indels
 * (Pindel/Strelka/Mutect/generic VCFs; variant_type `indel`) genotyped by realigning fragments
 * to ref/alt haplotypes. See docs/params.md and docs/outputs.md.
 */
nextflow.enable.dsl = 2

include { requireParams; targetBam; referenceFasta; samplesheetRows; publishVersions } from './subworkflows/local/utils'
include { SNV_GENOTYPING } from './workflows/snv_genotyping'

workflow {
    requireParams(['samplesheet', 'target_bam', 'sample_id'])
    def rows = samplesheetRows(params.samplesheet,
                               [snv: ['caveman', 'mutect', 'strelka', 'generic'],
                                indel: ['pindel', 'mutect', 'strelka', 'generic']], 'snv')
    if (!params.reference_fasta && rows.any { it[0].variant_type == 'indel' }) {
        error "Indel genotyping needs --reference_fasta (normalization and ref/alt haplotypes)"
    }
    SNV_GENOTYPING(Channel.fromList(rows), Channel.value(targetBam()), Channel.value(referenceFasta()),
                   params.sample_id)
    publishVersions(SNV_GENOTYPING.out.versions)
}
