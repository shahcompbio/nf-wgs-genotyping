#!/usr/bin/env nextflow
/*
 * nf-wgs-genotyping — SNV entry point.
 *   nextflow run <org>/nf-wgs-genotyping -main-script snv_genotyping.nf \
 *     --samplesheet ss.csv --target_bam T.bam --sample_id T --reference_fasta ref.fa --outdir out
 * Union of PASS SNVs (CaVEMan/Mutect/Strelka/generic VCFs in the samplesheet) genotyped in the
 * target BAM with fragment-level pysam counting. See docs/params.md and docs/outputs.md.
 */
nextflow.enable.dsl = 2

include { requireParams; targetBam; referenceFasta; samplesheetRows; publishVersions } from './subworkflows/local/utils'
include { SNV_GENOTYPING } from './workflows/snv_genotyping'

workflow {
    requireParams(['samplesheet', 'target_bam', 'sample_id'])
    def rows = samplesheetRows(params.samplesheet, ['caveman', 'mutect', 'strelka', 'generic'])
    SNV_GENOTYPING(Channel.fromList(rows), Channel.value(targetBam()), Channel.value(referenceFasta()),
                   params.sample_id)
    publishVersions(SNV_GENOTYPING.out.versions)
}
