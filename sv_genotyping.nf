#!/usr/bin/env nextflow
/*
 * nf-wgs-genotyping — SV entry point.
 *   nextflow run <org>/nf-wgs-genotyping -main-script sv_genotyping.nf \
 *     --samplesheet ss.csv --target_bam T.bam --sample_id T --reference_fasta ref.fa --outdir out
 * Union of PASS SVs (BRASS/GRIDSS/SvABA/generic VCFs in the samplesheet), clustered across callers
 * and samples, re-genotyped in the target BAM with svtyper (marcjwilliams1 fork). See docs/.
 */
nextflow.enable.dsl = 2

include { requireParams; targetBam; referenceFasta; samplesheetRows; publishVersions } from './subworkflows/local/utils'
include { SV_GENOTYPING } from './workflows/sv_genotyping'

workflow {
    requireParams(['samplesheet', 'target_bam', 'sample_id'])
    def rows = samplesheetRows(params.samplesheet, ['brass', 'gridss', 'svaba', 'generic'])
    SV_GENOTYPING(Channel.fromList(rows), Channel.value(targetBam()), Channel.value(referenceFasta()),
                  params.sample_id)
    publishVersions(SV_GENOTYPING.out.versions)
}
