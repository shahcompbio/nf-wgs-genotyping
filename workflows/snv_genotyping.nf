include { SNV_EXTRACT      } from '../modules/local/snv_extract/main'
include { SNV_UNION        } from '../modules/local/snv_union/main'
include { SNV_COUNT        } from '../modules/local/snv_count/main'
include { SNV_FINALIZE     } from '../modules/local/snv_finalize/main'
include { INDEL_GENOTYPING } from './indel_genotyping'

workflow SNV_GENOTYPING {
    take:
    ch_vcfs     // [meta, vcf]; meta.variant_type 'snv' or 'indel'
    ch_bam      // value: [bam, bai]
    ch_ref      // value: [fasta, fai]
    sample_id   // value

    main:
    ch_by_type = ch_vcfs.branch { meta, vcf ->
        snv: meta.variant_type == 'snv'
        indel: meta.variant_type == 'indel'
    }
    // SNVs: no rows of this type -> collect() emits nothing and no SNV step runs
    SNV_EXTRACT(ch_by_type.snv, ch_ref)
    SNV_UNION(SNV_EXTRACT.out.tsv.collect(), sample_id, ch_ref)
    SNV_COUNT(SNV_UNION.out.chunks.flatten(), ch_bam)
    // explicit gather: snv_finalize.py checks the gathered chunk ids against the manifest
    SNV_FINALIZE(SNV_UNION.out.union, SNV_UNION.out.vcf, SNV_COUNT.out.counts.collect(), SNV_UNION.out.manifest,
                 SNV_EXTRACT.out.stats.collect(), sample_id)

    INDEL_GENOTYPING(ch_by_type.indel, ch_bam, ch_ref, sample_id)

    ch_versions = SNV_EXTRACT.out.versions.mix(SNV_UNION.out.versions, SNV_COUNT.out.versions,
                                               SNV_FINALIZE.out.versions, INDEL_GENOTYPING.out.versions)

    emit:
    table       = SNV_FINALIZE.out.table
    qc          = SNV_FINALIZE.out.qc
    vcf         = SNV_FINALIZE.out.vcf
    indel_table = INDEL_GENOTYPING.out.table
    indel_qc    = INDEL_GENOTYPING.out.qc
    indel_vcf   = INDEL_GENOTYPING.out.vcf
    versions    = ch_versions
}
