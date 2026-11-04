include { SNV_EXTRACT  } from '../modules/local/snv_extract/main'
include { SNV_UNION    } from '../modules/local/snv_union/main'
include { SNV_COUNT    } from '../modules/local/snv_count/main'
include { SNV_FINALIZE } from '../modules/local/snv_finalize/main'

workflow SNV_GENOTYPING {
    take:
    ch_vcfs     // [meta, vcf]
    ch_bam      // value: [bam, bai]
    ch_ref      // value: [fasta, fai]
    sample_id   // value

    main:
    SNV_EXTRACT(ch_vcfs, ch_ref)
    SNV_UNION(SNV_EXTRACT.out.tsv.collect(), sample_id, ch_ref)
    SNV_COUNT(SNV_UNION.out.chunks.flatten(), ch_bam)
    // explicit gather: snv_finalize.py checks the gathered chunk ids against the manifest
    SNV_FINALIZE(SNV_UNION.out.union, SNV_COUNT.out.counts.collect(), SNV_UNION.out.manifest,
                 SNV_EXTRACT.out.stats.collect(), sample_id)

    ch_versions = SNV_EXTRACT.out.versions.mix(SNV_UNION.out.versions, SNV_COUNT.out.versions,
                                               SNV_FINALIZE.out.versions)

    emit:
    table    = SNV_FINALIZE.out.table
    qc       = SNV_FINALIZE.out.qc
    versions = ch_versions
}
