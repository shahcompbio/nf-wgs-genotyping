include { INDEL_EXTRACT  } from '../modules/local/indel_extract/main'
include { INDEL_UNION    } from '../modules/local/indel_union/main'
include { INDEL_COUNT    } from '../modules/local/indel_count/main'
include { INDEL_FINALIZE } from '../modules/local/indel_finalize/main'

workflow INDEL_GENOTYPING {
    take:
    ch_vcfs     // [meta, vcf] with meta.variant_type == 'indel'
    ch_bam      // value: [bam, bai]
    ch_ref      // value: [fasta, fai] (required: checked by the entry script)
    sample_id   // value

    main:
    INDEL_EXTRACT(ch_vcfs, ch_ref)
    INDEL_UNION(INDEL_EXTRACT.out.tsv.collect(), sample_id, ch_ref)
    INDEL_COUNT(INDEL_UNION.out.chunks.flatten(), ch_bam, ch_ref)
    // explicit gather: indel_finalize.py checks the gathered chunk ids against the manifest
    INDEL_FINALIZE(INDEL_UNION.out.union, INDEL_UNION.out.vcf, INDEL_COUNT.out.counts.collect(),
                   INDEL_UNION.out.manifest, INDEL_EXTRACT.out.stats.collect(), sample_id)

    ch_versions = INDEL_EXTRACT.out.versions.mix(INDEL_UNION.out.versions, INDEL_COUNT.out.versions,
                                                 INDEL_FINALIZE.out.versions)

    emit:
    table    = INDEL_FINALIZE.out.table
    qc       = INDEL_FINALIZE.out.qc
    vcf      = INDEL_FINALIZE.out.vcf
    versions = ch_versions
}
