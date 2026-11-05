include { SV_NORMALIZE    } from '../modules/local/sv_normalize/main'
include { SV_UNION        } from '../modules/local/sv_union/main'
include { SVTYPER_LIBINFO } from '../modules/local/svtyper_libinfo/main'
include { SVTYPER         } from '../modules/local/svtyper/main'
include { SV_INS_COUNT    } from '../modules/local/sv_ins_count/main'
include { SV_FINALIZE     } from '../modules/local/sv_finalize/main'

workflow SV_GENOTYPING {
    take:
    ch_vcfs     // [meta, vcf]
    ch_bam      // value: [bam, bai]
    ch_ref      // value: [fasta, fai]
    sample_id   // value

    main:
    SV_NORMALIZE(ch_vcfs, ch_ref)
    SV_UNION(SV_NORMALIZE.out.tsv.collect(), sample_id, ch_ref)
    SVTYPER_LIBINFO(ch_bam, sample_id)
    SVTYPER(SV_UNION.out.chunks.flatten(), ch_bam, SVTYPER_LIBINFO.out.lib, ch_ref)
    // insertions (svtyper skips SVTYPE=INS): realignment to ref/insertion haplotypes
    SV_INS_COUNT(SV_UNION.out.ins_chunks.flatten(), ch_bam, ch_ref)
    // explicit gather: sv_finalize.py checks both sets of chunk ids against their manifests
    SV_FINALIZE(SV_UNION.out.union, SV_UNION.out.adjacencies, SVTYPER.out.vcf.collect(), SV_UNION.out.manifest,
                SV_INS_COUNT.out.counts.collect(), SV_UNION.out.ins_manifest,
                SV_NORMALIZE.out.stats.collect(), sample_id, ch_ref)

    ch_versions = SV_NORMALIZE.out.versions.mix(SV_UNION.out.versions, SVTYPER_LIBINFO.out.versions,
                                                SVTYPER.out.versions, SV_INS_COUNT.out.versions,
                                                SV_FINALIZE.out.versions)

    emit:
    table           = SV_FINALIZE.out.table
    adjacency_table = SV_FINALIZE.out.adjacency_table
    bedpe           = SV_FINALIZE.out.bedpe
    qc              = SV_FINALIZE.out.qc
    versions        = ch_versions
}
