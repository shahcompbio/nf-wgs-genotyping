process SV_INS_COUNT {
    tag "${chunk.baseName}"
    label 'process_low'

    input:
    path chunk
    tuple path(bam), path(bai)
    tuple path(fasta), path(fai)

    output:
    path "ins_counts_${chunk.baseName}.tsv", emit: counts
    path "versions.yml"                    , emit: versions

    script:
    // without a FASTA every insertion is reported as skipped_no_reference
    def ref  = fasta.name.startsWith('NO_FASTA') ? '' : "--reference ${fasta}"
    def dups = params.count_duplicates ? '--count-duplicates' : ''
    """
    ${params.python} \$(command -v sv_ins_count.py) --bam ${bam} ${ref} --sites ${chunk} \\
        --out ins_counts_${chunk.baseName}.tsv --min-mapq ${params.min_mapq} --min-baseq ${params.min_baseq} \\
        --min-score-diff ${params.indel_min_score_diff} ${dups}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(${params.python} -c 'import pysam; print(pysam.__version__)')
        biopython: \$(${params.python} -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """

    stub:
    """
    printf 'adj_id\\tref_fragments\\talt_fragments\\tuninformative_fragments\\tconflict_fragments\\tins_status\\n' > ins_counts_${chunk.baseName}.tsv
    printf '"${task.process}":\\n    stub: true\\n' > versions.yml
    """
}
