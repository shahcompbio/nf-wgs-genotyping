process INDEL_COUNT {
    tag "${chunk.baseName}"
    label 'process_low'

    input:
    path chunk
    tuple path(bam), path(bai)
    tuple path(fasta), path(fai)

    output:
    path "counts_${chunk.baseName}.tsv", emit: counts
    path "versions.yml"                , emit: versions

    script:
    def dups = params.count_duplicates ? '--count-duplicates' : ''
    """
    ${params.python} \$(command -v indel_count.py) --bam ${bam} --reference ${fasta} --sites ${chunk} \\
        --out counts_${chunk.baseName}.tsv --min-mapq ${params.min_mapq} --min-baseq ${params.min_baseq} \\
        --min-score-diff ${params.indel_min_score_diff} ${dups}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(${params.python} -c 'import pysam; print(pysam.__version__)')
        biopython: \$(${params.python} -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """

    stub:
    """
    printf 'chromosome\\tposition\\tref\\talt\\tref_count\\talt_count\\tother_count\\tuninformative_count\\tconflict_count\\n' > counts_${chunk.baseName}.tsv
    printf '"${task.process}":\\n    stub: true\\n' > versions.yml
    """
}
