process SNV_COUNT {
    tag "${chunk.baseName}"
    label 'process_low'

    input:
    path chunk
    tuple path(bam), path(bai)

    output:
    path "counts_${chunk.baseName}.tsv", emit: counts
    path "versions.yml"                , emit: versions

    script:
    def dups = params.count_duplicates ? '--count-duplicates' : ''
    """
    snv_count.py --bam ${bam} --sites ${chunk} --out counts_${chunk.baseName}.tsv \\
        --min-mapq ${params.min_mapq} --min-baseq ${params.min_baseq} ${dups}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python3 -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """

    stub:
    """
    printf 'chromosome\\tposition\\tref\\talt\\tref_count\\talt_count\\tother_count\\tconflict_count\\n' > counts_${chunk.baseName}.tsv
    printf '"${task.process}":\\n    stub: true\\n' > versions.yml
    """
}
