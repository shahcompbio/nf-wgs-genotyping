process INDEL_UNION {
    tag "${sample_id}"
    label 'process_low'

    input:
    path tsvs
    val sample_id
    tuple path(fasta), path(fai)

    output:
    path "${sample_id}_union_indels.tsv.gz"        , emit: union
    path "${sample_id}_union_indels.vcf.gz{,.tbi}" , emit: vcf
    path "chunks/chunk_*.tsv"                      , emit: chunks
    path "chunks.manifest"                         , emit: manifest
    path "versions.yml"                            , emit: versions

    script:
    """
    ${params.python} \$(command -v indel_union.py) --inputs ${tsvs} --sample-id '${sample_id}' \\
        --reference ${fasta} --chunk-size ${params.indel_chunk_size} --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(${params.python} -c 'import platform; print(platform.python_version())')
        pandas: \$(${params.python} -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p chunks
    printf 'chromosome\\tposition\\tref\\talt\\tlocus_id\\n' | tee chunks/chunk_00000.tsv > chunks/chunk_00001.tsv
    printf 'chunk_id\\tn_items\\nchunk_00000\\t0\\nchunk_00001\\t0\\n' > chunks.manifest
    echo | gzip > ${sample_id}_union_indels.tsv.gz
    echo | gzip > ${sample_id}_union_indels.vcf.gz
    touch ${sample_id}_union_indels.vcf.gz.tbi
    printf '"${task.process}":\\n    stub: true\\n' > versions.yml
    """
}
