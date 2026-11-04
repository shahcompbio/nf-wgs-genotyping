process SV_UNION {
    tag "${sample_id}"
    label 'process_low'

    input:
    path tsvs
    val sample_id
    tuple path(fasta), path(fai)

    output:
    path "${sample_id}_union_svs.bedpe.gz"        , emit: union
    path "${sample_id}_union_svs.members.tsv.gz"  , emit: members
    path "${sample_id}_union_svs.vcf.gz{,.tbi}"   , emit: vcf
    path "chunks/chunk_*.vcf"                     , emit: chunks
    path "chunks.manifest"                        , emit: manifest
    path "versions.yml"                           , emit: versions

    script:
    def ref     = fasta.name.startsWith('NO_FASTA') ? '' : "--reference ${fasta}"
    def contigs = params.contigs ? "--contigs '${params.contigs}'" : ''
    """
    sv_union.py --inputs ${tsvs} --sample-id '${sample_id}' ${ref} ${contigs} \\
        --slop ${params.sv_slop} --rep-priority '${params.sv_rep_priority}' \\
        --chunk-size ${params.sv_chunk_size} --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p chunks
    printf '##fileformat=VCFv4.2\\n' | tee chunks/chunk_00000.vcf > chunks/chunk_00001.vcf
    printf 'chunk_id\\tn_items\\nchunk_00000\\t0\\nchunk_00001\\t0\\n' > chunks.manifest
    echo | gzip > ${sample_id}_union_svs.bedpe.gz
    echo | gzip > ${sample_id}_union_svs.members.tsv.gz
    echo | gzip > ${sample_id}_union_svs.vcf.gz
    touch ${sample_id}_union_svs.vcf.gz.tbi
    printf '"${task.process}":\\n    stub: true\\n' > versions.yml
    """
}
