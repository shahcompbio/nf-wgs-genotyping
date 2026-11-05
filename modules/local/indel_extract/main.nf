process INDEL_EXTRACT {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(vcf)
    tuple path(fasta), path(fai)

    output:
    path "${meta.id}.indels.tsv"       , emit: tsv
    path "${meta.id}.indels.stats.json", emit: stats
    path "versions.yml"                , emit: versions

    script:
    def contigs = params.contigs ? "--contigs '${params.contigs}'" : ''
    def missing = meta.caller == 'generic' ? '--missing-filter-is-pass' : ''
    """
    ${params.python} \$(command -v indel_extract.py) \\
        --vcf ${vcf} --caller ${meta.caller} \\
        --tumor-id '${meta.tumor_id}' --normal-id '${meta.normal_id}' \\
        --reference ${fasta} ${contigs} ${missing} \\
        --max-ref-mismatch-frac ${params.max_ref_mismatch_frac} --max-length ${params.indel_max_length} \\
        --out-tsv ${meta.id}.indels.tsv --out-stats ${meta.id}.indels.stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(${params.python} -c 'import platform; print(platform.python_version())')
        pysam: \$(${params.python} -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """

    stub:
    """
    printf 'chrom\\tpos\\tref\\talt\\ttumor_id\\tnormal_id\\tcaller\\n' > ${meta.id}.indels.tsv
    echo '{}' > ${meta.id}.indels.stats.json
    printf '"${task.process}":\\n    stub: true\\n' > versions.yml
    """
}
