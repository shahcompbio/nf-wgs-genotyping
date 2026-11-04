process SNV_EXTRACT {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(vcf)
    tuple path(fasta), path(fai)

    output:
    path "${meta.id}.snvs.tsv"       , emit: tsv
    path "${meta.id}.snvs.stats.json", emit: stats
    path "versions.yml"              , emit: versions

    script:
    def ref     = fasta.name.startsWith('NO_FASTA') ? '' : "--reference ${fasta}"
    def contigs = params.contigs ? "--contigs '${params.contigs}'" : ''
    def missing = meta.caller == 'generic' ? '--missing-filter-is-pass' : ''
    """
    ${params.python} \$(command -v snv_extract.py) \\
        --vcf ${vcf} --caller ${meta.caller} \\
        --tumor-id '${meta.tumor_id}' --normal-id '${meta.normal_id}' \\
        ${contigs} ${ref} ${missing} \\
        --max-ref-mismatch-frac ${params.max_ref_mismatch_frac} \\
        --out-tsv ${meta.id}.snvs.tsv --out-stats ${meta.id}.snvs.stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(${params.python} -c 'import platform; print(platform.python_version())')
        pysam: \$(${params.python} -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """

    stub:
    """
    printf 'chrom\\tpos\\tref\\talt\\ttumor_id\\tnormal_id\\tcaller\\n' > ${meta.id}.snvs.tsv
    echo '{}' > ${meta.id}.snvs.stats.json
    printf '"${task.process}":\\n    stub: true\\n' > versions.yml
    """
}
