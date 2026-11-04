process SV_NORMALIZE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(vcf)
    tuple path(fasta), path(fai)

    output:
    path "${meta.id}.svs.tsv"       , emit: tsv
    path "${meta.id}.svs.stats.json", emit: stats
    path "versions.yml"             , emit: versions

    script:
    def ref     = fasta.name.startsWith('NO_FASTA') ? '' : "--reference ${fasta}"
    def contigs = params.contigs ? "--contigs '${params.contigs}'" : ''
    def brass   = params.brass_missing_filter_is_pass ? '' : '--no-brass-missing-filter-is-pass'
    """
    sv_normalize.py \\
        --vcf ${vcf} --caller ${meta.caller} \\
        --tumor-id '${meta.tumor_id}' --normal-id '${meta.normal_id}' \\
        ${contigs} ${ref} ${brass} \\
        --min-size ${params.sv_min_size} --imprecise-ci ${params.sv_imprecise_ci} \\
        --gridss-max-normal-vf ${params.gridss_max_normal_vf} \\
        --out-tsv ${meta.id}.svs.tsv --out-stats ${meta.id}.svs.stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
    END_VERSIONS
    """

    stub:
    """
    printf 'chrom1\\n' > ${meta.id}.svs.tsv
    echo '{}' > ${meta.id}.svs.stats.json
    printf '"${task.process}":\\n    stub: true\\n' > versions.yml
    """
}
