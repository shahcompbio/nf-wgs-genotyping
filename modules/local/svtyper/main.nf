process SVTYPER {
    tag "${chunk.baseName}"
    label 'process_low'

    input:
    path chunk
    tuple path(bam), path(bai)
    path lib
    tuple path(fasta), path(fai)

    output:
    path "${chunk.baseName}.gt.vcf", emit: vcf
    path "versions.yml"            , emit: versions

    script:
    def ref       = fasta.name.startsWith('NO_FASTA') ? '' : "--reference ${fasta}"
    def max_reads = params.svtyper_max_reads ? " --max_reads ${params.svtyper_max_reads}" : ''
    """
    svtyper_run.py --svtyper-cmd '${params.svtyper_cmd}' --vcf ${chunk} --bam ${bam} \\
        --lib-info ${lib} --out ${chunk.baseName}.gt.vcf ${ref} \\
        --extra-args='${params.svtyper_args}${max_reads}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtyper: \$(python3 -c 'import svtyper.version as v; print(v.__version__)' 2>/dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    touch ${chunk.baseName}.gt.vcf
    printf '"${task.process}":\\n    stub: true\\n' > versions.yml
    """
}
