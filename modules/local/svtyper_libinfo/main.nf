process SVTYPER_LIBINFO {
    tag "${sample_id}"
    label 'process_single'

    input:
    tuple path(bam), path(bai)
    val sample_id

    output:
    path "${sample_id}.lib_info.json", emit: lib
    path "versions.yml"              , emit: versions

    script:
    """
    svtyper_run.py --svtyper-cmd '${params.svtyper_cmd}' --bam ${bam} \\
        --lib-info ${sample_id}.lib_info.json --lib-only

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtyper: \$(python3 -c 'import svtyper.version as v; print(v.__version__)' 2>/dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    echo '{}' > ${sample_id}.lib_info.json
    printf '"${task.process}":\\n    stub: true\\n' > versions.yml
    """
}
