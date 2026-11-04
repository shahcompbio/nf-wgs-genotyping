process SV_FINALIZE {
    tag "${sample_id}"
    label 'process_low'

    input:
    path union
    path genotyped
    path manifest
    path stats
    val sample_id
    tuple path(fasta), path(fai)

    output:
    path "${sample_id}_sv_genotyping.csv.gz"     , emit: table
    path "${sample_id}_sv_genotyping.csv.gz.yaml", emit: yaml
    path "${sample_id}_genotyping_qc.json"       , emit: qc
    path "${sample_id}_svtyper.vcf.gz{,.tbi}"    , emit: vcf
    path "versions.yml"                          , emit: versions

    script:
    def ref     = fasta.name.startsWith('NO_FASTA') ? '' : "--reference ${fasta}"
    def contigs = params.contigs ? "--contigs '${params.contigs}'" : ''
    """
    ${params.python} \$(command -v sv_finalize.py) --union ${union} --genotyped ${genotyped} --manifest ${manifest} \\
        --stats ${stats} --sample-id '${sample_id}' ${ref} ${contigs} --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nf-wgs-genotyping: ${workflow.manifest.version}
        pandas: \$(${params.python} -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${sample_id}_sv_genotyping.csv.gz
    touch ${sample_id}_sv_genotyping.csv.gz.yaml
    echo '{}' > ${sample_id}_genotyping_qc.json
    echo | gzip > ${sample_id}_svtyper.vcf.gz
    touch ${sample_id}_svtyper.vcf.gz.tbi
    printf '"${task.process}":\\n    stub: true\\n' > versions.yml
    """
}
