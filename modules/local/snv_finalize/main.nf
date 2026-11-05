process SNV_FINALIZE {
    tag "${sample_id}"
    label 'process_low'

    input:
    path union
    path union_vcf      // [vcf.gz, vcf.gz.tbi] from SNV_UNION; only the header is reused
    path counts
    path manifest
    path stats
    val sample_id

    output:
    path "${sample_id}_snv_genotyping.csv.gz"     , emit: table
    path "${sample_id}_snv_genotyping.csv.gz.yaml", emit: yaml
    path "${sample_id}_genotyping_qc.json"        , emit: qc
    path "${sample_id}_snv_genotyping.vcf.gz"     , emit: vcf
    path "${sample_id}_snv_genotyping.vcf.gz.tbi" , emit: tbi
    path "versions.yml"                           , emit: versions

    script:
    def vcf = [union_vcf].flatten().find { it.name.endsWith('.vcf.gz') }
    """
    ${params.python} \$(command -v snv_finalize.py) --union ${union} --union-vcf ${vcf} --counts ${counts} \\
        --manifest ${manifest} --stats ${stats} --sample-id '${sample_id}' --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nf-wgs-genotyping: ${workflow.manifest.version}
        pandas: \$(${params.python} -c 'import pandas; print(pandas.__version__)')
        pysam: \$(${params.python} -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${sample_id}_snv_genotyping.csv.gz
    touch ${sample_id}_snv_genotyping.csv.gz.yaml
    echo '{}' > ${sample_id}_genotyping_qc.json
    echo | gzip > ${sample_id}_snv_genotyping.vcf.gz
    touch ${sample_id}_snv_genotyping.vcf.gz.tbi
    printf '"${task.process}":\\n    stub: true\\n' > versions.yml
    """
}
