/*
 * Input validation and channel helpers shared by both entry scripts.
 */

def requireParams(List names) {
    def missing = names.findAll { !params[it] }
    if (missing) {
        error "Missing required parameter(s): ${missing.collect { '--' + it }.join(', ')}"
    }
    if (!(params.sample_id ==~ /[A-Za-z0-9._-]+/)) {
        error "--sample_id must match [A-Za-z0-9._-]+ (got '${params.sample_id}')"
    }
}

/* [bam, bai]: --target_bai, else <bam>.bai, else <bam without .bam>.bai */
def targetBam() {
    def bam = file(params.target_bam, checkIfExists: true)
    def bai = null
    if (params.target_bai) {
        bai = file(params.target_bai, checkIfExists: true)
    } else {
        def cands = [file("${bam}.bai"), file(bam.toString().replaceAll(/\.bam$/, '.bai'))]
        bai = cands.find { it.exists() }
        if (!bai) { error "No BAM index found for ${bam} (tried ${cands.join(', ')}); set --target_bai" }
    }
    return [bam, bai]
}

/* [fasta, fai] or the NO_FASTA placeholders when --reference_fasta is not set */
def referenceFasta() {
    if (!params.reference_fasta) {
        return [file("${projectDir}/assets/NO_FASTA"), file("${projectDir}/assets/NO_FASTA.fai")]
    }
    def fa = file(params.reference_fasta, checkIfExists: true)
    def fai = file("${fa}.fai")
    if (!fai.exists()) { error "Reference FASTA index not found: ${fai}" }
    return [fa, fai]
}

/* Samplesheet rows -> [meta, vcf]; relative VCF paths resolve against the samplesheet directory. */
def samplesheetRows(String path, List allowedCallers) {
    def ss = file(path, checkIfExists: true)
    def rows = ss.splitCsv(header: true)
    if (!rows) { error "Samplesheet ${ss} has no rows" }
    def seen = [] as Set
    return rows.collect { r ->
        ['tumor_id', 'normal_id', 'caller', 'vcf'].each { col ->
            if (!r[col]) { error "Samplesheet ${ss}: missing '${col}' in row ${r}" }
        }
        if (!(r.caller in allowedCallers)) {
            error "Samplesheet ${ss}: caller '${r.caller}' not supported by this workflow (allowed: ${allowedCallers.join(', ')})"
        }
        def id = "${r.tumor_id}_vs_${r.normal_id}.${r.caller}"
        if (!seen.add(id)) { error "Samplesheet ${ss}: duplicated tumor_id/normal_id/caller ${id}" }
        def vcf = r.vcf.startsWith('/') ? file(r.vcf) : ss.parent.resolve(r.vcf)
        if (!vcf.exists()) { error "Samplesheet ${ss}: VCF not found: ${vcf}" }
        [[id: id, tumor_id: r.tumor_id, normal_id: r.normal_id, caller: r.caller], vcf]
    }
}

def publishVersions(ch_versions) {
    ch_versions
        .unique { it.text }
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)
}
