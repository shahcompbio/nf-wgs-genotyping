#!/usr/bin/env nextflow
/*
 * nf-wgs-genotyping has two entry points and no default workflow:
 *   nextflow run <org>/nf-wgs-genotyping -main-script snv_genotyping.nf ...
 *   nextflow run <org>/nf-wgs-genotyping -main-script sv_genotyping.nf  ...
 */
nextflow.enable.dsl = 2

workflow {
    error "nf-wgs-genotyping has no default workflow. Run one of:\n" +
          "  nextflow run <repo> -main-script snv_genotyping.nf ...\n" +
          "  nextflow run <repo> -main-script sv_genotyping.nf ..."
}
