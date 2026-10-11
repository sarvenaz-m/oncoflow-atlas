nextflow.enable.dsl=2
include { VALIDATE_COHORT; ANALYSE_COHORT; RENDER_REPORT } from './modules/integrate'

workflow {
    expression = file(params.expression, checkIfExists: true)
    variants = file(params.variants, checkIfExists: true)
    metadata = file(params.metadata, checkIfExists: true)
    cna = params.cna ? [file(params.cna, checkIfExists: true)] : []
    source = file(params.source_manifest, checkIfExists: true)
    code = file("${projectDir}/src", checkIfExists: true)
    VALIDATE_COHORT(expression, variants, metadata, cna, code)
    ANALYSE_COHORT(expression, variants, metadata, cna, source, code,
                  VALIDATE_COHORT.out.validation, params.cohort_name,
                  params.min_vaf, params.min_purity, params.top_n)
    RENDER_REPORT(ANALYSE_COHORT.out.analysis, code)
}
