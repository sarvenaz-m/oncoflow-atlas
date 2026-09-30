nextflow.enable.dsl=2
include { MULTIOMIC_INTEGRATION } from './modules/integrate'

params.expression  = "${projectDir}/data/demo/expression.csv"
params.variants    = "${projectDir}/data/demo/variants.csv"
params.metadata    = "${projectDir}/data/demo/metadata.csv"
params.outdir      = "${projectDir}/results/nextflow"
params.cohort_name = 'demo-cohort'
params.min_vaf     = 0.05
params.top_n       = 12

workflow {
    expression_ch = Channel.fromPath(params.expression, checkIfExists: true)
    variants_ch   = Channel.fromPath(params.variants, checkIfExists: true)
    metadata_ch   = Channel.fromPath(params.metadata, checkIfExists: true)

    MULTIOMIC_INTEGRATION(expression_ch, variants_ch, metadata_ch)
}
