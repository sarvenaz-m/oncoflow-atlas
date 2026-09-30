process MULTIOMIC_INTEGRATION {
    tag "${params.cohort_name}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path expression
    path variants
    path metadata

    output:
    path "analysis/*", emit: results

    cpus 2
    memory '4 GB'
    time '30m'

    script:
    """
    mkdir -p analysis
    python ${projectDir}/bin/analyse_multiomics.py \\
      --expression ${expression} \\
      --variants ${variants} \\
      --metadata ${metadata} \\
      --outdir analysis \\
      --cohort-name '${params.cohort_name}' \\
      --min-vaf ${params.min_vaf} \\
      --top-n ${params.top_n}
    """
}
