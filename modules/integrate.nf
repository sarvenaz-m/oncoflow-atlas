def quote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

process VALIDATE_COHORT {
    tag 'structural and numeric validation'
    input:
    path expression, stageAs: 'expression.csv'
    path variants, stageAs: 'variants.csv'
    path metadata, stageAs: 'metadata.csv'
    path cna, stageAs: 'cna.csv'
    path code, stageAs: 'src'
    output:
    path 'validated/validation.json', emit: validation
    script:
    def cna_arg = cna ? '--cna cna.csv' : ''
    """
    PYTHONPATH="\$PWD/src" python -m oncoflow_atlas validate \\
      --expression expression.csv --variants variants.csv --metadata metadata.csv \\
      ${cna_arg} --outdir validated
    """
}

process ANALYSE_COHORT {
    tag 'matched multi-omic analysis'
    input:
    path expression, stageAs: 'expression.csv'
    path variants, stageAs: 'variants.csv'
    path metadata, stageAs: 'metadata.csv'
    path cna, stageAs: 'cna.csv'
    path source, stageAs: 'source_manifest.json'
    path code, stageAs: 'src'
    path validation, stageAs: 'validation.json'
    val cohort_name
    val min_vaf
    val min_purity
    val top_n
    output:
    path 'analysed', emit: analysis
    script:
    def cna_arg = cna ? '--cna cna.csv' : ''
    """
    OPENBLAS_NUM_THREADS=1 PYTHONPATH="\$PWD/src" python -m oncoflow_atlas analyse \\
      --expression expression.csv --variants variants.csv --metadata metadata.csv \\
      ${cna_arg} --source-manifest source_manifest.json --validation validation.json \\
      --cohort-name ${quote(cohort_name)} --min-vaf ${quote(min_vaf)} \\
      --min-purity ${quote(min_purity)} --top-n ${quote(top_n)} --outdir analysed
    """
}

process RENDER_REPORT {
    tag 'offline report and output manifest'
    publishDir params.outdir, mode: 'copy', overwrite: true
    input:
    path analysis, stageAs: 'analysed'
    path code, stageAs: 'src'
    output:
    path 'report', emit: report
    script:
    """
    mkdir report
    cp analysed/* report/
    PYTHONPATH="\$PWD/src" python -m oncoflow_atlas render \\
      --analysis report/analysis.json --outdir report
    """
}
