# Methodology

The included cohort is synthetic and designed only to exercise workflow behavior.

## Expression

Expression values are converted with `log1p` and used for exploratory PCA, sample correlation and a standardized heatmap of the most variable genes. This is intentionally downstream of assay-specific quantification and does not claim to replace count normalization or differential-expression methods.

## Somatic variants

Variants are filtered using a configurable VAF threshold, then summarized by sample, gene and annotated impact. A compact mutation landscape is generated from the most frequently altered genes.

## Integration

Per-sample expression signal, detected features, purity and variant burden are combined into a single table. Cohort summaries remain descriptive rather than inferential.

## Quality control

Hard schema/range failures stop execution. Cohort-level review checks are written to `qc_summary.csv` and surfaced in the HTML report.

## Provenance

Inputs and generated outputs are SHA-256 hashed. Runtime library versions and analysis parameters are written to `run_manifest.json`.
