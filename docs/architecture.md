# Architecture

![Workflow architecture](../assets/architecture.svg)

OncoFlow Atlas separates **orchestration**, **analysis**, **configuration**, **validation** and **presentation** so each layer can be replaced independently.

1. **Data contract** — metadata, expression and variant tables are validated before computation.
2. **Nextflow DSL2 orchestration** — executor selection and resource policy remain outside the analysis code.
3. **Python analysis layer** — performs exploratory expression structure, variant summaries, sample integration and QC.
4. **Provenance layer** — records runtime versions, input checksums, parameters and output checksums.
5. **Presentation layer** — generates both machine-readable tables and a standalone interactive HTML report.

## Execution targets

The same workflow can run with local, SLURM or AWS Batch executors. A Singularity profile is included for HPC environments where Docker is not appropriate.

## Extension points

The tabular input contract intentionally sits downstream of assay-specific preprocessing. FASTQ alignment/quantification, VCF annotation, single-cell ingestion or spatial preprocessing can therefore be added as upstream modules without rewriting the integration/reporting layer.
