# Changelog

## 2.0.0 — 2026-10-05

- Added an executed 36-patient public TCGA experiment with matched RNA, mutation and GISTIC data, pinned source bytes, open-data attribution and a complete exclusion ledger.
- Unified Python/browser log2 transformation, sample variance, heatmap ranking, VAF boundaries, selected-sample exports and missingness rules.
- Fixed loss of exact threshold records during decimal CSV parsing; source VAF 0.30 remains included at threshold 0.30.
- Added strict numeric/identifier validation, unknown-VAF handling, optional CNA, gene-level assay joins and independent validate/analyse/render commands.
- Replaced separate report implementations with one offline renderer, seven workspaces, native SVG charts and responsive layouts.
- Added independent sklearn/SciPy comparisons, source-byte reproduction, Python/JavaScript parity, browser interaction checks and recorded evidence.
- Converted the analyzer to an installable Python package while retaining `bin/analyse_multiomics.py` compatibility.
- Added staged Nextflow code and CI execution/resume checks. Local orchestration was blocked by a child-process runtime error; Docker/SLURM have not been executed.
- Removed old generated web-data copies and v1 screenshots. Removed unverified AWS Batch/Singularity placeholders; Docker and SLURM configuration remain explicit optional profiles.

## 1.0.0

- Nextflow DSL2 orchestration with local, SLURM and AWS Batch profiles
- deterministic synthetic multi-omic cohort
- configurable VAF filtering and top-feature selection
- PCA, expression heatmap and sample-correlation views
- mutation burden, impact composition and mutation landscape views
- sample-level multi-omic integration and subtype summary
- fail-fast validation plus cohort-level QC
- SHA-256 input/output provenance and runtime capture
- interactive research report and machine-readable output contract
- automated tests and multi-version CI quality gate
