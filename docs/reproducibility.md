# Reproducibility model

The repository treats reproducibility as an explicit output of the workflow.

- Demo data generation is deterministic when the seed is fixed.
- Primary input files are recorded with SHA-256 hashes.
- Runtime versions for Python, NumPy, pandas and Plotly are captured in `run_manifest.json`.
- Result files are hashed after creation and recorded in the same manifest.
- Nextflow trace, timeline and execution report files are enabled for orchestrated runs.
- CI executes the synthetic cohort on two supported Python versions.
- Local, SLURM and AWS Batch profiles keep executor-specific settings outside analysis logic.

The example cohort is synthetic and contains no patient data.
