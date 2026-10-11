# Reproduction and execution record

## Release environment

The public and synthetic Python workflows were executed on Python 3.12.14 with NumPy 2.3.5 and pandas 2.2.3. `requirements-lock.txt` pins the complete core dependency set used for those runs. Independent checks used SciPy 1.17.0 and scikit-learn 1.8.0. Node.js 24.19.0 runs cross-language checks; browser engine details are recorded in `validation/browser-checks.json`.

The input data and analytical protocol are fixed. Numeric outputs should agree within declared tolerances across supported platforms; SVD backends can differ in the last floating-point digits. Runtime/platform/timestamp metadata in `run_manifest.json` is intentionally not byte-identical across runs. `analysis.json` and the portable HTML contain no generated timestamp.

## From original source to report

```bash
python scripts/fetch_public_data.py --outdir data/public/source  # optional re-download
python scripts/prepare_public_data.py
python scripts/build_reports.py
OPENBLAS_NUM_THREADS=1 python scripts/verify_public_results.py
```

Preparation first verifies 15 source hashes and sizes. It writes matched inputs and a row-level exclusion ledger, records all prepared-file hashes, then validates the prepared cohort. The analysis checks that its prepared files match the source manifest. Each result directory contains:

- `validation.json`: structural/numeric checks and validated input hashes.
- `analysis.json`: complete portable analysis, matrices and records.
- `pca_coordinates.csv`, `pca_loadings.csv`, `sample_correlation.csv`.
- `sample_metrics.csv`, `top_genes.csv`, `variant_summary.csv`, `cna_summary.csv`, `subtype_summary.csv`, `qc_summary.csv`, `sensitivity.csv`.
- `cohort_summary.json`, `report.html`, `run_manifest.json`.

The run manifest hashes every listed output and all Python/web source files. It cannot hash itself. The source manifest and lock link each result back to the original dataset revision. The public dataset has no random sampling or model-training seed. Synthetic generator seed 19 is retained for reproducible software checks.

## Separate Python stages

```bash
oncoflow-atlas validate --expression expression.csv --variants variants.csv --metadata metadata.csv --outdir results/staged
oncoflow-atlas analyse --expression expression.csv --variants variants.csv --metadata metadata.csv --validation results/staged/validation.json --outdir results/staged
oncoflow-atlas render --analysis results/staged/analysis.json --outdir results/staged
```

Add `--cna cna.csv` to validate/analyse if supplied and `--source-manifest source_manifest.json` to analyse when available. These stages were tested against a direct run: all scientific outputs and the report match byte-for-byte. Altered input after validation is explicitly rejected.

## Nextflow status

The DSL2 graph separates `VALIDATE_COHORT`, `ANALYSE_COHORT` and `RENDER_REPORT`. Deep input caching includes the source tree; the final directory is published only by the render stage. To run with the bundled public cohort:

```bash
nextflow run main.nf -profile public,local --outdir results/nextflow-public
nextflow run main.nf -profile public,local --outdir results/nextflow-public -resume
```

Outputs are under `results/nextflow-public/report/`; trace, timeline and workflow execution reports are beside that directory. Nextflow requires a compatible Java runtime; the attempted release run used Nextflow 26.04.6 and Java 17.0.20. The launcher parsed the workflow and submitted validation, but child-process execution emitted `fatal library error, lookup self` and did not finish. A bounded retry failed the same way. No complete Nextflow run, cache hit count, memory benchmark or timing claim is presented as successful. CI contains an actual public-cohort run followed by a strict three-task cache-resume assertion.

Docker and SLURM are configuration targets, not verified deployments. Build the local image before the Docker profile:

```bash
docker build -t oncoflow-atlas:2.0.0 .
nextflow run main.nf -profile public,docker
```

The Docker image pins Python dependencies but uses the `python:3.12-slim` base tag, not an immutable image digest. The SLURM profile uses `SLURM_QUEUE` or `compute`; adapt the queue/resource settings for the target cluster. GitHub Actions has not been triggered by this prepared package; the user-controlled push will execute it.

## Rebuild the workbench

Edit `src/oncoflow_atlas/web/`, then run `python scripts/build_reports.py`. The renderer embeds both data and scripts with script hashes in the content security policy. Opening `docs/index.html` locally exercises the same source as the standalone Python report. `docs/demo.html` redirects to the combined workbench so the former report link stays useful.

Python/browser parity and browser interaction tests are distinct checks. The first recursively compares view values in Node; the second operates the rendered interface, exports files, imports invalid/hostile report labels and checks desktop/mobile layouts without data network requests.
