# Executed validation evidence — v2.0.0

Recorded for this prepared release. These are local measurements, not a claim that GitHub Actions has already run.

| Check | Outcome | Evidence |
|---|---|---|
| Python input, numerical-boundary, staged-run and regression tests | 27 passed | `pytest.txt` |
| Ruff static checks | Passed | `ruff.txt` |
| Fresh pinned source download | 15 files verified | `source-download-check.json` |
| Prepared CSV reproduction | All five files reproduced byte-for-byte | `independent-checks.json` |
| NumPy PCA vs scikit-learn full PCA | Maximum coordinate error 3.56 × 10⁻¹³ or less | `independent-checks.json` |
| Heatmap standardization vs SciPy | Maximum absolute error 5.56 × 10⁻¹⁶ or less | `independent-checks.json` |
| Sample correlations vs SciPy | Maximum absolute error 4.00 × 10⁻¹⁵ or less | `independent-checks.json` |
| Report VAF vs original tumor read-count ratio | Zero observed difference | `independent-checks.json` |
| Python/JavaScript view parity | 27 scenarios; 176,440 numeric comparisons | `independent-checks.json` |
| Offline browser interactions | 37 checks passed; desktop and 390/360 px mobile | `browser-checks.json` |
| Output/current-source integrity | Passed | `artifact-integrity.json` |
| Python distribution | Wheel build and installed-package smoke run | `package-check.json` |
| Nextflow execution | Blocked at local child-process launch | `nextflow-status.json` |
| Docker / SLURM | Not executed | `../docs/reproducibility.md` |

The browser check includes seven workspaces, both data modes, exact VAF 0.30, cohort and singleton/empty filters, a four-specimen CSV export, invalid imports, escaped imported HTML/unsafe URLs, theme persistence and overflow checks. It observed no console errors or external data requests in offline mode. Browser and Playwright versions are recorded in the JSON evidence.

The independent reference check is not a separate biological validation cohort. It verifies numerical results, source transformations and software consistency. See the [data card](../docs/data-card.md) for the scientific scope.

A floating-point parser defect was detected by the raw-count threshold check: 24 records exactly at VAF 0.30 were incorrectly excluded by the pandas numeric conversion path. The corrected parser preserves decimal round trips, and the final 0.30 count is 2,079. A dedicated regression test covers that boundary.

Reproduce these checks with the README commands. `scripts/verify_artifacts.py` detects stale report outputs, altered source tables, changed renderer code and mismatched browser evidence. Runtime metadata can differ between environments; compare scientific values using the declared tolerances.
