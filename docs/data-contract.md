# Input and report contract

Use UTF-8 comma-separated tables with unique, nonempty headers and a fixed number of fields per row. The supported input range is 3–500 specimens and 2–60,000 RNA features. Each input CSV is limited to 150 MiB; browser report import is limited to 64 MiB and 20 million cells per assay. Smaller subsets are recommended for interactive exploration.

| File | Required fields | Rules |
|---|---|---|
| `expression.csv` | `gene`, one column per specimen | Unique gene symbols; exactly the metadata specimen set; finite nonnegative linear-scale values; no all-zero specimen |
| `metadata.csv` | `sample_id`, `subtype` | Unique specimen IDs; nonempty source group label |
| Metadata optional fields | `patient_id`, `purity`, `batch`, other annotations | Patient IDs must be unique when supplied; otherwise specimen IDs stand in as identity, not a verified patient linkage; purity null or 0–1; missing batch remains “Not reported” |
| `variants.csv` | `sample_id`, `gene`, `variant_class`, `vaf`, `impact` | Known specimen IDs; VAF null or 0–1; duplicate records rejected |
| Variant optional fields | `variant_id`, positions, alleles, counts | Unique variant IDs if supplied; otherwise row IDs assigned |
| `cna.csv` (optional) | `gene`, one column per specimen | Exact specimen set; unique genes; discrete values −2, −1, 0, 1, 2 |

Specimen columns are aligned to metadata order. Gene rows are sorted deterministically. Leading/trailing identifier whitespace is rejected. Empty strings, `NA`, `NaN`, `nan` and `null` are accepted only for explicitly nullable numeric fields. Missing expression and nonnumeric strings cause an error; no silent imputation is performed.

One specimen per patient is the declared analysis unit. For repeated samples, select specimens under a documented rule before running the workflow; do not disguise repeated patients using fabricated IDs. Source group labels are not inferred from expression.

## Outputs

`analysis.json` is schema 2.0 and contains source information, matrices, records, full-cohort PCA/correlation, validation and a default view. `run_manifest.json` records input, software and output hashes plus runtime versions. All gene loadings, specimen coordinates, summaries and threshold sensitivity are also exported as CSV.

Browser filtering recomputes the view from the report matrices and records rather than trusting a stored summary. Import validates shape, IDs, bounds, dimensions and PCA coverage; it does not authenticate a third-party report's asserted origin or independently refit its PCA. Embedded strings are escaped, unsafe external URL schemes are rejected and the generated report has a content security policy with no data connections.

Current-view JSON carries applied parameters and source fingerprint. CSV contains the selected records and quotes spreadsheet formula prefixes in text fields. Full source matrices and unmodified identifiers remain in the analysis JSON; spreadsheet quoting is a presentation safeguard rather than a mutation of the scientific input.
