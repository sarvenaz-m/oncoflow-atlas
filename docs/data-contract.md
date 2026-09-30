# Data contract

OncoFlow Atlas uses three deliberately simple tabular inputs so that the workflow can be inspected without institution-specific infrastructure.

## Metadata

Required fields:

| Field | Type | Constraint |
|---|---|---|
| `sample_id` | string | unique; must match expression columns |
| `subtype` | string | cohort grouping label |
| `purity` | float | 0–1 |

Additional metadata columns are retained where possible and may be used for plot symbols or downstream extensions.

## Expression

A wide matrix with one `gene` column followed by one numeric column per sample. The current implementation performs `log1p` transformation for exploratory analysis; it does not perform raw-count normalization and therefore should not be interpreted as a replacement for an RNA-seq quantification pipeline.

## Somatic variants

Required fields: `sample_id`, `gene`, `variant_class`, `vaf`, and `impact`. The `--min-vaf` parameter controls the exploratory filtering threshold.

## Validation behavior

The analysis fails early for duplicated identifiers, sample mismatches, invalid purity/VAF ranges, or missing required columns. Review-level cohort QC is written separately to `qc_summary.csv`.
