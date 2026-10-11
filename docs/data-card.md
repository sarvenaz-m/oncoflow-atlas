# Public data card — TCGA cholangiocarcinoma

## Identity and scope

- Study: `chol_tcga_pan_can_atlas_2018`, TCGA PanCancer Atlas.
- Distributor: [cBioPortal Datahub](https://github.com/cBioPortal/datahub).
- Pinned revision: `dca75cb3f32b82d54a6f78bf0a6323e5b975aca1`.
- [Study overview](https://www.cbioportal.org/study/summary?id=chol_tcga_pan_can_atlas_2018).
- [Exact upstream folder](https://github.com/cBioPortal/datahub/tree/dca75cb3f32b82d54a6f78bf0a6323e5b975aca1/public/chol_tcga_pan_can_atlas_2018).
- 36 tumor specimens, each belonging to a different patient. No normal comparator.

The sample annotations contain 30 intrahepatic cholangiocarcinomas, 4 perihilar cholangiocarcinomas and 2 records labeled only “Cholangiocarcinoma”. The last group remains as supplied; it is not reassigned. These are source histology labels, not predictions from a trained model. Only the necessary public specimen-level clinical fields are bundled; the clinical patient table is not required or included.

## Assays and transformations

| Assay | Original file | Source rows | Retained rows | Interpretation |
|---|---|---:|---:|---|
| RNA | `data_mrna_seq_v2_rsem.txt` | 20,531 | 20,466 | Upstream batch-normalized RSEM abundance |
| CNA | `data_cna.txt` | 25,128 | 24,805 | Discrete GISTIC −2, −1, 0, 1, 2 |
| Mutation records | `data_mutations.txt` | 3,615 | 3,615 before VAF filtering | Source curated records, including synonymous and noncoding classes |

RNA is linear-scale normalized abundance, not raw counts, TPM or a newly generated alignment result. Analysis uses `log2(value + 1)` without repeating upstream batch normalization. See [cBioPortal's RNA-seq description](https://docs.cbioportal.org/user-guide/faq/).

For RNA and CNA, missing gene symbols, invalid Entrez identifiers and **all members of ambiguous symbol/Entrez groups** are excluded. No arbitrary first duplicate is kept and rows are not merged by a heuristic. All 65 RNA and 323 CNA exclusions are listed in `data/public/prepared/excluded_features.csv`. Source-data row numbers count records after the header.

Specimens are intersected by exact identifiers across clinical data, RNA matrix/case list, mutation coverage/case list and CNA matrix/case list. The WXS panel coverage distinguishes a measured specimen with no mutation record from an unmeasured specimen. All 36 are jointly profiled; none is excluded at this step. Missing/negative/nonfinite expression and invalid CNA states are errors, not imputed zeroes.

VAF is `t_alt_count / (t_ref_count + t_alt_count)`. Invalid read counts are rejected; unavailable or zero-depth VAF stays null. The selected source has no unknown VAF. Consequence and IMPACT fields retain their upstream annotations; no effect predictor is fitted here. Source records are not independently reclassified as somatic variants.

Purity and usable technical batch annotations are absent from the selected fields. Purity is null for all 36 specimens and batch is “Not reported”. Tissue collection sites are not substituted for technical batch. The 0.40 purity threshold is only a configured review flag for datasets with reported purity.

## What this experiment establishes

The package establishes that these processed public tables can be matched, validated, analyzed and reproduced with explicit numerical rules. It does not establish an independent biological discovery, predictive accuracy, treatment response or clinical utility. The cohort is small, anatomically imbalanced and has no independent validation cohort or normal controls. No differential-expression significance, survival association or diagnostic classifier is reported.

GISTIC counts are counts of affected profiled genes; one chromosomal event can affect many genes. Mutation counts include all source classes and are not normalized by callable megabases. Source read-level QC, alignment, variant calling and GISTIC estimation are not rerun.

## Provenance and access

The 15 pinned source files are inventoried in `config/public-source-lock.json`; SHA-256 and sizes are verified before preparation. `source_manifest.json` records all retained samples, exclusions, source hashes and prepared-file hashes. `scripts/fetch_public_data.py` can re-download the original bytes, including custom Git LFS storage. Network access is needed only for that optional fetch and package installation.

The complete adapted dataset is available in this repository under the applicable open-data terms. [Data license and attribution](../DATA_LICENSE.md). There is no account or paid API requirement for the bundled public source.
