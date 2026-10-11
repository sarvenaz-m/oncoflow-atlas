# Data license and attribution

The root MIT license applies to OncoFlow Atlas software and authored synthetic examples. It does not relicense the public TCGA/cBioPortal data.

The database in `data/public/`, its adapted tables and data-bearing exports are distributed under the [Open Data Commons Open Database License 1.0](https://opendatacommons.org/licenses/odbl/1-0/), following [cBioPortal Datahub's licensing notice](https://github.com/cBioPortal/datahub#license). Preserve attribution, keep shared adapted databases open, and apply the same database license to the adapted database. The upstream study's separate TCGA notice is retained byte-for-byte at [`data/public/source/LICENSE`](data/public/source/LICENSE).

**Attribution:** The Cancer Genome Atlas (TCGA), PanCancer Atlas, *Cell* (2018); cholangiocarcinoma study `chol_tcga_pan_can_atlas_2018`; distributed by cBioPortal Datahub. Original [PanCancer publications and data](https://gdc.cancer.gov/about-data/publications/pancanatlas), [study metadata](data/public/source/meta_study.txt), and [pinned source folder](https://github.com/cBioPortal/datahub/tree/dca75cb3f32b82d54a6f78bf0a6323e5b975aca1/public/chol_tcga_pan_can_atlas_2018).

**Adaptations:** exact specimen intersection; exclusion of missing/ambiguous gene-symbol or Entrez rows; selected specimen annotations; VAF calculated from tumor reference/alternate read counts; conversion to CSV; derived descriptive results. The original selected source files, complete adapted tables, transformation code and exclusion ledger are provided together. Public report JSON/HTML and public-data screenshots contain or visualize these attributed data.

The upstream study notice retains its original caution about preliminary TCGA data. The software documentation records the actual scope of this processed-data analysis.
