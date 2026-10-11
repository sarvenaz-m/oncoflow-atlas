# Analytical methods and interpretation

The machine-readable contract is [`config/analysis-protocol.json`](../config/analysis-protocol.json). Source selection is described in the [data card](data-card.md).

## Expression

Input values are finite, nonnegative, linear-scale normalized abundance. Transform every value with `log2(x + 1)`. No additional library-size normalization, batch correction or scaling across genes is applied. For the displayed heatmap, rank genes by descending variance across **selected** specimens with sample denominator `n − 1`; break ties by ascending gene symbol. Standardize each selected gene using that same mean and sample SD. Constant rows and singleton selections have displayed z-score zero. Empty selections display no heatmap.

PCA is fitted once to the full samples-by-genes transformed matrix after centering every gene, without scaling genes to unit variance. NumPy SVD provides coordinates and explained variance; each component's sign is oriented so its largest absolute loading is positive. Filtering selects existing coordinates and retains the full-cohort variance percentages. The report does not suggest that PCA has been refitted. A cohort with no between-sample variation is rejected. Pearson correlation uses all retained log-expression genes; a constant profile produces an unavailable correlation rather than zero.

## Mutation records

Match source specimen IDs exactly. Calculate VAF from tumor reference/alternate counts. Preserve unknown VAF as null and exclude it at every threshold, including zero. Keep a record when **VAF ≥ threshold**. Decimal parsing preserves round-tripped binary floats so a source ratio exactly at 0.30 is not rounded below the cutoff. All source consequence classes remain eligible; filtering on consequence would be a separate future protocol.

Gene tables count records and distinct affected specimens separately. The oncoprint can contain multiple records per cell; its color represents the highest upstream IMPACT among retained records and its tooltip lists them. No cell without a retained record is claimed to be biologically wild type.

The sensitivity series evaluates thresholds 0.00, 0.05, 0.10, 0.20, 0.30 and 0.50 on the same selected specimens. This illustrates threshold dependence; it is not an optimized predictive cutoff. In the public cohort the corresponding retained counts are 3,615, 3,592, 3,265, 2,701, 2,079 and 632.

## Copy number and matched genes

Preserve GISTIC calls: −2 deep deletion, −1 shallow deletion, 0 neutral, +1 gain and +2 amplification. Rank the overview by the number of specimens with +2 or −2, then gene symbol. Sample-level summaries count genes carrying each high-level call. Correlated neighboring genes can represent the same underlying event. The gene explorer joins by the retained gene symbol and exact specimen ID; an assay absent for that gene is shown as unavailable.

## Quality and missingness

Structural validation runs independently of calculation and rejects malformed tables, duplicate IDs, mismatched specimen sets, invalid numeric values, impossible VAF/purity values and non-discrete CNA. A prior validation artifact records input hashes; the analysis stage refuses a changed input. Unreported purity and unknown VAF remain explicit specimen-level review reasons. Low reported purity is a flag, not an automatic exclusion.

The machine-readable `qc_pass_rate_pct` summarizes only assessable **cohort checks**. It is not the percentage of biologically valid specimens. Unavailable fields and configured VAF filtering are excluded from that denominator and shown separately. All public specimens have a purity-review flag despite structurally valid tables. The user interface emphasizes availability and concrete flags instead of a headline QC score.

## Numerical verification

[`scripts/verify_public_results.py`](../scripts/verify_public_results.py) uses scikit-learn's full-SVD PCA and SciPy's z-score/Pearson implementations as independent references. It also reconstructs normalized tables from the original files and compares raw-count VAF to reported VAF. JavaScript is compared recursively with Python under 27 public/synthetic filter scenarios, including singleton and empty selections. Additional small hand-calculated fixtures check exact expected numbers and invalid inputs.

These checks establish numerical agreement and software behavior. A different library is not an independent biological cohort. No external predictive or clinical validation is claimed.
