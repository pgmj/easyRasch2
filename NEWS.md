# easyRasch2 0.8.0

## Function renaming for consistency

This release renames **22 exported functions** under a consistent
*domain-prefix → method → variant-suffix* scheme so that autocompletion
on a prefix (`RMdif`, `RMlocdep`, `RMitem`, `RMdim`, `RMplot`) surfaces
every related function. No semantic changes; only names.

A complete mapping is on the new `?easyRasch2-renaming` help page. The
short version:

| Domain          | Prefix       | Examples |
|-----------------|--------------|----------|
| DIF             | `RMdif`      | `RMdifGamma()`, `RMdifGammaCutoff()`, `RMdifGammaPlot()`, `RMdifLR()`, `RMdifTree()` |
| Local dependence| `RMlocdep`   | `RMlocdepQ3()`, `RMlocdepQ3Cutoff()`, `RMlocdepQ3Plot()`, `RMlocdepGamma()`, `RMlocdepGammaCutoff()`, `RMlocdepGammaPlot()` |
| Item statistics | `RMitem`     | `RMitemInfit()`, `RMitemInfitCutoff()`, `RMitemInfitCutoffPlot()`, `RMitemInfitMI()`, `RMitemInfitCutoffMI()`, `RMitemRestscore()`, `RMitemRestscoreBoot()`, `RMitemICCPlot()`, `RMitemHierarchy()` |
| Dimensionality  | `RMdim`      | `RMdimResidualPCA()`, `RMdimResidualPCACutoff()`, `RMdimCFACutoff()`, `RMdimCFAPlot()`, `RMdimMartinLof()`, `RMdimMartinLofResiduals()` |
| Descriptive plot| `RMplot`     | `RMplotTile()`, `RMplotBar()`, `RMplotStackedbar()` |

Unchanged: `RMreliability()`, `RMUreliability()`, `RMtargeting()`,
`RMscoreSE()`.

**Breaking change.** No deprecation aliases are shipped — existing
scripts will need a search-and-replace. See `?easyRasch2-renaming` for
the full table.

# easyRasch2 0.7.1

- New **`RMlocdepQ3plot()`** for visualising per-pair simulated Q3
  distributions against observed Q3 values, mirroring the design of
  `RMpgLDplot()` and `RMpgDIFplot()`. Supports `items` and `n_pairs`
  filters; with `data` supplied, ranks `n_pairs` by
  `|observed Q3 - median(simulated Q3 per pair)|`.
- `RMlocdepQ3cutoff()` now retains per-pair iteration-level data:
  the returned list gains `pair_results`, `pair_cutoffs`, `item_names`,
  `cutoff_method`, and `hdci_width`. Existing scalar outputs
  (`suggested_cutoff` etc.) are unchanged. New arguments `cutoff_method`
  (`"hdci"` / `"quantile"`) and `hdci_width` (default `0.99`) control
  how per-pair credible intervals are computed.
- `RMlocdepQ3()` can now be called with the full list returned by
  `RMlocdepQ3cutoff()` (it auto-extracts `$suggested_cutoff`), matching
  the partial-gamma family convention. Numeric scalar still works.
- `RMpartgamLD()` and `RMpgLDplot()` gain an `n_pairs` argument
  to keep only the top-N pairs by `|gamma|` (LD) or by
  `|observed - median(sim)|` deviation (plot).
- `RMpartgamLD()`, `RMpartgamDIF()`, and `RMitemrestscore()`
  now use consistent column headers: `Adj. p-value (BH)` and a new
  `p-value sign.` star-string column. `RMitemrestscore()`'s former
  `Absolute_difference` column is now a *signed* `Difference`
  (observed − expected), so over- and underfit are visually
  distinguishable at a glance.


# easyRasch2 0.7.0

- New **`RMbarplot()`** and **`RMstackedbarplot()`** join `RMtileplot()` for 
  easy visualization of item response data distributions.
- New **`RMciccPlot()`** for conditional item characteristic curves plot, also 
  includes DIF analysis for categorical DIF variables.
- New **`RMitemHierarchy()`**, which outputs a plot illustrating the item hierarchy
  with item thresholds and confidence intervals.
- Modified `RMcfaCutoff()` to only use .scaled metrics for RMSEA and CFI, for stability.
  See documentation for more details.
- Bug fix for `RMmartinLof()` with dichotomous items.
- Bug fix for `RMitemHierarchy()` making option `item_labels` work as intended. 
- Tests written for all functions, preparing for CRAN submission in the nearish
  future.
- Drafted a package intro/vignette

# easyRasch2 0.6.0

- New **`RMcfaCutoff()`** for testing unidimensionality against simulation-based
  cutoff values for model fit metrics.
  - **`RMcfaPlot()`** shows figure with distribution of simulation results and
  observed model fit values.

# easyRasch2 0.5.8

- New **`RMdifTree()`** for testing DIF with continuous variables, such as age
  in years, as well as combinations of DIF variables (DIF interactions).
  - Also implements use of `stablelearner`, as recommended by Henninger et al.
    (2025) to assess stability of tree-based results.
  - Borrows from
    * raschtreeMH (https://github.com/mirka-henninger/raschtreeMH)
    * effecttree  (https://github.com/mirka-henninger/effecttree)
  - See these two papers (and `?RMdifTree`) for more details:
    * Henninger, M., Debelak, R., & Strobl, C. (2023). A new stopping
      criterion for Rasch trees based on the Mantel-Haenszel effect size
      measure for DIF. Educational Psychological Measurement, 83, 181-212.
      doi:10.1177/00131644221077135
    * Henninger, M., Radek, J., Debelak, R., & Strobl, C. (2025). Partial
      credit trees meet the partial gamma coefficient for quantifying DIF
      and DSF in polytomous items. Behaviormetrika, 52, 221-257.
  
# easyRasch2 0.5.7

- New **`RMdifLR()`** for testing DIF of categorical variables using Andersens's
  Likelihood Ratio test as implemented in package `eRm`
- Visualization of response data using a tile plot with **`RMtileplot()`**.
  - Optional `group` faceting, which is useful for DIF-analyses.
- Added Martin-Löf test of dimensionality **`RMmartinLof()`** based on Christensen & Kreiner 
  (2007, doi: 10.1177/0146621605286204), for both dichotomous and polytomous (PCM) items.
  - Also added post-processing function **`RMmartinLofresiduals()`**.

# easyRasch2 0.5.6

- Added **`RMresidualPCA()`** for evaluating patterns in the standardized residuals 
  from the RM/PCM. Outputs either a table with eigenvalues and explained variance
  or a figure with standardized loadings on the first residual contrast and item
  locations.
- Added **`RMpcaCutoff()`** to determine a simulation-based critical value for the
  largest eigenvalue.
 
# easyRasch2 0.5.5

- Added **`RMbootRestscore()`** for use with large sample sizes. See <https://pgmj.github.io/rasch_itemfit/> for more details.
- Added **`RMscoreSE()`** that produces a transformation table (or figure) from ordinal sum
  scores to WLE (Weighted Likelihood Estimation) interval scores.
- Modified `RMlocdepQ3cutoff()` to improve speed.

# easyRasch2 0.5.4

- **`RMtargeting()`** for Wright map style plots. 
- New functions for handling missing data with conditional item infit MSQ by multiple imputation using package `mice`.
  - **`RMinfitcutoff_mi()`** uses the `mids` object containing multiple imputated datasets (output by `mice::mice()`) to run simulations on each datasets and combines the results. Splits iterations across imputations (e.g. 250 total / 5 imputations = 50 each).
  - **`RMiteminfit_mi()`** calculates and pools conditional infit from the imputated datasets and optionally uses the `RMinfitcutoff_mi()` output for cutoff values.
  - Works with existing **`RMinfitcutoffPlot()`**
  
# easyRasch2 0.5.3

* New functions **`RMpartgamLD()`** and **`RMpgLDcutoff()`** for evaluating local dependency of item pairs in both directions.
  * Plot function **`RMpgLDplot()`**, similar to **`RMpgDIFplot()`**
  * NOTE: The simulation-based cutoffs have not yet been evaluated in a systematic way.
  Default choices may not be sensible, please try different numbers of iterations and hdci_width.

# easyRasch2 0.5.2

* New functions **`RMpartgamDIF()`** and **`RMpgDIFcutoff()`** for evaluating DIF of categorical external variables.
  * Also with plot function **`RMpgDIFplot()`**, similar to `RMinfitcutoffPlot()`
  * NOTE: The simulation-based cutoffs have not yet been evaluated in a systematic way.
  Default choices may not be sensible, please try different numbers of iterations and hdci_width.

# easyRasch2 0.5.1

* New function **`RMinfitcutoffPlot()`** to illustrate distribution of simulated
  conditional item infit MSQ values together with the observed value.

# easyRasch2 0.5.0

* **`RMinfitcutoff()` gains `cutoff_method` and `hdci_width` parameters**:
  By default (`cutoff_method = "hdci"`, `hdci_width = 0.999`), per-item cutoff
  intervals are now computed using the Highest Density Interval via
  `ggdist::hdci()` (99.9% HDCI). Set `cutoff_method = "quantile"` to restore the
  previous behaviour (2.5th/97.5th percentiles). The `ggdist` package is only
  required when `cutoff_method = "hdci"` (added to Suggests). The returned list
  now also includes `cutoff_method` and `hdci_width` fields.

* **`RMiteminfit()` caption updated**: When `cutoff` is the return value of
  `RMinfitcutoff()`, the kable caption now states the cutoff method, e.g.
  `"Cutoff values based on 250 simulation iterations (99.9% HDCI)."` or
  `"Cutoff values based on 250 simulation iterations (2.5th/97.5th percentile)."`.


# easyRasch2 0.4.0

* **`RMiteminfit()` gains optional `cutoff` parameter**: Accepts the return
  value of `RMinfitcutoff()` (or its `$item_cutoffs` data.frame directly).
  When provided, per-item cutoff boundaries (`Infit_low`, `Infit_high`) and a
  logical `Flagged` column are added to both `"dataframe"` and `"kable"`
  output. The kable caption includes the number of simulation iterations when
  available.

* **New `RMinfitcutoff()`**: Simulation-based (parametric bootstrap) cutoff
  determination for [RMiteminfit()]. Supports both dichotomous and polytomous
  data. Optional parallel processing via `mirai` (falls back to sequential if
  not installed). Returns per-item 2.5th and 97.5th percentiles of simulated
  infit and outfit MSQ distributions (`$item_cutoffs`), together with the full
  iteration-level results (`$results`). Requires the `iarm` package (Suggests).


# easyRasch2 0.3.0

* **New `RMiteminfit()`**: Computes conditional infit MSQ statistics for each
  item via `iarm::out_infit()`, enriched with item locations relative to the
  sample mean person location. Supports both dichotomous (Rasch model via
  `eRm::RM()`) and polytomous (Partial Credit Model via `eRm::PCM()`) data.
  Requires the `iarm` package (Suggests). Output options: `"kable"` (default,
  plain-text `knitr::kable()`) or `"dataframe"`. Optional `sort = "infit"`
  sorts by infit MSQ descending. Only complete cases are used for conditional
  fit calculation.

* **New `RMitemrestscore()`**: Computes observed and model-expected
  item-restscore correlations via `iarm::item_restscore()`, enriched with
  absolute differences between observed and expected values, item average
  locations, and item locations relative to the sample mean person location.
  Supports both dichotomous (Rasch model via `eRm::RM()`) and polytomous
  (Partial Credit Model via `eRm::PCM()`) data. Requires the `iarm` package
  (Suggests). Output options: `"kable"` (default, plain-text `knitr::kable()`)
  or `"dataframe"`. Optional `sort = "diff"` sorts by absolute difference
  descending.


# easyRasch2 0.2.0

* **`RMlocdepQ3()`**: `cutoff` is now optional (default `NULL`). When omitted,
  the raw Q3 residual correlation matrix is returned without any dynamic
  cut-off applied. When provided, the dynamic cut-off (mean Q3 + cutoff) is
  shown in the kable caption as before. Supersedes the requirement to always
  supply a cutoff value (PRs #2 and #3).
* **New `RMlocdepQ3cutoff()`**: Simulation-based (parametric bootstrap) cutoff
  determination for `RMlocdepQ3()`. Supports both dichotomous and polytomous
  data. Optional parallel processing via `mirai` (falls back to sequential if
  not installed). The `$suggested_cutoff` is the 99th percentile of the
  simulated Q3 max–mean distribution.
* **New internal helpers** in `R/utils-simulation.R`: `sim_poly_item()`,
  `sim_partial_score()`, and `extract_item_thresholds()`.
* **Dependencies**: `eRm`, `psychotools (>= 0.7-3)`, `parallel`, and `utils`
  moved/added to Imports. `mirai` added to Suggests.


# easyRasch2 0.1.0

* Initial package structure.
* Added `RMlocdepQ3()` for Yen's Q3 residual correlation analysis of local dependence.
