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
