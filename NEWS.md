# easyRasch2 0.1.0

* Initial package structure.
* Added `RMlocdepQ3()` for Yen's Q3 residual correlation analysis of local dependence.

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
