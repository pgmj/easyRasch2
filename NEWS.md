# easyRasch2 0.1.0

* Initial package structure.
* Added `RMlocdepQ3()` for Yen's Q3 residual correlation analysis of local
  dependence.
* Refactored `RMlocdepQ3()`: the `cutoff` parameter is now optional
  (default `NULL`). When omitted, the raw Q3 matrix is returned without any
  dynamic cut-off; when provided, the caption includes the dynamic cut-off
  value and a recommendation to use `RMlocdepQ3cutoff()`.
* Added `RMlocdepQ3cutoff()` for simulation-based determination of an
  appropriate Q3 cut-off value. Uses CML estimation (`eRm`) and parametric
  bootstrap simulation, with optional parallel processing via **mirai** and an
  optional text progress bar. The recommended workflow is to first run
  `RMlocdepQ3cutoff()` to obtain a simulation-based cutoff (the 99th
  percentile of the max–mean Q3 difference distribution), then pass that value
  to `RMlocdepQ3()`.
