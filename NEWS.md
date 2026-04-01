# easyRasch2 0.1.0

* Initial package structure.
* Added `RMlocdepQ3()` for Yen's Q3 residual correlation analysis of local dependence.
* Refactored `RMlocdepQ3()`: `cutoff` is now optional (default `NULL`). When
  `NULL`, the raw Q3 matrix is returned without any cut-off annotation. When
  provided, a dynamic cut-off (mean Q3 + cutoff) is shown in the caption.
* Added `RMlocdepQ3cutoff()`: simulation-based (parametric bootstrap) cutoff
  determination for `RMlocdepQ3()`. Supports both dichotomous and polytomous
  data. Uses `foreach` + `doParallel` + `doRNG` for reproducible parallel
  simulation. Returns a list including percentile cutoffs (p95–p999) and a
  `suggested_cutoff` (95th percentile) ready for use with `RMlocdepQ3()`.
* Recommended workflow: first run `RMlocdepQ3cutoff()` to obtain a
  simulation-based cutoff, then pass `suggested_cutoff` to `RMlocdepQ3()`.
