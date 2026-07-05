# Simulation-Based \\Q_3\\ Cutoff Determination

Uses parametric bootstrap simulation to determine an appropriate cutoff
value for
[`RMlocdepQ3`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3.md).
Under a correctly fitting Rasch model, \\Q_3\\ residuals have an unknown
distribution; this function simulates that distribution and returns
empirical percentiles.

## Usage

``` r
RMlocdepQ3Cutoff(
  data,
  iterations = 500,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.99,
  estimator = c("CML", "MML"),
  dgp = c("resample", "conditional")
)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers).

- iterations:

  Integer. Number of simulation iterations (default 500).

- parallel:

  Logical. Use parallel processing via `mirai` if available (default
  `TRUE`).

- n_cores:

  Integer or `NULL`. Number of parallel workers. When `NULL`,
  `getOption("mc.cores")` is checked first. If neither is set and
  `parallel = TRUE`, a warning is issued and execution falls back to
  sequential (single core) processing.

- verbose:

  Logical. Show a progress bar (default `FALSE`).

- seed:

  Integer or `NULL`. Random seed for reproducibility.

- cutoff_method:

  Character. Method used to compute per-pair \\Q_3\\ credible intervals
  in `pair_cutoffs`. One of `"hdci"` (the default, Highest Density
  Continuous Interval via
  [`ggdist::hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html))
  or `"quantile"` (symmetric 2.5th / 97.5th percentiles). Only affects
  `pair_cutoffs`; the global `$suggested_cutoff` (99th percentile of
  `max(Q3) - mean(Q3)`) is unaffected.

- hdci_width:

  Numeric in (0, 1). Width of the HDCI when `cutoff_method = "hdci"`.
  Default `0.99`. Ignored when `cutoff_method = "quantile"`.

- estimator:

  Character. Estimation engine for the simulated \\Q_3\\ values, passed
  through to the per-iteration computation. `"CML"` (default) uses CML
  item parameters and WLE person locations; `"MML"` uses `mirt`. This
  must match the `estimator` later given to
  [`RMlocdepQ3`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3.md);
  the value is stored in the returned object's `$estimator` and reused
  automatically.

- dgp:

  Character. Data-generating process for the parametric bootstrap.
  `"resample"` (default) draws person locations by resampling the WLE
  estimates with replacement and simulates responses under the model – a
  *marginal* null. `"conditional"` instead simulates each respondent's
  pattern from the exact Rasch conditional distribution given their
  observed total score (and answered items), with item parameters fixed
  – a *conditional* null that fixes the score margin and needs no latent
  distribution, avoiding the over-dispersion of resampled point
  estimates. The two give different cut-offs; see the package's
  comparison study. **Experimental.**

## Value

A list with components:

- `results`:

  data.frame with columns `mean`, `max`, `diff` (one row per successful
  iteration).

- `pair_results`:

  Long data.frame with columns `Item1`, `Item2`, `Q3`, `iteration` — one
  row per item pair per successful iteration. Used by
  [`RMlocdepQ3Plot`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3plot.md).

- `pair_cutoffs`:

  data.frame with per-pair cutoff summaries: `Item1`, `Item2`, `Q3_low`,
  `Q3_high`. Boundaries are computed via the method specified by
  `cutoff_method`.

- `actual_iterations`:

  Number of successful iterations.

- `sample_n`:

  Number of persons in the original data.

- `sample_n_total`:

  Equal to `sample_n`: no respondents are dropped (incomplete responses
  are retained). Stored for consistency with the other `*Cutoff()`
  objects.

- `sample_has_na`:

  Logical. Whether the data contained any missing values.

- `sample_summary`:

  Summary statistics of estimated person parameters.

- `item_names`:

  Character vector of item names from `data`.

- `max_diff`, `sd_diff`:

  Max and SD of the `diff` distribution.

- `p95`, `p99`, `p995`, `p999`:

  Empirical percentiles of `diff`.

- `suggested_cutoff`:

  The 99th percentile (`p99`) — recommended scalar cutoff for
  [`RMlocdepQ3`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3.md).

- `cutoff_method`:

  The method used for `pair_cutoffs` (`"hdci"` or `"quantile"`).

- `hdci_width`:

  The HDCI width used (only meaningful when `cutoff_method = "hdci"`).

- `estimator`:

  The estimator used for the simulated \\Q_3\\ (`"CML"` or `"MML"`);
  reused by
  [`RMlocdepQ3`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3.md)
  and
  [`RMlocdepQ3Plot`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3plot.md).

- `dgp`:

  The data-generating process used (`"resample"` or `"conditional"`).

## Details

The generating model is fitted once: CML item parameters (via
`psychotools`) and WLE person locations. For each simulation iteration,
those WLE thetas are resampled with replacement, response data are
simulated under the Rasch / Partial Credit model, the model is refitted,
and \\Q_3\\ residuals are computed under `estimator`. The distribution
of `max(Q3) - mean(Q3)` across iterations provides empirical critical
values. Failed iterations (e.g., due to convergence issues) are silently
discarded.

Supports both **dichotomous** data (simulated via
[`psychotools::rrm()`](https://rdrr.io/pkg/psychotools/man/rrm.html))
and **polytomous** data (via an internal partial credit score
simulator).

Parallel processing is provided by the `mirai` package (optional).
Install it with `install.packages("mirai")` to enable parallelisation.

## See also

[`RMlocdepQ3`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3.md)

## Examples

``` r
# \donttest{
if (requireNamespace("ggdist", quietly = TRUE)) {
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
  )
  colnames(sim_data) <- paste0("Item", 1:10)

  # Few iterations for a fast example; use 500+ in real analyses
  cutoff_res <- RMlocdepQ3Cutoff(sim_data, iterations = 50, parallel = FALSE,
                                 seed = 42)
  cutoff_res$suggested_cutoff  # 99th percentile

  # Use the cutoff in RMlocdepQ3()
  RMlocdepQ3(sim_data, cutoff = cutoff_res$suggested_cutoff)
}
#> 
#> 
#> Table: Dynamic cut-off: 0.112 (mean Q3 -0.11 + 0.222). Correlations exceeding the cut-off may indicate local dependence; see the per-pair table for detail. n = 200 respondents.
#> 
#> |       |Item1 |Item2 |Item3 |Item4 |Item5 |Item6 |Item7 |Item8 |Item9 |Item10 |above_cutoff |
#> |:------|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:------|:------------|
#> |Item1  |      |      |      |      |      |      |      |      |      |       |             |
#> |Item2  |-0.02 |      |      |      |      |      |      |      |      |       |             |
#> |Item3  |-0.17 |-0.2  |      |      |      |      |      |      |      |       |             |
#> |Item4  |-0.13 |-0.2  |-0.1  |      |      |      |      |      |      |       |             |
#> |Item5  |-0.19 |-0.06 |-0.02 |-0.12 |      |      |      |      |      |       |             |
#> |Item6  |-0.1  |-0.1  |-0.17 |0.03  |-0.19 |      |      |      |      |       |             |
#> |Item7  |-0.17 |-0.13 |-0.07 |0.03  |-0.04 |-0.11 |      |      |      |       |             |
#> |Item8  |-0.03 |-0.06 |-0.11 |-0.28 |-0.05 |-0.16 |-0.18 |      |      |       |             |
#> |Item9  |-0.03 |-0.08 |-0.16 |-0.19 |-0.22 |-0.03 |-0.1  |-0.02 |      |       |             |
#> |Item10 |-0.16 |-0.14 |0     |-0.09 |-0.1  |-0.18 |-0.15 |-0.12 |-0.09 |       |             |
# }
```
