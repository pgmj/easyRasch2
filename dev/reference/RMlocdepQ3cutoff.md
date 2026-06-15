# Simulation-Based Q3 Cutoff Determination

Uses parametric bootstrap simulation to determine an appropriate cutoff
value for
[`RMlocdepQ3`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md).
Under a correctly fitting Rasch model, Q3 residuals have an unknown
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
  hdci_width = 0.99
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

  Character. Method used to compute per-pair Q3 credible intervals in
  `pair_cutoffs`. One of `"hdci"` (the default, Highest Density
  Continuous Interval via
  [`ggdist::hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html))
  or `"quantile"` (symmetric 2.5th / 97.5th percentiles). Only affects
  `pair_cutoffs`; the global `$suggested_cutoff` (99th percentile of
  `max(Q3) - mean(Q3)`) is unaffected.

- hdci_width:

  Numeric in (0, 1). Width of the HDCI when `cutoff_method = "hdci"`.
  Default `0.99`. Ignored when `cutoff_method = "quantile"`.

## Value

A list with components:

- `results`:

  data.frame with columns `mean`, `max`, `diff` (one row per successful
  iteration).

- `pair_results`:

  Long data.frame with columns `Item1`, `Item2`, `Q3`, `iteration` — one
  row per item pair per successful iteration. Used by
  [`RMlocdepQ3Plot`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md).

- `pair_cutoffs`:

  data.frame with per-pair cutoff summaries: `Item1`, `Item2`, `Q3_low`,
  `Q3_high`. Boundaries are computed via the method specified by
  `cutoff_method`.

- `actual_iterations`:

  Number of successful iterations.

- `sample_n`:

  Number of persons in the original data.

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
  [`RMlocdepQ3`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md).

- `cutoff_method`:

  The method used for `pair_cutoffs` (`"hdci"` or `"quantile"`).

- `hdci_width`:

  The HDCI width used (only meaningful when `cutoff_method = "hdci"`).

## Details

For each simulation iteration, person parameters (thetas) are resampled
with replacement from ML estimates, response data are simulated under
the Rasch model, a `mirt` Rasch model is fitted to the simulated data,
and Q3 residuals are extracted. The distribution of `max(Q3) - mean(Q3)`
across iterations provides empirical critical values. Failed iterations
(e.g., due to convergence issues) are silently discarded.

Supports both **dichotomous** data (via
[`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html) and
[`psychotools::rrm()`](https://rdrr.io/pkg/psychotools/man/rrm.html))
and **polytomous** data (via
[`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html) and an internal
partial credit score simulator).

Parallel processing is provided by the `mirai` package (optional).
Install it with `install.packages("mirai")` to enable parallelisation.

## See also

[`RMlocdepQ3`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)

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
#> Table: Dynamic cut-off: 0.237 (mean Q3 = -0.012 + 0.249). Correlations exceeding the cut-off may indicate local dependence.
#> 
#> |       |Item1 |Item2 |Item3 |Item4 |Item5 |Item6 |Item7 |Item8 |Item9 |Item10 |above_cutoff |
#> |:------|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:------|:------------|
#> |Item1  |      |      |      |      |      |      |      |      |      |       |             |
#> |Item2  |0.08  |      |      |      |      |      |      |      |      |       |             |
#> |Item3  |-0.07 |-0.09 |      |      |      |      |      |      |      |       |             |
#> |Item4  |-0.07 |-0.12 |-0.03 |      |      |      |      |      |      |       |             |
#> |Item5  |-0.09 |0.02  |0.08  |-0.04 |      |      |      |      |      |       |             |
#> |Item6  |-0.03 |-0.02 |-0.07 |0.09  |-0.1  |      |      |      |      |       |             |
#> |Item7  |-0.05 |-0.03 |0.04  |0.13  |0.08  |0     |      |      |      |       |             |
#> |Item8  |0.07  |0.05  |0     |-0.19 |0.05  |-0.06 |-0.05 |      |      |       |             |
#> |Item9  |0.08  |0.05  |-0.02 |-0.07 |-0.09 |0.07  |0.04  |0.12  |      |       |             |
#> |Item10 |-0.07 |-0.05 |0.08  |-0.04 |-0.02 |-0.12 |-0.04 |-0.03 |0     |       |             |
# }
```
