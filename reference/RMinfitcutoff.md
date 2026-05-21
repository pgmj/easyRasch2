# Simulation-Based Infit MSQ Cutoff Determination

Uses parametric bootstrap simulation to determine appropriate cutoff
values for
[`RMiteminfit`](https://pgmj.github.io/easyRasch2/reference/RMiteminfit.md).
This function simulates data from a correctly fitting Rasch model that
mimics your data and returns per-item empirical cutoffs.

## Usage

``` r
RMinfitcutoff(
  data,
  iterations = 250,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.999
)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Only complete cases (rows
  without any `NA`) are used.

- iterations:

  Integer. Number of simulation iterations (default 250).

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

  Character string specifying how cutoff intervals are computed. Either
  `"hdci"` (default) for the Highest Density Interval via
  [`ggdist::hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html),
  or `"quantile"` for the 2.5th/97.5th percentiles via
  [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html).

- hdci_width:

  Numeric. Width of the HDCI when `cutoff_method = "hdci"`. Default is
  `0.999` (99.9% HDCI). Ignored when `cutoff_method = "quantile"`.

## Value

A list with components:

- `results`:

  data.frame with columns `iteration`, `Item`, `InfitMSQ`, `OutfitMSQ`
  (one row per item per successful iteration).

- `item_cutoffs`:

  data.frame with per-item cutoff summaries: `Item`, `infit_low`,
  `infit_high`, `outfit_low`, `outfit_high`. Bounds are computed using
  the method specified by `cutoff_method`.

- `actual_iterations`:

  Number of successful iterations.

- `sample_n`:

  Number of complete cases used.

- `sample_summary`:

  Summary statistics of estimated person parameters.

- `item_names`:

  Character vector of item names from data.

- `cutoff_method`:

  The method used to compute cutoffs (`"hdci"` or `"quantile"`).

- `hdci_width`:

  The HDCI width used (only meaningful when `cutoff_method = "hdci"`).

## Details

For each simulation iteration, person parameters (thetas) are resampled
with replacement from ML estimates, response data are simulated under
the Rasch model, the model is refitted, and conditional infit and outfit
MSQ statistics are computed via
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html). The
distribution of these statistics across iterations provides empirical
critical values per item. Failed iterations (e.g., due to convergence
issues or degenerate data) are silently discarded.

Supports both **dichotomous** data (via
[`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html) and
[`psychotools::rrm()`](https://rdrr.io/pkg/psychotools/man/rrm.html))
and **polytomous** data (via
[`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html) and an internal
partial credit score simulator).

Parallel processing is provided by the `mirai` package (optional).
Install it with `install.packages("mirai")` to enable parallelisation.

The `iarm` package must be installed (it is in Suggests, not Imports).

## See also

[`RMiteminfit`](https://pgmj.github.io/easyRasch2/reference/RMiteminfit.md)

## Examples

``` r
# \donttest{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)

# Run 100 iterations sequentially for a quick demo
cutoff_res <- RMinfitcutoff(sim_data, iterations = 100, parallel = FALSE,
                            seed = 42)
cutoff_res$item_cutoffs
#>      Item infit_low infit_high outfit_low outfit_high
#> 1   Item1     0.902      1.152      0.851       1.324
#> 2   Item2     0.873      1.124      0.844       1.161
#> 3   Item3     0.872      1.101      0.822       1.181
#> 4   Item4     0.883      1.131      0.843       1.186
#> 5   Item5     0.884      1.122      0.870       1.191
#> 6   Item6     0.878      1.172      0.808       1.252
#> 7   Item7     0.863      1.132      0.813       1.237
#> 8   Item8     0.922      1.095      0.864       1.173
#> 9   Item9     0.862      1.105      0.866       1.169
#> 10 Item10     0.901      1.129      0.877       1.181

# Use the cutoffs in RMiteminfit()
RMiteminfit(sim_data)
#> 
#> 
#> Table: MSQ values based on conditional estimation (n = 200 complete cases).
#> 
#> |Item   | Infit MSQ| Relative location|
#> |:------|---------:|-----------------:|
#> |Item1  |     1.008|             -0.22|
#> |Item2  |     0.999|              0.14|
#> |Item3  |     0.994|              0.00|
#> |Item4  |     1.055|             -0.04|
#> |Item5  |     1.004|              0.10|
#> |Item6  |     1.032|             -0.08|
#> |Item7  |     0.943|             -0.04|
#> |Item8  |     0.988|             -0.02|
#> |Item9  |     0.933|              0.12|
#> |Item10 |     1.044|              0.39|
# }
```
