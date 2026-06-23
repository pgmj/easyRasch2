# Simulation-Based Infit MSQ Cutoff Determination

Uses parametric bootstrap simulation to determine appropriate cutoff
values for
[`RMitemInfit`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md).
This function simulates data from a correctly fitting Rasch model that
mimics your data and returns per-item empirical cutoffs.

## Usage

``` r
RMitemInfitCutoff(
  data,
  iterations = 250,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.999,
  dgp = c("resample", "conditional")
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

- dgp:

  Character. Data-generating process for the parametric bootstrap.
  `"resample"` (default) resamples WLE person locations with replacement
  and simulates responses under the model (a *marginal* null).
  `"conditional"` simulates each respondent's pattern from the exact
  Rasch conditional distribution given their observed total score, item
  parameters fixed (a *conditional* null). Because the conditional
  infit/outfit statistic is itself conditional on the total score,
  `"conditional"` is its naturally matched null. **Experimental.**

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

- `dgp`:

  The data-generating process used (`"resample"` or `"conditional"`).

## Details

The generating model is CML item parameters (via `psychotools`) with WLE
person locations. For each iteration a dataset is simulated under the
chosen `dgp`, the model is refitted by CML
([`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html),
which handles dichotomous and polytomous data and is accepted by
`iarm`), and conditional infit and outfit MSQ are computed via
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html). The
distribution of these statistics across iterations provides empirical
critical values per item. Failed iterations (e.g., degenerate simulated
data) are silently discarded.

Parallel processing is provided by the `mirai` package (optional).
Install it with `install.packages("mirai")` to enable parallelisation.

The `iarm` package must be installed (it is in Suggests, not Imports).

## See also

[`RMitemInfit`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)

## Examples

``` r
# \donttest{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)

# Run 100 iterations sequentially for a quick demo
cutoff_res <- RMitemInfitCutoff(sim_data, iterations = 100, parallel = FALSE,
                            seed = 42)
cutoff_res$item_cutoffs
#>      Item infit_low infit_high outfit_low outfit_high
#> 1   Item1     0.819      1.180      0.787       1.360
#> 2   Item2     0.873      1.133      0.813       1.229
#> 3   Item3     0.872      1.139      0.819       1.210
#> 4   Item4     0.865      1.122      0.814       1.180
#> 5   Item5     0.868      1.125      0.835       1.165
#> 6   Item6     0.837      1.121      0.786       1.185
#> 7   Item7     0.883      1.147      0.828       1.248
#> 8   Item8     0.873      1.117      0.836       1.205
#> 9   Item9     0.890      1.157      0.829       1.270
#> 10 Item10     0.873      1.112      0.788       1.233

# Use the cutoffs in RMitemInfit()
RMitemInfit(sim_data)
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
#> |Item10 |     1.044|              0.38|
# }
```
