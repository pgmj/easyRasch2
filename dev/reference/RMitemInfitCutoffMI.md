# Simulation-Based Infit MSQ Cutoff Determination for Multiply Imputed Data

Extends
[`RMitemInfitCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md)
to work with multiply imputed datasets produced by the `mice` package.
Runs the parametric bootstrap simulation on each imputed dataset and
stacks the resulting distributions, so that the final cutoff intervals
reflect both sampling variability and imputation uncertainty.

## Usage

``` r
RMitemInfitCutoffMI(
  mids_object,
  iterations = 500,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.999
)
```

## Arguments

- mids_object:

  A `mids` object (multiply imputed dataset) as returned by
  [`mice::mice()`](https://amices.org/mice/reference/mice.html). Each
  completed dataset must contain only the item response columns to be
  analysed (i.e., no ID or grouping variables). Items must be scored
  starting at 0 (non-negative integers).

- iterations:

  Integer. Total number of simulation iterations to run across all
  imputations. These are distributed approximately evenly across the `m`
  imputed datasets (default 500).

- parallel:

  Logical. Use parallel processing via `mirai` within each imputed
  dataset (default `TRUE`). Passed to
  [`RMitemInfitCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md).

- n_cores:

  Integer or `NULL`. Number of parallel workers. Passed to
  [`RMitemInfitCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md).

- verbose:

  Logical. Show progress messages (default `FALSE`).

- seed:

  Integer or `NULL`. Master random seed for reproducibility. A unique
  per-imputation seed is derived from this value.

- cutoff_method:

  Character string specifying how cutoff intervals are computed from the
  stacked distribution. Either `"hdci"` (default) for the Highest
  Density Interval via
  [`ggdist::hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html),
  or `"quantile"` for the 2.5th/97.5th percentiles via
  [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html).

- hdci_width:

  Numeric. Width of the HDCI when `cutoff_method = "hdci"`. Default is
  `0.999` (99.9% HDCI). Ignored when `cutoff_method = "quantile"`.

## Value

A list with the same structure as
[`RMitemInfitCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md),
so that the result can be passed directly to
[`RMitemInfit`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md),
[`RMitemInfitMI`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md),
and
[`RMitemInfitPlot`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md):

- `results`:

  data.frame with columns `iteration`, `imputation`, `Item`, `InfitMSQ`,
  `OutfitMSQ` — the stacked simulation results from all imputed
  datasets.

- `item_cutoffs`:

  data.frame with per-item cutoff summaries: `Item`, `infit_low`,
  `infit_high`, `outfit_low`, `outfit_high`. Computed from the stacked
  distribution.

- `actual_iterations`:

  Total number of successful iterations across all imputations.

- `sample_n`:

  Number of rows (respondents) per imputed dataset.

- `sample_summary`:

  Summary statistics of estimated person parameters from the first
  imputed dataset.

- `item_names`:

  Character vector of item names.

- `cutoff_method`:

  The method used to compute cutoffs.

- `hdci_width`:

  The HDCI width used.

- `n_imputations`:

  Number of imputed datasets used.

- `iterations_per_imputation`:

  Integer vector of requested iterations per imputed dataset.

- `actual_iterations_per_imputation`:

  Integer vector of successful iterations per imputed dataset.

## Details

The function completes each of the `m` imputed datasets via
[`mice::complete()`](https://tidyr.tidyverse.org/reference/complete.html),
then calls
[`RMitemInfitCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md)
on each one. The total number of iterations is split approximately
evenly across imputations (i.e., each imputed dataset receives
`ceiling(iterations / m)` or `floor(iterations / m)` iterations). The
per-imputation simulation results are stacked into a single distribution
from which cutoff intervals are computed, naturally incorporating
imputation uncertainty.

Imputed datasets that cause model convergence failures are dropped with
a warning. If all imputations fail, the function stops with an error.

The `mice` package must be installed (it is in Suggests, not Imports).

## See also

[`RMitemInfitCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md),
[`RMitemInfitMI`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md),
[`RMitemInfitPlot`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md)

## Examples

``` r
# \donttest{
if (requireNamespace("mice", quietly = TRUE) &&
    requireNamespace("iarm", quietly = TRUE) &&
    requireNamespace("ggdist", quietly = TRUE)) {
  # Create example data with ~10% MCAR missingness
  set.seed(42)
  mat <- matrix(sample(0:1, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
  mat[sample(length(mat), round(0.10 * length(mat)))] <- NA
  sim_data <- as.data.frame(mat)
  colnames(sim_data) <- paste0("Item", 1:8)

  # mice's ordinal method (`polr`) requires the items to be ordered
  # factors, so code them as such before imputing. RMitemInfitCutoffMI()
  # converts the completed factors back to numeric internally.
  sim_data[] <- lapply(sim_data, function(x) factor(x, ordered = TRUE))

  # Impute (use more imputations, e.g. m = 5+, in real analyses)
  imp <- mice::mice(sim_data, m = 2, method = "polr", seed = 123,
                    printFlag = FALSE)

  # Compute simulation-based cutoffs across imputations
  # (use more iterations, e.g. 250+, in real analyses)
  cutoff_mi <- RMitemInfitCutoffMI(imp, iterations = 50, parallel = FALSE,
                                seed = 42)
  cutoff_mi$item_cutoffs

  # Use with RMitemInfitMI()
  RMitemInfitMI(imp, cutoff = cutoff_mi)
}
#> 
#> 
#> Table: Pooled MSQ values from 2 imputations (Rubin's rules). n = 200 respondents (missing values imputed, m = 2) per imputed dataset. Cutoff values based on 50 total simulation iterations across 2 imputations (99.9% HDCI). Flagged: overfit = infit below range (more predictable); underfit = above range (noisier).
#> 
#> |Item  | Infit MSQ| Infit SE| Infit low| Infit high|Flagged | Relative location|
#> |:-----|---------:|--------:|---------:|----------:|:-------|-----------------:|
#> |Item1 |     0.950|    0.063|     0.889|      1.160|        |             -0.25|
#> |Item2 |     1.022|    0.050|     0.889|      1.094|        |              0.05|
#> |Item3 |     1.008|    0.055|     0.835|      1.138|        |             -0.01|
#> |Item4 |     1.034|    0.058|     0.880|      1.144|        |             -0.01|
#> |Item5 |     0.995|    0.059|     0.873|      1.160|        |              0.10|
#> |Item6 |     1.016|    0.053|     0.825|      1.124|        |             -0.08|
#> |Item7 |     0.942|    0.050|     0.858|      1.190|        |             -0.08|
#> |Item8 |     1.030|    0.050|     0.869|      1.132|        |              0.02|
# }
```
