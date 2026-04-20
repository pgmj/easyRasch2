# Simulation-Based Infit MSQ Cutoff Determination for Multiply Imputed Data

Extends
[`RMinfitcutoff`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md)
to work with multiply imputed datasets produced by the `mice` package.
Runs the parametric bootstrap simulation on each imputed dataset and
stacks the resulting distributions, so that the final cutoff intervals
reflect both sampling variability and imputation uncertainty.

## Usage

``` r
RMinfitcutoff_mi(
  mids_object,
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

- mids_object:

  A `mids` object (multiply imputed dataset) as returned by
  `mice::mice()`. Each completed dataset must contain only the item
  response columns to be analysed (i.e., no ID or grouping variables).
  Items must be scored starting at 0 (non-negative integers).

- iterations:

  Integer. Total number of simulation iterations to run across all
  imputations. These are distributed approximately evenly across the `m`
  imputed datasets (default 250).

- parallel:

  Logical. Use parallel processing via `mirai` within each imputed
  dataset (default `TRUE`). Passed to
  [`RMinfitcutoff`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md).

- n_cores:

  Integer or `NULL`. Number of parallel workers. Passed to
  [`RMinfitcutoff`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md).

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
[`RMinfitcutoff`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md),
so that the result can be passed directly to
[`RMiteminfit`](https://pgmj.github.io/easyRasch2/reference/RMiteminfit.md),
[`RMiteminfit_mi`](https://pgmj.github.io/easyRasch2/reference/RMiteminfit_mi.md),
and
[`RMinfitcutoffPlot`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoffPlot.md):

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
`mice::complete()`, then calls
[`RMinfitcutoff`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md)
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

[`RMinfitcutoff`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md),
[`RMiteminfit_mi`](https://pgmj.github.io/easyRasch2/reference/RMiteminfit_mi.md),
[`RMinfitcutoffPlot`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoffPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(mice)

# Create example data with missing values
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
)
colnames(sim_data) <- paste0("Item", 1:8)
# Introduce ~10% MCAR missingness
sim_data[sample(length(sim_data), 0.10 * length(sim_data))] <- NA

# Impute
imp <- mice(sim_data, m = 5, method = "polr", seed = 123, printFlag = FALSE)

# Compute simulation-based cutoffs across imputations
cutoff_mi <- RMinfitcutoff_mi(imp, iterations = 250, parallel = FALSE,
                              seed = 42)
cutoff_mi$item_cutoffs

# Use with RMiteminfit_mi()
RMiteminfit_mi(imp, cutoff = cutoff_mi)
} # }
```
