# Simulation-based infit MSQ cutoffs for multiply imputed data

The counterpart of
[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md)
for data imputed with the mice package. Runs the simulations in each
imputed dataset and stacks them, so that the interval per item reflects
both sampling variability and the uncertainty from imputation. The
result is used by
[`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md).

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

  A `mids` object from
  [`mice::mice()`](https://amices.org/mice/reference/mice.html). Each
  completed dataset must contain only the item response columns, with no
  ID or grouping variables. Items are scored from 0. Items imputed as
  factors, as the `polr` method needs, are converted back to numbers.
  Their levels must be the response codes (`"0"`, `"1"`, ...).

- iterations:

  Integer. Total number of simulated datasets across all imputations,
  split about evenly between the `m` imputed datasets. Default `500`.
  See the section on interval width.

- parallel:

  Logical. Run the simulations within each imputed dataset in parallel
  via mirai. Default `TRUE`. Passed to
  [`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md).

- n_cores:

  Integer or `NULL`. Number of parallel workers. Passed to
  [`RMitemInfitCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md).

- verbose:

  Logical. Show progress messages (default `FALSE`).

- seed:

  Integer or `NULL`. Random seed for reproducibility, from which a seed
  for each imputation is derived. See
  [easyRasch2-reproducibility](https://pgmj.github.io/easyRasch2/dev/reference/easyRasch2-reproducibility.md)
  for what this guarantees and how it interacts with `parallel`.

- cutoff_method:

  Character. How the interval per item is computed from the stacked
  simulations: `"hdci"` (default), the highest-density continuous
  interval from
  [`ggdist::hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html)
  (needs ggdist), or `"quantile"`, the 2.5th to 97.5th percentiles.

- hdci_width:

  Numeric. Width of the interval with `cutoff_method = "hdci"`. Default
  `0.999`, see the section on interval width. Ignored with `"quantile"`.

## Value

A list like that of
[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md),
which
[`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md)
and
[`RMitemInfitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md)
accept, with:

- `results`:

  data.frame with the columns `iteration`, `imputation`, `Item`,
  `InfitMSQ` and `OutfitMSQ`, the stacked simulations from all imputed
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

Each of the `m` imputed datasets is completed with
[`mice::complete()`](https://tidyr.tidyverse.org/reference/complete.html)
and passed to
[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md),
with `ceiling(iterations / m)` or `floor(iterations / m)` iterations
each. The simulations from all imputations are stacked, and the interval
per item is computed from the stack.

Imputed datasets in which the model cannot be fitted are dropped with a
warning. If all imputations fail, the function stops with an error.

The mice package must be installed (it is in Suggests).

## Interval width under multiple imputation

`hdci_width` defaults to `0.999` here, not to the `0.95` of
[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md).
[`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md)
has no corrected p-value yet, so under imputation the interval is still
the decision rule, and it keeps a width chosen for that. Its family-wise
error rate is `1 - 0.999^k` over `k` items. An interval this wide needs
roughly 5000 iterations in all before it reaches its stated width
(Johansson, 2026). Combining bootstrap p-values across imputations is
planned but needs its own calibration, since Johansson (2026) covers
complete data only.

## References

Johansson, M. (2026). Simulation-based cutoffs for conditional item fit
in Rasch models: Iterations, multiplicity correction, and decision
stability. *PsyArXiv*.
[doi:10.31234/osf.io/7pqz4_v2](https://doi.org/10.31234/osf.io/7pqz4_v2)

## See also

[`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md)
for the pooled test,
[`RMitemInfitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md)
to plot the simulated values,
[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md)
for complete data.

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

  # Simulation-based cutoffs across imputations (use far more
  # iterations in real analyses, see the section on interval width)
  cutoff_mi <- RMitemInfitCutoffMI(imp, iterations = 50, parallel = FALSE,
                                seed = 42)
  cutoff_mi$item_cutoffs

  # Use with RMitemInfitMI()
  RMitemInfitMI(imp, cutoff = cutoff_mi)
}
#> 
#> 
#> Table: Pooled MSQ values from 2 imputations (Rubin's rules). n = 200 respondents (missing values imputed, m = 2) per imputed dataset. Cutoff values based on 50 total simulation iterations across 2 imputations (99.9% HDCI). Flagged: overfit = infit below range; underfit = above range.
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
