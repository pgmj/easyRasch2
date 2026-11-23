# Simulation-based infit MSQ cutoffs

Simulates datasets from a Rasch or partial credit model that fits and
resembles your data, and computes the conditional infit and outfit MSQ
of every item in each. Passed to
[`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md),
the result gives bootstrap p-values corrected across items, and an
interval per item that describes where a fitting item's statistic is
expected to fall.

## Usage

``` r
RMitemInfitCutoff(
  data,
  iterations = 400,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.95,
  dgp = c("resample", "conditional")
)
```

## Arguments

- data:

  A data.frame or matrix of item responses, items scored from 0. Only
  complete cases (rows without `NA`) are used.

- iterations:

  Integer. Number of simulated datasets. Default `400`, the recommended
  minimum for reproducible decisions. Use 1000 to 2000 for a final
  analysis (Johansson, 2026).

- parallel:

  Logical. Run the simulations in parallel via mirai when it is
  installed. Default `TRUE`.

- n_cores:

  Integer or `NULL`. Number of parallel workers. When `NULL`,
  `getOption("mc.cores")` is used. If that is unset too, a warning is
  given and the simulations run sequentially.

- verbose:

  Logical. Show a progress bar. Default `FALSE`.

- seed:

  Integer or `NULL`. Random seed for reproducibility. See
  [easyRasch2-reproducibility](https://pgmj.github.io/easyRasch2/dev/reference/easyRasch2-reproducibility.md)
  for what this guarantees and how it interacts with `parallel`.

- cutoff_method:

  Character. How the interval per item is computed: `"hdci"` (default),
  the highest-density continuous interval from
  [`ggdist::hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html)
  (needs ggdist), or `"quantile"`, the 2.5th to 97.5th percentiles.

- hdci_width:

  Numeric. Width of the interval with `cutoff_method = "hdci"`. Default
  `0.95`. Ignored with `"quantile"`.

  The interval describes where a fitting item's statistic is expected to
  fall. It is not a decision rule. Flagging every item outside a
  width-`w` interval tests all `k` items at once, with a family-wise
  error rate of `1 - w^k`, 37% for a 95% interval over nine items.
  Decisions come from the corrected p-value in
  [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md).
  Up to version 1.1.1 the default was `0.999`, which needs roughly 5000
  iterations before the interval reaches its stated width. A 95%
  interval reaches it by about 1000 (Johansson, 2026).

- dgp:

  Character. How datasets are simulated. `"resample"` (default) draws
  person locations from the WLE estimates with replacement and simulates
  responses from the model (a marginal null). `"conditional"` simulates
  each respondent's response pattern from the model given their observed
  total score, with the item parameters fixed (a conditional null). The
  conditional infit is itself conditional on the total score, so
  `"conditional"` is its natural null. **Experimental.**

## Value

A list with:

- `results`:

  data.frame with the columns `iteration`, `Item`, `InfitMSQ` and
  `OutfitMSQ`, one row per item and simulated dataset.

- `item_cutoffs`:

  data.frame with the interval per item: `Item`, `infit_low`,
  `infit_high`, `outfit_low`, `outfit_high`, computed by
  `cutoff_method`.

- `actual_iterations`:

  Number of simulated datasets that could be used. Everything downstream
  rests on this rather than on `iterations`, so it is the number to
  report.

- `requested_iterations`:

  The `iterations` argument, so that callers can tell how many simulated
  datasets were left out.

- `sample_n`:

  Number of complete cases used.

- `sample_n_total`:

  Number of respondents in `data`, before the complete-case filter.

- `sample_has_na`:

  Logical. Whether `data` contained missing values.

- `sample_summary`:

  Summary of the WLE person locations.

- `item_names`:

  Item names from `data`.

- `cutoff_method`:

  `"hdci"` or `"quantile"`.

- `hdci_width`:

  The HDCI width (used only with `"hdci"`).

- `dgp`:

  `"resample"` or `"conditional"`.

## Details

The generating model has the CML item parameters (`psychotools`) and the
WLE person locations of `data`. In each iteration a dataset is simulated
by `dgp`, the model is refitted by CML
([`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)),
and the conditional infit and outfit MSQ are computed with
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html). A
simulated dataset that cannot be refitted, for example because a
response category is empty, is left out.
[`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)
reports how many were left out in its caption. A category used only by
respondents with a zero or perfect total score stops the function before
it simulates, since CML cannot estimate its threshold.

The method and the default settings come from simulation studies of
detection (Johansson, 2025) and of iterations and multiplicity
correction (Johansson, 2026).

Parallel processing needs the optional mirai package. The iarm package
must be installed (it is in Suggests).

## References

Johansson, M. (2025). Detecting item misfit in Rasch models.
*Educational Methods & Psychometrics, 3*(18).
[doi:10.61186/emp.2025.5](https://doi.org/10.61186/emp.2025.5)

Johansson, M. (2026). Simulation-based cutoffs for conditional item fit
in Rasch models: Iterations, multiplicity correction, and decision
stability. *PsyArXiv*.
[doi:10.31234/osf.io/7pqz4_v2](https://doi.org/10.31234/osf.io/7pqz4_v2)

## See also

[`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)
for the test,
[`RMitemInfitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md)
to plot the simulated values,
[`RMitemInfitCutoffMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md)
for multiply imputed data.

## Examples

``` r
# \donttest{
if (requireNamespace("iarm", quietly = TRUE) &&
    requireNamespace("ggdist", quietly = TRUE)) {
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
  )
  colnames(sim_data) <- paste0("Item", 1:10)

  # Run 100 iterations sequentially for a quick demo
  cutoff_res <- RMitemInfitCutoff(sim_data, iterations = 100,
                                  parallel = FALSE, seed = 42)
  cutoff_res$item_cutoffs

  # Use the simulations in RMitemInfit()
  RMitemInfit(sim_data, cutoff = cutoff_res)
}
#> Bootstrap p-values are based on 100 iterations, below the recommended 400.
#> ℹ The family-wise error rate is controlled at any number of iterations. With few iterations the p-values are coarse, and which comparisons get flagged depends on the seed.
#> ℹ See Johansson (2026), doi:10.31234/osf.io/7pqz4_v2.
#> ℹ Raise `iterations` in RMitemInfitCutoff().
#> This message is displayed once per session.
#> 
#> 
#> Table: MSQ values based on conditional estimation. n = 200 respondents. Two-sided bootstrap p-values from 100 iterations, multiplicity correction: Westfall-Young step-down (FWER). Flagged on the corrected p-value at alpha = 0.05. p-values cannot be smaller than 1/(100+1) = 0.0099. This is below the recommended 400. The family-wise error rate is controlled at any number of iterations. Decisions depend on the seed, so use 1000 to 2000 for a final analysis (Johansson, 2026). Flagged: underfit (infit > 1) / overfit (infit < 1).
#> 
#> |Item   | Infit MSQ| Infit low| Infit high|      p| p (adj)|Flagged | Relative location|
#> |:------|---------:|---------:|----------:|------:|-------:|:-------|-----------------:|
#> |Item1  |     1.008|     0.912|      1.117| 0.8119|  1.0000|        |             -0.22|
#> |Item2  |     0.999|     0.928|      1.112| 0.8614|  1.0000|        |              0.14|
#> |Item3  |     0.994|     0.888|      1.109| 0.9010|  1.0000|        |              0.00|
#> |Item4  |     1.055|     0.931|      1.118| 0.2277|  0.9109|        |             -0.04|
#> |Item5  |     1.004|     0.917|      1.111| 0.9604|  1.0000|        |              0.10|
#> |Item6  |     1.032|     0.882|      1.087| 0.4752|  1.0000|        |             -0.08|
#> |Item7  |     0.943|     0.902|      1.127| 0.3564|  0.9703|        |             -0.04|
#> |Item8  |     0.988|     0.885|      1.098| 0.9703|  1.0000|        |             -0.02|
#> |Item9  |     0.933|     0.902|      1.125| 0.1485|  0.8119|        |              0.12|
#> |Item10 |     1.044|     0.910|      1.112| 0.4851|  1.0000|        |              0.38|
# }
```
