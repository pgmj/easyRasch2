# Simulation-based Q3 cutoffs

Simulates datasets from a Rasch or partial credit model that fits and
resembles your data, and computes \\Q_3\\ for every item pair in each.
Returns percentiles of the maximum minus mean \\Q_3\\ for a global
cutoff (Christensen et al., 2017) and an interval per pair. Passed to
[`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md),
the result gives bootstrap p-values corrected across pairs.

## Usage

``` r
RMlocdepQ3Cutoff(
  data,
  iterations = 400,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.95,
  estimator = c("CML", "MML"),
  dgp = c("resample", "conditional")
)
```

## Arguments

- data:

  A data.frame or matrix of item responses, items scored from 0.
  Respondents with no responses are dropped, incomplete ones are kept.

- iterations:

  Integer. Number of simulated datasets. Default `400`, the recommended
  minimum for reproducible Westfall-Young decisions and the number a 95%
  interval needs to converge (Johansson, 2026). Use 1000 to 2000 for a
  final analysis.

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

  Character. How the interval per pair in `pair_cutoffs` is computed:
  `"hdci"` (default), the highest-density continuous interval from
  [`ggdist::hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html)
  (needs ggdist), or `"quantile"`, the 2.5th to 97.5th percentiles. The
  global `suggested_cutoff` is not affected.

- hdci_width:

  Numeric in (0, 1). Width of the interval with
  `cutoff_method = "hdci"`. Default `0.95` (`0.99` before 1.2.0). The
  interval describes where a fitting pair's \\Q_3\\ is expected to fall
  and is no longer the default decision rule, so the width is chosen to
  converge at the default number of iterations. Ignored with
  `"quantile"`.

- estimator:

  Character. `"CML"` (default, CML item parameters and WLE person
  locations) or `"MML"` (`mirt`). It is stored in the result and used by
  [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
  and
  [`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md).

- dgp:

  Character. How datasets are simulated. `"resample"` (default) draws
  person locations from the WLE estimates with replacement and simulates
  responses from the model (a marginal null). Each drawn respondent
  keeps their pattern of missing responses, so the simulated \\Q_3\\ of
  a pair rests on as many respondents as the observed one.
  `"conditional"` simulates each respondent's pattern from the model
  given their observed total score and answered items, with the item
  parameters fixed (a conditional null). It needs no latent distribution
  and avoids the over-dispersion of resampled point estimates. The two
  give different cutoffs. **Experimental.**

## Value

A list with:

- `results`:

  data.frame with `mean`, `max` and `diff` (maximum minus mean \\Q_3\\),
  one row per simulated dataset.

- `pair_results`:

  data.frame with `Item1`, `Item2`, `Q3` and `iteration`, one row per
  item pair and simulated dataset.

- `pair_cutoffs`:

  data.frame with the interval per pair: `Item1`, `Item2`, `Q3_low`,
  `Q3_high`, computed by `cutoff_method`.

- `actual_iterations`:

  Number of simulated datasets that could be used.

- `requested_iterations`:

  The `iterations` argument.

- `sample_n`:

  Number of respondents used (those with at least one response).

- `sample_n_total`:

  Number of respondents in `data`.

- `sample_has_na`:

  Logical. Whether `data` contained missing values.

- `sample_summary`:

  Summary of the WLE person locations.

- `item_names`:

  Item names.

- `max_diff`, `sd_diff`:

  Maximum and SD of `diff`.

- `p95`, `p99`, `p995`, `p999`:

  Percentiles of `diff`.

- `suggested_cutoff`:

  The 99th percentile (`p99`), the global cutoff for
  [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md).

- `cutoff_method`:

  `"hdci"` or `"quantile"`.

- `hdci_width`:

  The HDCI width (used only with `"hdci"`).

- `estimator`:

  `"CML"` or `"MML"`.

- `dgp`:

  `"resample"` or `"conditional"`.

## Details

The generating model is fitted once, with CML item parameters
(`psychotools`) and WLE person locations. In each iteration a dataset is
simulated by `dgp`, the model is refitted, and \\Q_3\\ is computed with
`estimator`. A simulated dataset that cannot be refitted is left out,
and
[`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
reports how many in its caption.

Parallel processing needs the optional mirai package.

## References

Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values
for Yen's \\Q_3\\: Identification of local dependence in the Rasch
model. *Applied Psychological Measurement, 41*(3), 178–194.
[doi:10.1177/0146621616677520](https://doi.org/10.1177/0146621616677520)

Johansson, M. (2026). Simulation-based cutoffs for conditional item fit
in Rasch models: Iterations, multiplicity correction, and decision
stability. *PsyArXiv*.
[doi:10.31234/osf.io/7pqz4_v2](https://doi.org/10.31234/osf.io/7pqz4_v2)

## See also

[`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
for the tables,
[`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md)
for the plots.

## Examples

``` r
# \donttest{
if (requireNamespace("ggdist", quietly = TRUE)) {
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
  )
  colnames(sim_data) <- paste0("Item", 1:10)

  # Few iterations for a fast example. Use 1000 or more in real analyses.
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
