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

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Only complete cases (rows
  without any `NA`) are used.

- iterations:

  Integer. Number of simulation iterations (default 400).

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

  Integer or `NULL`. Random seed for reproducibility. See
  [easyRasch2-reproducibility](https://pgmj.github.io/easyRasch2/dev/reference/easyRasch2-reproducibility.md)
  for what this guarantees and how it interacts with `parallel`.

- cutoff_method:

  Character string specifying how cutoff intervals are computed. Either
  `"hdci"` (default) for the Highest Density Interval via
  [`ggdist::hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html),
  or `"quantile"` for the 2.5th/97.5th percentiles via
  [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html).

- hdci_width:

  Numeric. Width of the HDCI when `cutoff_method = "hdci"`. Default is
  `0.95` (95% HDCI). Ignored when `cutoff_method = "quantile"`.

  The interval is a **description** of where a fitting item's statistic
  is expected to fall, not a decision rule. Flagging every item outside
  a width-`w` interval tests all `k` items at once, so the family-wise
  error rate is `1 - w^k`, which is 37% for a 95% interval over nine
  items. Decisions should come from the corrected p-value instead
  ([`RMitemInfit`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)
  with `p_value = NULL` and the full object returned here). The default
  was `0.999` up to and including version 1.1.1, which needs roughly
  5000 iterations before the interval reaches its stated width; `0.95`
  reaches it by about 1000, so the band shown to readers means close to
  what it says (Johansson, 2026).

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

  Number of successful iterations. Everything downstream rests on this
  rather than on `iterations`, so it is the number to report.

- `requested_iterations`:

  The `iterations` argument, kept so callers can tell how many simulated
  datasets were discarded.

- `sample_n`:

  Number of complete cases used.

- `sample_n_total`:

  Number of respondents in the raw input data, before the complete-case
  filter.

- `sample_has_na`:

  Logical. Whether the raw input data contained any missing values.

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

## References

Johansson, M. (2025). Detecting item misfit in Rasch models.
*Educational Methods & Psychometrics, 3*(18).
[doi:10.61186/emp.2025.5](https://doi.org/10.61186/emp.2025.5)

Johansson, M. (2026). Simulation-based cutoffs for conditional item fit
in Rasch models: Iterations, multiplicity correction, and decision
stability. *PsyArXiv*.
[doi:10.31234/osf.io/7pqz4_v2](https://doi.org/10.31234/osf.io/7pqz4_v2)

## See also

[`RMitemInfit`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)

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

  # Use the cutoffs in RMitemInfit()
  RMitemInfit(sim_data)
}
#> 
#> 
#> Table: MSQ values based on conditional estimation. n = 200 respondents.
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
