# Simulation-Based Item-Restscore Null Distribution

Uses parametric bootstrap simulation to build the null distribution of
the item-restscore statistic for
[`RMitemRestscore`](https://pgmj.github.io/easyRasch2/reference/RMitemrestscore.md).
This function simulates data from a correctly fitting Rasch model that
mimics your data and returns, per item, the simulated difference between
observed and expected item-restscore gamma.

## Usage

``` r
RMitemRestscoreCutoff(
  data,
  iterations = 400,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.95,
  dgp = c("conditional", "resample")
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
  [easyRasch2-reproducibility](https://pgmj.github.io/easyRasch2/reference/easyRasch2-reproducibility.md)
  for what this guarantees and how it interacts with `parallel`.

- cutoff_method:

  Character string specifying how the intervals are computed. Either
  `"hdci"` (default) for the Highest Density Interval via
  [`ggdist::hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html),
  or `"quantile"` for the 2.5th/97.5th percentiles via
  [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html).

- hdci_width:

  Numeric. Width of the HDCI when `cutoff_method = "hdci"`. Default is
  `0.95` (95% HDCI). Ignored when `cutoff_method = "quantile"`.

  The interval is a **description** of where a fitting item's difference
  is expected to fall, not a decision rule. Flagging every item outside
  a width-`w` interval tests all `k` items at once, so the family-wise
  error rate is `1 - w^k`. Decisions should come from the corrected
  p-value instead
  ([`RMitemRestscore`](https://pgmj.github.io/easyRasch2/reference/RMitemrestscore.md)
  with `p_value = NULL` and the full object returned here).

- dgp:

  Character. Data-generating process for the parametric bootstrap.
  `"conditional"` (default) simulates each respondent's pattern from the
  exact Rasch conditional distribution given their observed total score,
  item parameters fixed (a *conditional* null). The expected gamma is
  computed from the observed score distribution, which the conditional
  null holds fixed. `"resample"` resamples WLE person locations with
  replacement and simulates responses under the model (a *marginal*
  null). In simulation under a true Rasch model, the conditional null
  gave a family-wise rate of 5.0 percent with `correction = "fwer"` and
  the resample null 6.4 percent, pooled over four designs. The
  conditional null takes about 1.2 to 1.6 times as long.
  **Experimental.**

## Value

A list with components:

- `results`:

  data.frame with columns `iteration`, `Item`, `Observed`, `Expected`,
  `Difference` (one row per item per successful iteration). `Difference`
  is `Observed - Expected`, the statistic
  [`RMitemRestscore`](https://pgmj.github.io/easyRasch2/reference/RMitemrestscore.md)
  tests.

- `item_cutoffs`:

  data.frame with per-item interval bounds for the difference: `Item`,
  `diff_low`, `diff_high`. Bounds are computed using the method
  specified by `cutoff_method`.

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

  The method used to compute the intervals (`"hdci"` or `"quantile"`).

- `hdci_width`:

  The HDCI width used (only meaningful when `cutoff_method = "hdci"`).

- `dgp`:

  The data-generating process used (`"resample"` or `"conditional"`).

## Details

The asymptotic test in
[`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html)
divides the observed minus expected gamma by the standard error of the
observed gamma alone. The expected gamma is estimated from the same data
and correlates with the observed one, so that standard error is too
large for the difference, and the difference is also biased upwards in
small samples. Under a true Rasch model the resulting test is liberal
for dichotomous items in small or mistargeted samples and conservative
for polytomous items, in every case flagging too few underfitting items.
The bootstrap null replaces the asymptotic reference distribution and
absorbs both problems.

The generating model is CML item parameters (via `psychotools`) with WLE
person locations. For each iteration a dataset is simulated under the
chosen `dgp`, the model is **refitted** by CML
([`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)),
and the observed and expected item-restscore gamma are computed as in
[`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html),
by a faster internal routine that skips the standard errors and the
rounding of the printed values. The refit matters: the expected gamma
varies from sample to sample because the thresholds do, and holding them
fixed would reproduce the problem the bootstrap exists to solve. Failed
iterations (e.g., degenerate simulated data) are silently discarded.

Parallel processing is provided by the `mirai` package (optional).
Install it with `install.packages("mirai")` to enable parallelisation.

The `iarm` package must be installed (it is in Suggests, not Imports).

## References

Kreiner, S. (2011). A Note on Item-Restscore Association in Rasch
Models. *Applied Psychological Measurement, 35*(7), 557-561.
[doi:10.1177/0146621611410227](https://doi.org/10.1177/0146621611410227)

Johansson, M. (2026). Simulation-based cutoffs for conditional item fit
in Rasch models: Iterations, multiplicity correction, and decision
stability. *PsyArXiv*.
[doi:10.31234/osf.io/7pqz4_v2](https://doi.org/10.31234/osf.io/7pqz4_v2)

## See also

[`RMitemRestscore`](https://pgmj.github.io/easyRasch2/reference/RMitemrestscore.md),
[`RMitemRestscorePlot`](https://pgmj.github.io/easyRasch2/reference/RMitemRestscorePlot.md)

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
  cutoff_res <- RMitemRestscoreCutoff(sim_data, iterations = 100,
                                      parallel = FALSE, seed = 42)
  cutoff_res$item_cutoffs

  # Flag on bootstrap p-values in RMitemRestscore()
  RMitemRestscore(sim_data, cutoff = cutoff_res)
}
#> Bootstrap p-values are based on 100 iterations, below the calibrated floor of 400.
#> ℹ Below 400 the Westfall-Young correction is mildly liberal under the null, so the family-wise error rate is above the nominal level.
#> ℹ See Johansson (2026), doi:10.31234/osf.io/7pqz4_v2.
#> ℹ Raise `iterations` in RMitemRestscoreCutoff().
#> This message is displayed once per session.
#> 
#> 
#> 
#> Table: Item-restscore associations. n = 200 respondents. Two-sided parametric bootstrap p-values for the difference from 100 iterations, conditional DGP, multiplicity correction: Westfall-Young step-down (FWER). p-values cannot be smaller than 1/(100+1) = 0.0099. This is below the calibrated floor of 400, where the correction is mildly liberal (Johansson, 2026). Flagged (adj. p < 0.05): overfit = observed above expected (over-discrimination, often local dependence); underfit = below (under-discrimination, often multidimensionality/noise).
#> 
#> |Item   | Observed| Expected| Difference| Diff low| Diff high|      p| p (adj)|Flagged | Rel. location|
#> |:------|--------:|--------:|----------:|--------:|---------:|------:|-------:|:-------|-------------:|
#> |Item1  |     0.02|     0.03|     -0.006|   -0.174|     0.176| 1.0000|  1.0000|        |         -0.22|
#> |Item2  |     0.04|     0.03|      0.008|   -0.155|     0.174| 0.9703|  1.0000|        |          0.14|
#> |Item3  |     0.05|     0.03|      0.013|   -0.195|     0.151| 0.9010|  1.0000|        |          0.00|
#> |Item4  |    -0.08|     0.03|     -0.110|   -0.135|     0.183| 0.2574|  0.8614|        |         -0.04|
#> |Item5  |     0.02|     0.03|     -0.012|   -0.175|     0.167| 0.9703|  1.0000|        |          0.10|
#> |Item6  |    -0.05|     0.03|     -0.083|   -0.172|     0.142| 0.1980|  0.8614|        |         -0.08|
#> |Item7  |     0.15|     0.03|      0.122|   -0.211|     0.186| 0.1683|  0.8317|        |         -0.04|
#> |Item8  |     0.08|     0.03|      0.052|   -0.149|     0.157| 0.6238|  1.0000|        |         -0.02|
#> |Item9  |     0.17|     0.03|      0.139|   -0.151|     0.169| 0.0990|  0.5743|        |          0.12|
#> |Item10 |    -0.06|     0.03|     -0.095|   -0.186|     0.170| 0.2574|  0.8614|        |          0.38|
# }
```
