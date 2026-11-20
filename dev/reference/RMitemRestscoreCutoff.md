# Simulation-based item-restscore cutoffs

Simulates datasets from a Rasch or partial credit model that fits and
resembles your data, and computes each item's observed minus expected
item-restscore gamma (Kreiner, 2011) in each. Passed to
[`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md),
the result gives bootstrap p-values corrected across items, and an
interval per item that describes where a fitting item's difference is
expected to fall.

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

  A data.frame or matrix of item responses, items scored from 0. Only
  complete cases (rows without `NA`) are used.

- iterations:

  Integer. Number of simulated datasets. Default `400`. Use 1000 or more
  for a final analysis (see Reproducibility in
  [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)).

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

  The interval describes where a fitting item's difference is expected
  to fall. It is not a decision rule. Flagging every item outside a
  width-`w` interval tests all `k` items at once, with a family-wise
  error rate of `1 - w^k`. Decisions come from the corrected p-value in
  [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md),
  as shown for infit in Johansson (2026).

- dgp:

  Character. How datasets are simulated. `"conditional"` (default)
  simulates each respondent's response pattern from the model given
  their observed total score, with the item parameters fixed (a
  conditional null). The expected gamma is computed from the observed
  score distribution, which the conditional null keeps. `"resample"`
  draws person locations from the WLE estimates with replacement and
  simulates responses from the model (a marginal null). With the Rasch
  model true, pooled over four designs, the family-wise rate with
  `correction = "fwer"` was 5.0% with the conditional null and 6.4% with
  the resample null. The conditional null takes about 1.2 to 1.6 times
  as long. **Experimental.**

## Value

A list with:

- `results`:

  data.frame with the columns `iteration`, `Item`, `Observed`,
  `Expected` and `Difference`, one row per item and simulated dataset.
  `Difference` is `Observed - Expected`, the statistic
  [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
  tests.

- `item_cutoffs`:

  data.frame with the interval per item: `Item`, `diff_low`,
  `diff_high`, computed by `cutoff_method`.

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

  `"conditional"` or `"resample"`.

## Details

**Why simulate.** The asymptotic test in
[`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html)
divides the observed minus expected gamma by the standard error of the
observed gamma alone. The expected gamma is estimated from the same data
and correlates with the observed one, so that standard error is too
large for the difference, and the difference is also biased upwards in
small samples. With the Rasch model true, the test is liberal for
dichotomous items in small or mistargeted samples and conservative for
polytomous items, and it flags too few underfitting items in every case.
The simulated reference distribution replaces the asymptotic one and
absorbs both problems. The simulation study of the asymptotic test is
reported at
<https://pgmj.github.io/rasch_book/chapters/restscore-asymptotic-null.html>.

**Simulation.** The generating model has the CML item parameters
(`psychotools`) and the WLE person locations of `data`. In each
iteration a dataset is simulated by `dgp`, the model is refitted by CML
([`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)),
and the observed and expected gamma are computed as in
[`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html),
by a faster internal routine that skips the standard errors and the
rounding of the printed values. The refit matters. The expected gamma
varies from sample to sample because the thresholds do, and holding them
fixed would reproduce the problem the simulation is meant to solve. A
simulated dataset that cannot be refitted, for example because a
response category is empty, is left out, and
[`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
reports how many in its caption.

Parallel processing needs the optional mirai package. The iarm package
must be installed (it is in Suggests).

## References

Johansson, M. (2026). Simulation-based cutoffs for conditional item fit
in Rasch models: Iterations, multiplicity correction, and decision
stability. *PsyArXiv*.
[doi:10.31234/osf.io/7pqz4_v2](https://doi.org/10.31234/osf.io/7pqz4_v2)

Kreiner, S. (2011). A Note on Item-Restscore Association in Rasch
Models. *Applied Psychological Measurement, 35*(7), 557-561.
[doi:10.1177/0146621611410227](https://doi.org/10.1177/0146621611410227)

## See also

[`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
for the test,
[`RMitemRestscorePlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscorePlot.md)
to plot the simulated values.

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
#> Bootstrap p-values are based on 100 iterations, below the recommended 400.
#> ℹ The family-wise error rate is controlled at any number of iterations. With few iterations the p-values are coarse, and which comparisons get flagged depends on the seed.
#> ℹ See Johansson (2026), doi:10.31234/osf.io/7pqz4_v2.
#> ℹ Raise `iterations` in RMitemRestscoreCutoff().
#> This message is displayed once per session.
#> 
#> 
#> 
#> Table: Item-restscore associations. n = 200 respondents. Equal-tailed parametric bootstrap p-values for the difference from 100 iterations, conditional DGP, multiplicity correction: Westfall-Young step-down (FWER). p-values cannot be smaller than 2/(100+1) = 0.02. This is below the recommended 400, where decisions depend on the seed (Johansson, 2026). Flagged (adj. p < 0.05): overfit = observed above expected (over-discrimination, often local dependence); underfit = below (under-discrimination, often multidimensionality/noise).
#> 
#> |Item   | Observed| Expected| Difference| Diff low| Diff high|      p| p (adj)|Flagged | Rel. location|
#> |:------|--------:|--------:|----------:|--------:|---------:|------:|-------:|:-------|-------------:|
#> |Item1  |     0.02|     0.03|     -0.006|   -0.174|     0.176| 0.9901|  1.0000|        |         -0.22|
#> |Item2  |     0.04|     0.03|      0.008|   -0.155|     0.174| 0.9703|  1.0000|        |          0.14|
#> |Item3  |     0.05|     0.03|      0.013|   -0.195|     0.151| 0.9307|  1.0000|        |          0.00|
#> |Item4  |    -0.08|     0.03|     -0.110|   -0.135|     0.183| 0.1980|  1.0000|        |         -0.04|
#> |Item5  |     0.02|     0.03|     -0.012|   -0.175|     0.167| 0.9307|  1.0000|        |          0.10|
#> |Item6  |    -0.05|     0.03|     -0.083|   -0.172|     0.142| 0.2376|  1.0000|        |         -0.08|
#> |Item7  |     0.15|     0.03|      0.122|   -0.211|     0.186| 0.1584|  1.0000|        |         -0.04|
#> |Item8  |     0.08|     0.03|      0.052|   -0.149|     0.157| 0.6535|  1.0000|        |         -0.02|
#> |Item9  |     0.17|     0.03|      0.139|   -0.151|     0.169| 0.1188|  0.7921|        |          0.12|
#> |Item10 |    -0.06|     0.03|     -0.095|   -0.186|     0.170| 0.2376|  1.0000|        |          0.38|
# }
```
