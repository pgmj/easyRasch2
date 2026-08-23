# Simulation-Based Partial Gamma LD Cutoff Determination

Uses parametric bootstrap simulation to determine appropriate cutoff
values for partial gamma Local Dependence analysis via
[`partgam_LD`](https://rdrr.io/pkg/iarm/man/partgam_LD.html). Under a
correctly fitting Rasch model where items are locally independent, this
function generates the expected distribution of partial gamma values per
item pair, providing empirical critical values.

## Usage

``` r
RMlocdepGammaCutoff(
  data,
  iterations = 400,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.95
)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Only complete cases (rows
  without any `NA`) are used.

- iterations:

  Integer. Number of simulation iterations (default 400, was 250 before
  1.2.0). 400 is the calibrated floor for the Westfall-Young correction
  (Johansson, 2026) and the count a 95\\ converge. Use 1000 to 2000 for
  a final analysis.

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

  Character string specifying how cutoff intervals are computed. Either
  `"hdci"` (default) for the Highest Density Interval via
  [`ggdist::hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html),
  or `"quantile"` for the 2.5th/97.5th percentiles via
  [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html).

- hdci_width:

  Numeric. Width of the HDCI when `cutoff_method = "hdci"`. Default is
  `0.95` (95\\ describes where a fitting pair's coefficient is expected
  to fall and is no longer the default decision rule, so the width is
  chosen to converge at the default iteration count rather than to imply
  an error rate. Ignored when `cutoff_method = "quantile"`.

## Value

A list with components:

- `results`:

  data.frame with columns `iteration`, `Item1`, `Item2`, and `gamma`
  (one row per item pair per successful iteration). Contains results
  from direction 1 only (rest score = total - Item2), which is the
  conventional direction.

- `pair_cutoffs`:

  data.frame with per-pair cutoff summaries: `Item1`, `Item2`,
  `gamma_low`, `gamma_high`. Bounds are computed using the method
  specified by `cutoff_method`.

- `actual_iterations`:

  Number of successful iterations.

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

## Details

For each simulation iteration the function:

1.  Resamples person parameters (thetas) with replacement from the WLE
    person locations.

2.  Simulates item response data under a Rasch model (dichotomous via
    [`psychotools::rrm()`](https://rdrr.io/pkg/psychotools/man/rrm.html)
    or polytomous via an internal partial credit simulator).

3.  Computes partial gamma for every item pair in the canonical
    rest-score direction. The coefficients are identical to those of
    [`iarm::partgam_LD()`](https://rdrr.io/pkg/iarm/man/partgam_LD.html),
    but are computed by a vectorised internal, since `iarm` also derives
    the asymptotic standard error and confidence interval that a
    simulated null does not need and costs roughly two orders of
    magnitude more per iteration.

Because the data are simulated under the Rasch model, items are locally
independent by construction. The distribution of partial gamma values
across iterations provides empirical critical values per item pair.
Values from real data that fall outside these bounds suggest local
dependence that exceeds what would be expected by chance. Failed
iterations (e.g., due to convergence issues or degenerate data) are
silently discarded.

The generating model uses CML item thresholds via
[`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)
(a dichotomous item is a 2-category PCM) and WLE person locations,
consistent with the rest of the package; responses are simulated with
[`psychotools::rrm()`](https://rdrr.io/pkg/psychotools/man/rrm.html)
(dichotomous) or an internal partial credit score simulator
(polytomous).

Parallel processing is provided by the `mirai` package (optional).
Install it with `install.packages("mirai")` to enable parallelisation.

The `iarm` package must be installed (it is in Suggests, not Imports).

## References

Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) (2013). *Rasch
Models in Health*, pp. 133–135. ISTE & Wiley.
[doi:10.1002/9781118574454](https://doi.org/10.1002/9781118574454)

## See also

[`partgam_LD`](https://rdrr.io/pkg/iarm/man/partgam_LD.html),
[`RMlocdepGamma`](https://pgmj.github.io/easyRasch2/reference/RMlocdepGamma.md),
[`RMlocdepGammaPlot`](https://pgmj.github.io/easyRasch2/reference/RMlocdepGammaPlot.md)

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
  cutoff_res <- RMlocdepGammaCutoff(sim_data, iterations = 100,
                                    parallel = FALSE, seed = 42)
  cutoff_res$pair_cutoffs
}
#>    Item1  Item2  gamma_low gamma_high
#> 1  Item1  Item2 -0.2777321  0.3198758
#> 2  Item1  Item3 -0.2571429  0.3843537
#> 3  Item1  Item4 -0.2265372  0.3262195
#> 4  Item1  Item5 -0.3109244  0.3537519
#> 5  Item1  Item6 -0.2680723  0.3297003
#> 6  Item1  Item7 -0.2089041  0.2947559
#> 7  Item1  Item8 -0.3167702  0.3052264
#> 8  Item1  Item9 -0.2711268  0.3800000
#> 9  Item1 Item10 -0.2714777  0.2960000
#> 10 Item2  Item3 -0.2605364  0.3693694
#> 11 Item2  Item4 -0.2141653  0.3376906
#> 12 Item2  Item5 -0.1967865  0.4179104
#> 13 Item2  Item6 -0.3203540  0.2933912
#> 14 Item2  Item7 -0.3580264  0.4016620
#> 15 Item2  Item8 -0.2684564  0.3580705
#> 16 Item2  Item9 -0.3421927  0.2688498
#> 17 Item2 Item10 -0.2885662  0.2979177
#> 18 Item3  Item4 -0.2500000  0.3216561
#> 19 Item3  Item5 -0.2911392  0.2585670
#> 20 Item3  Item6 -0.3357143  0.3543860
#> 21 Item3  Item7 -0.3416537  0.3750000
#> 22 Item3  Item8 -0.3453355  0.3333333
#> 23 Item3  Item9 -0.2324723  0.3308271
#> 24 Item3 Item10 -0.2238806  0.4440559
#> 25 Item4  Item5 -0.2923077  0.3902848
#> 26 Item4  Item6 -0.2087912  0.4170854
#> 27 Item4  Item7 -0.3563636  0.2890365
#> 28 Item4  Item8 -0.2226402  0.3563579
#> 29 Item4  Item9 -0.3491311  0.2613241
#> 30 Item4 Item10 -0.2740741  0.3381555
#> 31 Item5  Item6 -0.2629969  0.3025335
#> 32 Item5  Item7 -0.2527472  0.3256151
#> 33 Item5  Item8 -0.3303965  0.3283358
#> 34 Item5  Item9 -0.2671233  0.3651877
#> 35 Item5 Item10 -0.2531876  0.3196481
#> 36 Item6  Item7 -0.2467190  0.3983051
#> 37 Item6  Item8 -0.2130584  0.3502235
#> 38 Item6  Item9 -0.2437886  0.3966102
#> 39 Item6 Item10 -0.2382609  0.3498350
#> 40 Item7  Item8 -0.3460803  0.2953216
#> 41 Item7  Item9 -0.2693498  0.3362319
#> 42 Item7 Item10 -0.2040201  0.3078261
#> 43 Item8  Item9 -0.2593918  0.3160813
#> 44 Item8 Item10 -0.2667877  0.3750000
#> 45 Item9 Item10 -0.2635379  0.3019197
# }
```
