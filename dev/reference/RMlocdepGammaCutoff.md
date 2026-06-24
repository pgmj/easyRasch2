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
  iterations = 250,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.99
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
  `0.99` (99\\ `cutoff_method = "quantile"`.

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

3.  Computes partial gamma LD statistics via
    [`iarm::partgam_LD()`](https://rdrr.io/pkg/iarm/man/partgam_LD.html).

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
[`RMlocdepGamma`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md),
[`RMlocdepGammaPlot`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaPlot.md)

## Examples

``` r
# \donttest{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)

# Run 100 iterations sequentially for a quick demo
cutoff_res <- RMlocdepGammaCutoff(sim_data, iterations = 100, parallel = FALSE,
                           seed = 42)
cutoff_res$pair_cutoffs
#>    Item1  Item2  gamma_low gamma_high
#> 1  Item1  Item2 -0.4246238  0.3198758
#> 2  Item1  Item3 -0.4725291  0.4400000
#> 3  Item1  Item4 -0.4455579  0.3262195
#> 4  Item1  Item5 -0.4192593  0.4463083
#> 5  Item1  Item6 -0.3890363  0.3372093
#> 6  Item1  Item7 -0.3239512  0.4533333
#> 7  Item1  Item8 -0.4260870  0.4747872
#> 8  Item1  Item9 -0.3201970  0.3996248
#> 9  Item1 Item10 -0.3210702  0.3171954
#> 10 Item2  Item3 -0.3423763  0.3693694
#> 11 Item2  Item4 -0.2822086  0.4506687
#> 12 Item2  Item5 -0.3933887  0.4179104
#> 13 Item2  Item6 -0.3816425  0.3362769
#> 14 Item2  Item7 -0.5614692  0.3977456
#> 15 Item2  Item8 -0.4305085  0.3895771
#> 16 Item2  Item9 -0.4109589  0.3747739
#> 17 Item2 Item10 -0.3531353  0.3608790
#> 18 Item3  Item4 -0.3698551  0.3174603
#> 19 Item3  Item5 -0.4548193  0.4344904
#> 20 Item3  Item6 -0.4166667  0.4460015
#> 21 Item3  Item7 -0.4143335  0.4424779
#> 22 Item3  Item8 -0.4346727  0.4184874
#> 23 Item3  Item9 -0.3793103  0.3470952
#> 24 Item3 Item10 -0.3508501  0.4581142
#> 25 Item4  Item5 -0.3603667  0.4174067
#> 26 Item4  Item6 -0.3335321  0.4556213
#> 27 Item4  Item7 -0.3969336  0.3440059
#> 28 Item4  Item8 -0.3689320  0.3484576
#> 29 Item4  Item9 -0.3971429  0.3337701
#> 30 Item4 Item10 -0.3641851  0.3455481
#> 31 Item5  Item6 -0.3874426  0.3315698
#> 32 Item5  Item7 -0.4620253  0.4496595
#> 33 Item5  Item8 -0.3712256  0.4384670
#> 34 Item5  Item9 -0.3704415  0.4069529
#> 35 Item5 Item10 -0.4049501  0.3196481
#> 36 Item6  Item7 -0.4504065  0.3748283
#> 37 Item6  Item8 -0.2495922  0.3668176
#> 38 Item6  Item9 -0.3822401  0.3897638
#> 39 Item6 Item10 -0.3070326  0.4840983
#> 40 Item7  Item8 -0.4427245  0.3660965
#> 41 Item7  Item9 -0.3649123  0.3639609
#> 42 Item7 Item10 -0.3774802  0.2980132
#> 43 Item8  Item9 -0.4102012  0.3160813
#> 44 Item8 Item10 -0.4171975  0.3995774
#> 45 Item9 Item10 -0.3581972  0.3491311
# }
```
