# Simulation-Based Partial Gamma LD Cutoff Determination

Uses parametric bootstrap simulation to determine appropriate cutoff
values for partial gamma Local Dependence analysis via
[`partgam_LD`](https://rdrr.io/pkg/iarm/man/partgam_LD.html). Under a
correctly fitting Rasch model where items are locally independent, this
function generates the expected distribution of partial gamma values per
item pair, providing empirical critical values.

## Usage

``` r
RMpgLDcutoff(
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

1.  Resamples person parameters (thetas) with replacement from ML
    estimates.

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

Supports both **dichotomous** data (via
[`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html) and
[`psychotools::rrm()`](https://rdrr.io/pkg/psychotools/man/rrm.html))
and **polytomous** data (via
[`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html) and an internal
partial credit score simulator).

Parallel processing is provided by the `mirai` package (optional).
Install it with `install.packages("mirai")` to enable parallelisation.

The `iarm` package must be installed (it is in Suggests, not Imports).

## References

Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) *Rasch Models in
Health*. Iste and Wiley (2013), pp. 133–135.

## See also

[`partgam_LD`](https://rdrr.io/pkg/iarm/man/partgam_LD.html),
[`RMpartgamLD`](https://pgmj.github.io/easyRasch2/reference/RMpartgamLD.md),
[`RMpgLDplot`](https://pgmj.github.io/easyRasch2/reference/RMpgLDplot.md)

## Examples

``` r
# \donttest{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)

# Run 100 iterations sequentially for a quick demo
cutoff_res <- RMpgLDcutoff(sim_data, iterations = 100, parallel = FALSE,
                           seed = 42)
cutoff_res$pair_cutoffs
#>    Item1  Item2  gamma_low gamma_high
#> 1  Item1  Item2 -0.4306569  0.4306761
#> 2  Item1  Item3 -0.3582424  0.3249651
#> 3  Item1  Item4 -0.5086744  0.4304419
#> 4  Item1  Item5 -0.3945250  0.3072829
#> 5  Item1  Item6 -0.4643478  0.5241633
#> 6  Item1  Item7 -0.4729345  0.4150514
#> 7  Item1  Item8 -0.5020921  0.3392174
#> 8  Item1  Item9 -0.5251263  0.3824885
#> 9  Item1 Item10 -0.3755858  0.3699346
#> 10 Item2  Item3 -0.4204398  0.3994601
#> 11 Item2  Item4 -0.3644315  0.3534209
#> 12 Item2  Item5 -0.4473970  0.4878331
#> 13 Item2  Item6 -0.3730728  0.3733529
#> 14 Item2  Item7 -0.4377953  0.4086538
#> 15 Item2  Item8 -0.4017972  0.4588756
#> 16 Item2  Item9 -0.3388430  0.3944811
#> 17 Item2 Item10 -0.3537332  0.3087994
#> 18 Item3  Item4 -0.3931570  0.4112150
#> 19 Item3  Item5 -0.4140481  0.3298013
#> 20 Item3  Item6 -0.3815603  0.4300039
#> 21 Item3  Item7 -0.4352428  0.3839662
#> 22 Item3  Item8 -0.4105651  0.3771429
#> 23 Item3  Item9 -0.3846154  0.3843663
#> 24 Item3 Item10 -0.3751171  0.4253835
#> 25 Item4  Item5 -0.3895447  0.3688490
#> 26 Item4  Item6 -0.4051667  0.4102190
#> 27 Item4  Item7 -0.4657110  0.3573265
#> 28 Item4  Item8 -0.4961380  0.4708393
#> 29 Item4  Item9 -0.4414125  0.4671026
#> 30 Item4 Item10 -0.3295775  0.2959772
#> 31 Item5  Item6 -0.3552068  0.4885720
#> 32 Item5  Item7 -0.4122038  0.3920705
#> 33 Item5  Item8 -0.4357542  0.4475191
#> 34 Item5  Item9 -0.3253012  0.3421405
#> 35 Item5 Item10 -0.4073171  0.3440970
#> 36 Item6  Item7 -0.3982301  0.4649695
#> 37 Item6  Item8 -0.3762811  0.3933511
#> 38 Item6  Item9 -0.3669915  0.4166667
#> 39 Item6 Item10 -0.4095238  0.4507380
#> 40 Item7  Item8 -0.3514916  0.4971444
#> 41 Item7  Item9 -0.3357228  0.4801166
#> 42 Item7 Item10 -0.3532645  0.3842173
#> 43 Item8  Item9 -0.4445880  0.4311649
#> 44 Item8 Item10 -0.3395225  0.4713264
#> 45 Item9 Item10 -0.3807474  0.4108911
# }
```
