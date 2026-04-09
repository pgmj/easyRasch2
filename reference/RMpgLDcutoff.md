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
if (FALSE) { # \dontrun{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)

# Run 100 iterations sequentially for a quick demo
cutoff_res <- RMpgLDcutoff(sim_data, iterations = 100, parallel = FALSE,
                           seed = 42)
cutoff_res$pair_cutoffs
} # }
```
