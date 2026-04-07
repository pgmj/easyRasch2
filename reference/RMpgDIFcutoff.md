# Simulation-Based Partial Gamma DIF Cutoff Determination

Uses parametric bootstrap simulation to determine appropriate cutoff
values for partial gamma DIF analysis via
[`partgam_DIF`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html). Under a
correctly fitting Rasch model where the DIF variable is unrelated to
item responses (i.e., no true DIF), this function generates the expected
distribution of absolute partial gamma values per item, providing
empirical critical values.

## Usage

``` r
RMpgDIFcutoff(
  data,
  dif_var,
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

- dif_var:

  A vector (factor, character, or integer) defining group membership for
  DIF analysis. Must have the same length as `nrow(data)`. The actual
  group labels are used to determine the number of groups and their
  relative sizes; during simulation, respondents are randomly assigned
  to groups with the same proportions, so there is no true DIF by
  construction.

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

  data.frame with columns `iteration`, `Item`, and `gamma` (one row per
  item per successful iteration).

- `item_cutoffs`:

  data.frame with per-item cutoff summaries: `Item`, `gamma_low`,
  `gamma_high`. Bounds are computed using the method specified by
  `cutoff_method`.

- `actual_iterations`:

  Number of successful iterations.

- `sample_n`:

  Number of complete cases used.

- `sample_summary`:

  Summary statistics of estimated person parameters.

- `item_names`:

  Character vector of item names from data.

- `dif_group_sizes`:

  Named integer vector of group sizes used in the simulation (matches
  proportions in the observed `dif_var`).

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

3.  Creates a random DIF variable by sampling group labels with the same
    proportions as the observed `dif_var`, so there is **no true DIF**
    by construction.

4.  Computes partial gamma DIF statistics via
    [`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html).

The distribution of partial gamma values across iterations provides
empirical critical values per item. Values from real data that fall
outside these bounds suggest DIF that exceeds what would be expected by
chance under a correctly fitting Rasch model. Failed iterations (e.g.,
due to convergence issues or degenerate data) are silently discarded.

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

Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., & Bech, P.
(1998). Differential item functioning in the Danish translation of the
SF-36. *Journal of Clinical Epidemiology*, 51(11), 1189–1202.

Henninger, M., Radek, J., Sengewald, M.-A., & Strobl, C. (2024). Partial
credit trees meet the partial gamma coefficient for quantifying DIF and
DSF in polytomous items. OSF Preprints.
[doi:10.31234/osf.io/47sah](https://doi.org/10.31234/osf.io/47sah)

## See also

[`partgam_DIF`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)
dif_sex <- sample(c("male", "female"), 200, replace = TRUE)

# Run 100 iterations sequentially for a quick demo
cutoff_res <- RMpgDIFcutoff(sim_data, dif_var = dif_sex,
                                  iterations = 100, parallel = FALSE,
                                  seed = 42)
cutoff_res$item_cutoffs
} # }
```
