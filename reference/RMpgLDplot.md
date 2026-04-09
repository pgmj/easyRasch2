# Plot Distribution of Simulated Partial Gamma LD Values

Visualises the distribution of simulation-based partial gamma LD values
from
[`RMpgLDcutoff`](https://pgmj.github.io/easyRasch2/reference/RMpgLDcutoff.md),
optionally overlaying observed partial gamma values computed from real
data via [`partgam_LD`](https://rdrr.io/pkg/iarm/man/partgam_LD.html).

## Usage

``` r
RMpgLDplot(simfit, data, items = NULL)
```

## Arguments

- simfit:

  The return value of
  [`RMpgLDcutoff`](https://pgmj.github.io/easyRasch2/reference/RMpgLDcutoff.md)
  (a list with components `results`, `pair_cutoffs`,
  `actual_iterations`, `sample_n`, and `item_names`).

- data:

  Optional. A data.frame or matrix of item responses for computing and
  overlaying observed partial gamma values. Items must be scored
  starting at 0 (non-negative integers). When provided, the plot
  includes orange diamond markers for the observed partial gamma
  alongside the simulated distribution, plus segment summaries from the
  cutoff intervals.

- items:

  Optional character vector of item names to include in the plot. Only
  item pairs where **both** items are in this vector will be shown. When
  `NULL` (default), all item pairs are plotted.

## Value

A `ggplot` object.

## Details

Uses
[`ggdist::stat_dotsinterval()`](https://mjskay.github.io/ggdist/reference/stat_dotsinterval.html)
(when `data` is not supplied) or
[`ggdist::stat_dots()`](https://mjskay.github.io/ggdist/reference/stat_dots.html)
(when `data` is supplied) with `point_interval = "median_hdci"` and
`.width = c(0.66, 0.95, 0.99)`.

The plot shows one row per item pair (labelled as "Item1 - Item2"). Only
direction 1 (rest score = total - Item2) is plotted, matching the
convention used in the simulation.

When `data` is **not** supplied, the function plots the simulated
partial gamma distributions as dot-interval plots using
[`ggdist::stat_dotsinterval()`](https://mjskay.github.io/ggdist/reference/stat_dotsinterval.html)
with median and Highest Density Continuous Interval (HDCI) summaries.

When `data` **is** supplied, the function:

1.  Computes observed partial gamma values via
    [`iarm::partgam_LD()`](https://rdrr.io/pkg/iarm/man/partgam_LD.html).

2.  Overlays observed gamma values as orange diamond markers on the
    simulated distributions.

3.  Shows per-pair cutoff intervals (from `simfit$pair_cutoffs`) as
    black line segments, with thicker segments for the 66\\ black dots
    for the median.

The `ggplot2`, `ggdist`, and optionally `iarm` packages must be
installed (they are in Suggests, not Imports).

## See also

[`RMpgLDcutoff`](https://pgmj.github.io/easyRasch2/reference/RMpgLDcutoff.md),
[`RMpartgamLD`](https://pgmj.github.io/easyRasch2/reference/RMpartgamLD.md)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)

# Run simulation
cutoff_res <- RMpgLDcutoff(sim_data, iterations = 100, parallel = FALSE,
                           seed = 42)

# Simulated distribution only
RMpgLDplot(cutoff_res)

# With observed partial gamma overlaid
RMpgLDplot(cutoff_res, data = sim_data)

# Plot only a subset of items
RMpgLDplot(cutoff_res, data = sim_data,
           items = c("Item1", "Item2", "Item3"))
} # }
```
