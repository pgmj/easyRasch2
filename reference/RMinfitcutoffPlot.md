# Plot Distribution of Simulated Infit and Outfit MSQ Values

Visualises the distribution of simulation-based conditional item fit
values from
[`RMinfitcutoff`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md),
optionally overlaying observed item fit from the original data.

## Usage

``` r
RMinfitcutoffPlot(simfit, data, output = "infit")
```

## Arguments

- simfit:

  The return value of
  [`RMinfitcutoff`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md)
  (a list with components `results`, `item_cutoffs`,
  `actual_iterations`, `sample_n`, and `item_names`).

- data:

  Optional. A data.frame or matrix of item responses for computing and
  overlaying observed conditional item fit values. Items must be scored
  starting at 0 (non-negative integers). When provided, the plot
  includes orange diamond markers for the observed infit/outfit MSQ
  alongside the simulated distribution, plus segment summaries from the
  cutoff intervals.

- output:

  Character string. Either `"infit"` (default) to show only the infit
  panel, `"outfit"` to show only the outfit panel, or `"both"` to show
  infit and outfit side by side (requires the `patchwork` package when
  `data` is supplied).

## Value

A `ggplot` object (or a `patchwork` object when `output = "both"` and
`data` is supplied).

## Details

Uses
[`ggdist::stat_dotsinterval()`](https://mjskay.github.io/ggdist/reference/stat_dotsinterval.html)
(when `data` is not supplied) or
[`ggdist::stat_dots()`](https://mjskay.github.io/ggdist/reference/stat_dots.html)
(when `data` is supplied) with `point_interval = "median_hdci"` and
`.width = c(0.66, 0.999)`.

When `data` is **not** supplied, the function plots the simulated MSQ
distributions as dot-interval plots using
[`ggdist::stat_dotsinterval()`](https://mjskay.github.io/ggdist/reference/stat_dotsinterval.html)
with median and Highest Density Continuous Interval (HDCI) summaries,
faceted by statistic (InfitMSQ / OutfitMSQ).

When `data` **is** supplied, the function:

1.  Fits a Rasch model
    ([`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html) for dichotomous
    data or [`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html) for
    polytomous data) and computes observed conditional infit and outfit
    MSQ via
    [`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html).

2.  Overlays observed fit values as orange diamond markers on the
    simulated distributions.

3.  Shows per-item cutoff intervals (from `simfit$item_cutoffs`) as
    black line segments, with thicker segments for the 66% HDCI range
    and black dots for the median.

The `ggplot2`, `ggdist`, and optionally `patchwork` packages must be
installed (they are in Suggests, not Imports).

## See also

[`RMinfitcutoff`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md),
[`RMiteminfit`](https://pgmj.github.io/easyRasch2/reference/RMiteminfit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)

# Run simulation
cutoff_res <- RMinfitcutoff(sim_data, iterations = 100, parallel = FALSE,
                            seed = 42)

# Simulated distribution only (infit + outfit faceted)
RMinfitcutoffPlot(cutoff_res)

# With observed fit overlaid (infit only, the default)
RMinfitcutoffPlot(cutoff_res, data = sim_data)

# Both infit and outfit panels side by side
RMinfitcutoffPlot(cutoff_res, data = sim_data, output = "both")
} # }
```
