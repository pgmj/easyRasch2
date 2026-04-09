# Plot Distribution of Simulated Partial Gamma DIF Values

Visualises the distribution of simulation-based partial gamma DIF values
from
[`RMpgDIFcutoff`](https://pgmj.github.io/easyRasch2/reference/RMpgDIFcutoff.md),
optionally overlaying observed partial gamma values computed from real
data via [`partgam_DIF`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html).

## Usage

``` r
RMpgDIFplot(simfit, data, dif_var)
```

## Arguments

- simfit:

  The return value of
  [`RMpgDIFcutoff`](https://pgmj.github.io/easyRasch2/reference/RMpgDIFcutoff.md)
  (a list with components `results`, `item_cutoffs`,
  `actual_iterations`, `sample_n`, and `item_names`).

- data:

  Optional. A data.frame or matrix of item responses for computing and
  overlaying observed partial gamma values. Items must be scored
  starting at 0 (non-negative integers). When provided, the plot
  includes orange diamond markers for the observed partial gamma
  alongside the simulated distribution, plus segment summaries from the
  cutoff intervals.

- dif_var:

  Required when `data` is supplied. A vector (factor, character, or
  integer) defining group membership for the DIF analysis. Must have the
  same length as `nrow(data)`.

## Value

A `ggplot` object.

## Details

Uses
[`ggdist::stat_dotsinterval()`](https://mjskay.github.io/ggdist/reference/stat_dotsinterval.html)
(when `data` is not supplied) or
[`ggdist::stat_dots()`](https://mjskay.github.io/ggdist/reference/stat_dots.html)
(when `data` is supplied) with `point_interval = "median_hdci"` and
`.width = c(0.66, 0.95, 0.99)`.

When `data` is **not** supplied, the function plots the simulated
partial gamma distributions as dot-interval plots using
[`ggdist::stat_dotsinterval()`](https://mjskay.github.io/ggdist/reference/stat_dotsinterval.html)
with median and Highest Density Continuous Interval (HDCI) summaries.

When `data` **is** supplied (along with `dif_var`), the function:

1.  Computes observed partial gamma values via
    [`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html).

2.  Overlays observed gamma values as orange diamond markers on the
    simulated distributions.

3.  Shows per-item cutoff intervals (from `simfit$item_cutoffs`) as
    black line segments, with thicker segments for the 66\\ black dots
    for the median.

The `ggplot2`, `ggdist`, and optionally `iarm` packages must be
installed (they are in Suggests, not Imports).

## See also

[`RMpgDIFcutoff`](https://pgmj.github.io/easyRasch2/reference/RMpgDIFcutoff.md),
[`RMpartgamDIF`](https://pgmj.github.io/easyRasch2/reference/RMpartgamDIF.md)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)
dif_group <- factor(sample(c("A", "B"), 200, replace = TRUE))

# Run simulation
cutoff_res <- RMpgDIFcutoff(sim_data, dif_var = dif_group,
                             iterations = 100, parallel = FALSE, seed = 42)

# Simulated distribution only
RMpgDIFplot(cutoff_res)

# With observed partial gamma overlaid
RMpgDIFplot(cutoff_res, data = sim_data, dif_var = dif_group)
} # }
```
