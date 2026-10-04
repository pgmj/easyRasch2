# Plot the Simulated Item-Restscore Null Distribution

Visualises the per-item null distribution of the observed minus expected
item-restscore gamma from
[`RMitemRestscoreCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md),
optionally overlaying the observed differences from the original data.

## Usage

``` r
RMitemRestscorePlot(simfit, data)
```

## Arguments

- simfit:

  The return value of
  [`RMitemRestscoreCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md)
  (a list with components `results`, `item_cutoffs`,
  `actual_iterations`, `sample_n`, and `item_names`).

- data:

  Optional. A data.frame or matrix of item responses for computing and
  overlaying the observed item-restscore differences. Items must be
  scored starting at 0 (non-negative integers). When provided, the plot
  includes orange diamond markers for the observed difference alongside
  the simulated distribution, plus segment summaries of the intervals.

## Value

A `ggplot` object.

## Details

Uses
[`ggdist::stat_dotsinterval()`](https://mjskay.github.io/ggdist/reference/stat_dotsinterval.html)
(when `data` is not supplied) or
[`ggdist::stat_dots()`](https://mjskay.github.io/ggdist/reference/stat_dots.html)
(when `data` is supplied) with `point_interval = "median_hdci"`. The
outer `.width` follows `simfit$hdci_width`, so the shaded interval
matches the one
[`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
tabulates.

The x-axis is the difference between observed and model-expected gamma,
the statistic
[`RMitemRestscore`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
tests. A dashed line marks zero. The simulated distributions are
generally not centred on zero in small samples, which is one of the
reasons the asymptotic test is miscalibrated (see
[`RMitemRestscoreCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md)).
Positive values indicate over-discrimination (overfit), negative values
under-discrimination (underfit).

When `data` **is** supplied, the observed differences are computed with
[`RMitemRestscore`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
and overlaid as orange diamonds, with per-item intervals drawn as black
line segments (thicker for the 66% range) and black dots for the
simulated median.

The `ggplot2` and `ggdist` packages must be installed (they are in
Suggests, not Imports).

## See also

[`RMitemRestscoreCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md),
[`RMitemRestscore`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)

## Examples

``` r
# \donttest{
if (requireNamespace("iarm", quietly = TRUE) &&
    requireNamespace("ggdist", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
  )
  colnames(sim_data) <- paste0("Item", 1:10)

  cutoff_res <- RMitemRestscoreCutoff(sim_data, iterations = 100,
                                      parallel = FALSE, seed = 42)

  # Simulated distribution only
  RMitemRestscorePlot(cutoff_res)

  # With the observed differences overlaid
  RMitemRestscorePlot(cutoff_res, data = sim_data)
}
#> 

# }
```
