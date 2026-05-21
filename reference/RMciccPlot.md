# Conditional Item Characteristic Curves

Faceted panel of Conditional Item Characteristic Curves (cICCs) for all
items in `data`. Each panel shows the model-expected ICC together with
average observed item scores within class intervals (or for every
possible total score, with `method = "score"`). When a `dif_var` is
supplied, observed averages are computed separately per group, so each
panel becomes a visual DIF check.

## Usage

``` r
RMciccPlot(
  data,
  class_intervals = 5L,
  method = c("cut", "score"),
  dif_var = NULL,
  output = c("patchwork", "list")
)
```

## Arguments

- data:

  A data.frame or matrix of item responses (non-negative integers,
  0-based). One column per item, one row per person.

- class_intervals:

  Integer. Number of class intervals used to group respondents when
  `method = "cut"`. Default `5`. Each class interval needs enough
  respondents to compute a stable average observed item score; with
  small samples or many possible total scores, lower this. Ignored when
  `method = "score"`.

- method:

  One of `"cut"` (default) or `"score"`. `"cut"` groups respondents into
  `class_intervals` adjacent bins; `"score"` uses every possible total
  raw score.

- dif_var:

  Optional vector of length `nrow(data)` (factor or coercible to factor)
  defining a DIF grouping variable. When supplied, observed averages are
  shown separately per group. Default `NULL` (no DIF).

- output:

  One of `"patchwork"` (default) – composite patchwork figure – or
  `"list"` – a named list of per-item `ggplot` objects.

## Value

Either a `ggplot` / `patchwork` composite (default), or a list of
`ggplot` objects (`output = "list"`), one per item.

## Details

Wraps [`iarm::ICCplot()`](https://rdrr.io/pkg/iarm/man/ICCplot.html)
(Santiago, Mueller, & Müller, 2022) for the per-item computation and
composes the panels with
[`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html).

**Class intervals vs every total score.** With `method = "cut"` (the
default) respondents are divided into `class_intervals` bins of roughly
equal size based on their total raw score; each bin's average observed
item score is plotted as a point. With `method = "score"`, the average
observed item score is computed for every observed total raw score. The
cut method is the conventional Rasch-analysis default and works well at
moderate sample sizes; the score method is more granular but can be
unstable when raw-score totals are sparsely populated.

**DIF mode.** With `dif_var` supplied,
[`iarm::ICCplot()`](https://rdrr.io/pkg/iarm/man/ICCplot.html) computes
average observed item scores separately for each level of the grouping
variable. Each item's panel then shows the model-expected ICC (one black
line) plus one set of observed points per DIF group. Stark visual
divergence between groups suggests DIF.

**Dependencies.** Requires `iarm` for the per-item curve, and
`patchwork` + `ggplot2` for the figure composition. All three are in
Suggests with runtime guards.

## References

Mueller, M., & Santiago, P. H. R. (2022). *iarm: Item Analysis in Rasch
Models*. R package version 0.4.3.
<https://CRAN.R-project.org/package=iarm>

## See also

[`iarm::ICCplot()`](https://rdrr.io/pkg/iarm/man/ICCplot.html),
[`RMtargeting()`](https://pgmj.github.io/easyRasch2/reference/RMtargeting.md),
[`RMtileplot()`](https://pgmj.github.io/easyRasch2/reference/RMtileplot.md)

## Examples

``` r
# \donttest{
data("pcmdat2", package = "eRm")
RMciccPlot(pcmdat2, class_intervals = 5)


# With DIF grouping
set.seed(1)
grp <- factor(sample(c("A", "B"), nrow(pcmdat2), replace = TRUE))
RMciccPlot(pcmdat2, dif_var = grp)
#> Error: iarm::ICCplot() failed on item I1: There are not enough subjects in each total score to produce this number of class intervals. This often means `class_intervals` (4) is too large for the available total scores. Try a smaller value or use `method = "score"`.
# }
```
