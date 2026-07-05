# Conditional Item Characteristic Curves

Faceted panel of Conditional Item Characteristic Curves (cICCs), one per
item. Each panel shows the model-expected conditional item score
\\E\[X_i \mid R = r\]\\ against the total score \\r\\, together with the
average observed item score within class intervals of the total score
and their confidence intervals. When a `dif_var` is supplied, observed
averages are computed separately per group, turning each panel into a
visual DIF check, and the partial-gamma DIF magnitude is reported per
item.

## Usage

``` r
RMitemICCPlot(
  data,
  dif_var = NULL,
  method = c("cut", "score"),
  class_intervals = 4,
  ci = TRUE,
  error_band = FALSE,
  conf_level = 0.95,
  min_n = 8L,
  items = NULL,
  output = c("patchwork", "list")
)
```

## Arguments

- data:

  A data.frame or matrix of item responses (non-negative integers,
  0-based). One column per item, one row per person.

- dif_var:

  Optional vector of length `nrow(data)` (factor or coercible to factor)
  defining a DIF grouping variable. When supplied, observed averages are
  shown separately per group and the partial-gamma DIF coefficient
  ([`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html),
  as in
  [`RMdifGamma()`](https://pgmj.github.io/easyRasch2/reference/RMdifGamma.md))
  is annotated per panel. Default `NULL` (no DIF).

- method:

  One of `"cut"` (default) or `"score"`. `"cut"` groups respondents into
  `class_intervals` total-score bins (equal-count quantile bins, common
  across groups); `"score"` uses every observed total score as its own
  group.

- class_intervals:

  Integer \>= 2. Number of class intervals when `method = "cut"`.
  Default `4`. Ignored when `method = "score"`.

- ci:

  Logical. Draw confidence intervals (error bars) on the observed
  class-interval means. Default `TRUE`.

- error_band:

  Logical. Add a shaded band for the model's interval of the *observed
  mean* at each total score, \\E \pm z \sqrt{\mathrm{Var}/n_r}\\ (`n_r`
  = number of respondents at that total score), as in RASCHplot's
  `CICCplot`. Default `FALSE`.

- conf_level:

  Numeric in (0, 1). Confidence level for the observed error bars and
  the model band. Default `0.95`.

- min_n:

  Integer. A group-by-interval cell needs at least this many respondents
  to contribute an observed point + CI; sparser cells are skipped.
  Default `8` (the package's per-cell stability floor).

- items:

  Optional character or integer vector selecting which items to plot.
  The model is always fitted on all items; only the rendering is
  filtered. Default `NULL` (all items).

- output:

  One of `"patchwork"` (default) – composite patchwork figure – or
  `"list"` – a named list of per-item `ggplot` objects.

## Value

Either a `patchwork`/`ggplot` composite (default) or a named list of
per-item `ggplot` objects (`output = "list"`). In DIF mode the per-item
partial-gamma table is attached as `attr(., "dif_gamma")`.

## Details

The model curve and its conditional variance come from the shared CML
engine (CML item thresholds via `psychotools` and the exact conditional
distribution of each item score given the total score, via elementary
symmetric functions). The conditional-ICC approach follows Buchardt,
Christensen & Jensen (2023) and their `RASCHplot` package; the
implementation here is native to easyRasch2's CML/WLE engine.

**Conditioning.** The expected curve is the exact conditional
expectation of the item score given the total score (it accounts for the
item being part of the total), not a marginal ICC.

**Class intervals.** With `method = "cut"`, total scores are split into
`class_intervals` equal-count bins using common boundaries (so groups
share the x-axis); each bin contributes one observed point at its mean
total score. With `method = "score"`, every observed total score is a
point. If the score distribution is too sparse to form `class_intervals`
distinct bins, the function falls back to score-level points.

**Confidence intervals.** Observed error bars use the normal
approximation \\\bar{x}\_l \pm z \sqrt{\mathrm{var}(x_l) / n_l}\\ within
each (group, interval) cell, clamped to the item's score range; cells
with fewer than `min_n` respondents are dropped. In DIF mode this makes
sparse group differences visibly uncertain rather than over-interpreted.

**DIF magnitude.** The annotated partial gamma (Bjorner et al., 1998) is
the association between the item score and the group conditional on the
rest score – a non-parametric effect size, complementary to the
total-score-conditional visual. It is the same statistic as
[`RMdifGamma()`](https://pgmj.github.io/easyRasch2/reference/RMdifGamma.md).

## References

Buchardt, A.-S., Christensen, K. B., & Jensen, S. N. (2023). Visualizing
Rasch item fit using conditional item characteristic curves in R.
*Psychological Test and Assessment Modeling, 65*(2), 206-219.

Andersen, E. B. (1995). Polytomous Rasch models and their estimation. In
G. H. Fischer & I. W. Molenaar (Eds.), *Rasch Models: Foundations,
Recent Developments, and Applications* (pp. 271-291). Springer.
(Conditional expectation of item scores given the total score; formula
15.22.)

Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., & Bech, P.
(1998). Differential item functioning in the Danish translation of the
SF-36. *Journal of Clinical Epidemiology, 51*(11), 1189-1202.

## See also

[`RMdifGamma`](https://pgmj.github.io/easyRasch2/reference/RMdifGamma.md),
[`RMitemRestscore`](https://pgmj.github.io/easyRasch2/reference/RMitemrestscore.md)
