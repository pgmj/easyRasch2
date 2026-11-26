# Category functioning given the total score

Four diagnostic plots of how each item's response categories are used,
all built on the probability of each category given the respondent's
total score, \\P(X_i = j \mid r)\\. Under the partial credit model this
probability depends on the item parameters alone, so no person locations
are estimated and the comparisons are not affected by their measurement
error. Items are shown in rows.
[`RMitemCatFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFit.md)
gives the test that the calibration plot shows.

## Usage

``` r
RMitemCatFitPlot(
  data,
  type = c("residuals", "adjacent", "restscore", "calibration"),
  fit = NULL,
  dif_var = NULL,
  grouping = c("quantile", "width", "score", "manual"),
  class_intervals = 5L,
  score_breaks = NULL,
  iterations = 200L,
  conf_level = 0.95,
  item_labels = NULL,
  label_wrap = 15L,
  parallel = FALSE,
  n_cores = NULL,
  seed = NULL,
  output = c("ggplot", "dataframe")
)
```

## Arguments

- data:

  A data.frame or matrix of item responses, items scored from 0.
  Respondents with missing responses are left out (complete cases),
  since the total score must cover the same items for everyone.

- type:

  Character. Which plot to draw. Default `"residuals"`, or
  `"calibration"` when `fit` is given, since that plot shows the test in
  [`RMitemCatFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFit.md).

  - `"residuals"`: standardised residuals of the number of respondents
    in each category, within groups of total scores, against simulated
    ranges when `fit` is given.

  - `"adjacent"`: for each pair of adjacent categories, the share of
    respondents in the upper one, against the total score.

  - `"restscore"`: the mean and SD of the rest score among respondents
    in each category, observed minus expected, against simulated
    intervals.

  - `"calibration"`: PAV-adjusted calibration of the cumulative events
    \\X_i \ge k\\ against the model's predicted probability.

- fit:

  Optional result of
  [`RMitemCatFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFit.md)
  with `output = "list"`, made from the same `data` (checked). Its
  simulations are used for the rest-score and calibration plots, so
  `iterations` and `seed` are ignored, and give the ranges in the
  residual plot. Each item's corrected p from
  [`RMitemCatFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFit.md)
  is shown under its label, in bold when below .05 (marked with an
  asterisk when ggtext is not installed), and the fit's item labels are
  used unless `item_labels` is given. Not available with `dif_var`.

- dif_var:

  Optional vector with one value per row of `data` (factor, character or
  numeric, treated as groups), at most four groups. Each group's curves
  are then drawn against predictions from one fit to all respondents,
  with a band showing the range expected without differential item
  functioning (DIF). Available for `type = "calibration"`, `"adjacent"`
  and `"residuals"`. Rows with a missing group are left out with a
  message. See "Groups" in Details.

- grouping:

  How respondents are grouped by total score for `type = "residuals"`,
  as in
  [`RMitemICCPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemICCPlot.md):
  `"quantile"` (default), `class_intervals` groups with about equal
  numbers of respondents, `"width"`, `class_intervals` intervals of
  equal width over the observed total scores, `"score"`, every observed
  total score its own group, or `"manual"`, groups set by
  `score_breaks`. Groups without respondents are left out. The caption
  names the grouping used.

- class_intervals:

  Integer, at least 2. Number of groups for `grouping = "quantile"` and
  `"width"`. Default `5L`.

- score_breaks:

  Integer vector for `grouping = "manual"`: the total scores at which a
  new group starts, in increasing order. Default `NULL`.

- iterations:

  Integer, at least 20. Number of simulated datasets for
  `type = "restscore"` and `type = "calibration"` when `fit` is not
  given. Default `200L`. The ends of each pointwise range are extreme
  quantiles of the simulated values, so below about 1000 datasets they,
  and which points fall outside them, change with the seed (the caption
  says so). With `dif_var`, the number of permutations of the group
  labels behind the bands, default `1000L`.

- conf_level:

  Numeric in (0, 1). Coverage of the pointwise bands and intervals, and
  with `dif_var` of the band for each curve. Default `0.95`.

- item_labels:

  Optional character vector of unique item labels for the row strips,
  one per column of `data`. Default: column names.

- label_wrap:

  Integer. Characters per line in the row strips. Default `15L`.

- parallel, n_cores:

  Logical / integer. Run the simulations in parallel via mirai when
  available. Default sequential. When `n_cores` is `NULL`,
  `getOption("mc.cores")` is used. If that is unset too, a warning is
  given and the simulations run sequentially. Not used with `dif_var`,
  where the permutations run sequentially.

- seed:

  Optional integer for reproducible simulations. See
  [easyRasch2-reproducibility](https://pgmj.github.io/easyRasch2/dev/reference/easyRasch2-reproducibility.md).

- output:

  Character. `"ggplot"` (default) or `"dataframe"` for the data behind
  the plot.

## Value

- `output = "ggplot"`: a
  [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object, one row of panels per item.

- `output = "dataframe"`: the data behind the plot, with columns `item`
  and, by `type`:

  - `"residuals"`: `score_group`, `category`, `observed`, `expected`,
    `z`, `lower`, `upper`, `iterations_used` (the last three `NA`
    without `fit`).

  - `"adjacent"`: `pair`, `score`, `n_pair`, `observed`, `expected`,
    `lower`, `upper`.

  - `"restscore"`: `statistic`, `category`, `observed`, `lower`,
    `upper`, `iterations_used`.

  - `"calibration"`: `event`, `predicted`, `cep`, `n`, `lower`, `upper`.

  With `dif_var`, each type also has `group` (the group with its n),
  `perm_mean` (the mean over the permutations), `lower` and `upper` (the
  band), `outside` (the point is outside its band) and `iterations_used`
  (the number of permutations). The non-DIF reference columns are
  replaced: `adjacent` keeps `expected` (the model's share),
  `calibration` and `adjacent` have no binomial or simulated range, and
  `residuals` has no cells with zero expected variance.

## Details

Respondents with the lowest or highest possible total score are left
out, since their responses are fixed by the score. The item parameters
are estimated by conditional maximum likelihood (`psychotools`). As in
[`RMitemCatFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFit.md),
a category that is empty, or used only by respondents with an extreme
total score, stops the function with an error that names the item.

A plain comparison of observed and expected category counts would show
nothing. The category counts are sufficient statistics for the item
parameters, so at the CML estimates the expected counts equal the
observed ones. Differences show only once respondents are split by total
score, which every plot here does.

**Residuals.** For each item, category and score group, the observed
number of respondents minus the expected number, divided by its standard
deviation under the model. This comparison of observed and expected
counts within score groups is the one behind Glas's (1988) R1c test. The
standard deviation treats the item parameters as known, but they are
estimated from the same counts, so the residuals vary less than a
standard normal: with 5 score groups about 2.5% of cells fall beyond
\\\pm 1.96\\ when the model holds, fewer with fewer groups. In cells
with few expected respondents the residual can only be large in the
positive direction. Fixed lines at \\\pm 1.96\\ would therefore mislead,
and none are drawn. With `fit`, each cell gets the pointwise range of
the residual in the datasets
[`RMitemCatFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFit.md)
simulated (regenerated for the chosen grouping, without refitting),
which accounts for both. The plot has many cells, so some fall outside
their range even when the model holds.

**Adjacent categories.** Among respondents in category \\j - 1\\ or
\\j\\ at total score \\r\\, the share in \\j\\. Under the model this
share rises with \\r\\. Observed shares that do not rise suggest that
respondents do not distinguish the two categories. The band is the
pointwise binomial range around the expected share, given the number of
respondents in the pair.

**Rest score.** The mean and SD of the rest score (total score minus the
item's score) among respondents in each category, observed minus
expected. A category used by respondents with too low or too high rest
scores, or with too wide a spread, stands out. Intervals come from
datasets simulated from the model given each respondent's total score,
with the item parameters re-estimated in each, so the uncertainty in the
estimates is included. An item that discriminates less than the model
expects (underfit) gives a rest-score SD above its expectation in every
category, and compressed means: low categories above and high categories
below their expectation. An item that discriminates more (overfit, often
from local dependence) gives the reverse. The SD panel shows this as
clearly as the mean panel, and in the middle category only the SD
changes. This is a property of the item, not of its categories (see
"What a flag means" in
[`RMitemCatFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFit.md)).
A raised SD in a single category points to that category instead, for
example a middle category chosen by respondents at all levels.

**Calibration.** For each cumulative event \\X_i \ge k\\, the model
predicts \\p = P(X_i \ge k \mid r)\\. The observed events are regressed
on \\p\\ by pool-adjacent-violators (PAV), which gives the conditional
event probability, the best non-decreasing fit to the observed events
(Dimitriadis et al., 2021). This is the PAV-adjusted reliability diagram
recommended for ordinal outcomes by Säilynoja et al. (2025), with the
prediction taken from the total score instead of a person estimate.
Agreement puts the curve on the diagonal. The band is the pointwise
range of the same curve over datasets simulated given the total scores
at the estimated item parameters. A curve that is too flat means the
event separates respondents less well than the model predicts, and too
steep means better. Because CML reproduces each category's count, the
curve cannot lie above or below the diagonal overall, only cross it.

Without `dif_var`, all bands and intervals are pointwise. They show
where to look, not whether an item fits, and with many points some fall
outside even when the model holds.

**Groups (`dif_var`).** The item parameters are estimated once, from all
respondents used, and every group's curves are drawn against the same
predictions. Under the Rasch model without DIF, an item's response is
independent of the group given the total score (Kreiner & Christensen,
2011), whatever the groups' ability distributions. Curves that separate
between groups therefore point to DIF, while curves that depart from the
model together point to misfit of the item for everyone.

The band is a simultaneous consistency band (Dimitriadis et al., 2021)
for each curve with coverage `conf_level` (default 95%): the range the
whole curve keeps in that share of permutations of the group labels
within total scores, a permutation that keeps the number of respondents
of each group at every score (an exact conditional test, Kreiner, 1987).
For the residual plot, a "curve" is all cells of one item and group. The
band is centred on the pooled data, not on the model, so an item that
misfits for everyone moves all groups' bands together. Each band is cut
to the values the permutations reached, which matters for cells the
model expects to be almost empty. In the residual plot, the y-axis
covers the bars and at least -4 to 4, and a band that reaches further is
cut at the panel edge. Points (bars for residuals) outside their band
are marked. The bands are not corrected across curves, so with many
curves some cross their band by chance. In simulations without DIF, at
most about 5% of curves crossed their band (dev study in the package
repository, `dev/catfit_dif/`).

Two conditions limit the reading. First, the independence given the
total score holds under the Rasch model, not when items differ in
discrimination (Kreiner & Christensen, 2011). Then, if the groups differ
in ability, curves can separate without DIF, so check item fit first
([`RMitemCatFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFit.md),
[`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)).
Second, the total score includes the other items: DIF in several items
in the same direction shifts the matching, and DIF-free items can then
look different between groups.
[`RMdifAnchor()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifAnchor.md)
(item and threshold locations) and
[`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md)
test DIF.

## References

Dimitriadis, T., Gneiting, T., & Jordan, A. I. (2021). Stable
reliability diagrams for probabilistic classifiers. *Proceedings of the
National Academy of Sciences, 118*(8), e2016191118.
[doi:10.1073/pnas.2016191118](https://doi.org/10.1073/pnas.2016191118)

Glas, C. A. W. (1988). The derivation of some tests for the Rasch model
from the multinomial distribution. *Psychometrika, 53*(4), 525–546.
[doi:10.1007/BF02294405](https://doi.org/10.1007/BF02294405)

Kreiner, S. (1987). Analysis of multidimensional contingency tables by
exact conditional tests: Techniques and strategies. *Scandinavian
Journal of Statistics, 14*(2), 97–112.

Kreiner, S., & Christensen, K. B. (2011). Item screening in graphical
loglinear Rasch models. *Psychometrika, 76*(2), 228–256.
[doi:10.1007/s11336-011-9203-y](https://doi.org/10.1007/s11336-011-9203-y)

Säilynoja, T., Johnson, A. R., Martin, O. A., & Vehtari, A. (2025).
Recommendations for visual predictive checks in Bayesian workflow.
*arXiv*.
[doi:10.48550/arXiv.2503.01509](https://doi.org/10.48550/arXiv.2503.01509)

## See also

[`RMitemCatFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFit.md)
for the test,
[`RMitemCatProb()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatProb.md)
for the model-implied category curves,
[`RMitemICCPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemICCPlot.md)
for conditional item characteristic curves.

## Examples

``` r
# \donttest{
if (requireNamespace("eRm", quietly = TRUE)) {
  data(pcmdat2, package = "eRm")
  RMitemCatFitPlot(pcmdat2)
  RMitemCatFitPlot(pcmdat2, type = "adjacent")
  RMitemCatFitPlot(pcmdat2, type = "calibration", iterations = 100,
                     seed = 1)
  # By group (here a random split, so no DIF is expected)
  set.seed(1)
  group <- sample(c("A", "B"), nrow(pcmdat2), replace = TRUE)
  RMitemCatFitPlot(pcmdat2, type = "adjacent", dif_var = group, seed = 1)
}

# }
```
