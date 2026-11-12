# DIF analysis via Andersen's likelihood-ratio test

Fits a Rasch model by conditional maximum likelihood in the full sample
and in each group of an external grouping variable, tests the difference
with Andersen's likelihood-ratio test, and reports per-group item
locations (or per-group threshold locations) together with their
standard errors. `level` (item or threshold) and `output` (data.frame,
kable or ggplot) give six views of the same fits.

## Usage

``` r
RMdifLR(
  data,
  dif_var,
  model = c("auto", "PCM", "RM"),
  level = c("item", "threshold"),
  output = c("kable", "dataframe", "ggplot"),
  cutoff = 0.5,
  conf = 0.95,
  sort = FALSE
)
```

## Arguments

- data:

  A data.frame or matrix of item responses (non-negative integers,
  0-based). One column per item, one row per person. Person IDs and
  grouping variables must not be included – pass the grouping variable
  separately via `dif_var`.

- dif_var:

  Vector of length `nrow(data)` (factor, character, or numeric) defining
  the DIF grouping variable. Coerced to factor; unused levels are
  dropped. Rows where `dif_var` is `NA` are dropped with a message. Must
  result in at least 2 groups after cleaning.

- model:

  One of `"auto"` (default), `"PCM"`, or `"RM"`. `"auto"` fits the Rasch
  model when the data are dichotomous (max response = 1) and the partial
  credit model otherwise (the two coincide for dichotomous data). `"RM"`
  errors on polytomous data.

- level:

  One of `"item"` (default) or `"threshold"`. `"item"` reports each
  item's mean threshold location per group; `"threshold"` reports each
  individual threshold per group. For dichotomous (RM) data the two
  views are equivalent (one threshold per item).

- output:

  One of `"kable"` (default), `"dataframe"`, or `"ggplot"`. The
  data.frame has a logical `Flagged` column (unless `cutoff = NULL`),
  the kable bolds flagged rows, and the ggplot shows confidence
  intervals.

- cutoff:

  Numeric or `NULL`. Threshold (in logits) for the `Flagged` column: a
  row is flagged when `MaxDiff > cutoff`, where `MaxDiff` is the
  difference between the largest and smallest per-group location for
  that item (or threshold). Set to `NULL` to suppress flagging. Default
  `0.5`. This default is a convention in applied Rasch work (Jin et al.,
  2013; Kim et al., 2017; Pickard et al., 2006), not a value with a
  statistical or substantive justification. See Details.

- conf:

  Numeric in (0, 1). Confidence level used for the ggplot error bars.
  Default `0.95`.

- sort:

  Logical. Sort the rows of the data.frame and kable output by `MaxDiff`
  (descending). Default `FALSE`.

## Value

A data.frame, a `knitr_kable` object, or a `ggplot` object, depending on
`output`.

The data.frame has one row per item (`level = "item"`) or per item x
threshold (`level = "threshold"`), with columns `Item` (and `Threshold`
at threshold level), one numeric column per group level, an `All` column
for the unsplit fit, `MaxDiff`, `Flagged` (when `cutoff` is non-`NULL`),
and matching `SE_*` columns.

The Andersen LR test result is attached as `attr(result, "lr_test")` on
the data.frame, in the kable footnote, and in the ggplot caption (LR
\\\chi^2\\, df, p-value).

## Details

The Partial Credit Model (PCM) is fitted by default for polytomous data
and the dichotomous Rasch Model (RM) is fitted when all responses are
0/1; this can be overridden via `model`.

Each fit is a CML partial credit model from
[`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)
(a dichotomous item is a 2-category PCM). Within each fit the thresholds
sum to zero, as in
[`RMitemParameters`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemParameters.md),
so each group's locations are centred on that group's own mean
threshold. The `All` column is the fit to the full dataset. An item's
location is the mean of its thresholds, and its SE is the SE of that
mean, computed from the thresholds' covariance matrix.

The Andersen LR statistic is twice the difference between the sum of the
group conditional log-likelihoods and the full-sample one, with df equal
to the extra free parameters of the group fits. It equals
[`eRm::LRtest()`](https://rdrr.io/pkg/eRm/man/LRtest.html) on the same
data. Respondents with missing responses are kept (CML uses the items
they answered).

**Response categories.** Every category an item uses must be used in
every group, or that item's thresholds cannot be compared across groups.
The same holds when a category is used in a group only by respondents
with a zero or perfect score there, since CML leaves them out and the
threshold cannot be estimated. The function then stops and names the
items, categories and groups
([`eRm::LRtest()`](https://rdrr.io/pkg/eRm/man/LRtest.html) instead left
items with an unused category out of the test without saying so, and
returned diverging estimates in the second case). With small groups this
is common for sparse categories: merge them with an adjacent one,
combine groups, or leave the item out.

**DIF size and the cutoff.** The `Flagged` column applies one fixed
cutoff to `MaxDiff`, and no single value suits every scale. An
appropriate cutoff depends on the number of items, the number of
response categories and the type of DIF (Van de Winckel et al., 2022).
Treat a flag as a prompt for a closer look, not as a verdict. Judging
DIF means weighing its size against its practical consequences
(Henninger et al., 2025). These are best evaluated by simulation:
estimate person locations with and without accounting for DIF, for
instance by splitting the affected item by group, and check whether the
difference changes substantive conclusions (Johansson et al., 2026). For
polytomous items, DIF can differ between an item's thresholds, which
`level = "threshold"` shows (Henninger et al., 2025). `MaxDiff` is a
point estimate, and with small groups its uncertainty is large (see the
`SE_*` columns). The LR test addresses a different question: whether any
DIF is present, not how large it is.

**Real and artificial DIF.** Real DIF in one item induces artificial DIF
in the other items, in the opposite direction (Hagquist & Andrich,
2017). Here the mechanism is direct: each group's locations are centred
on that group's own mean threshold, so a real shift in one item moves
the other items' locations the other way. An item can therefore be
flagged only because another item has DIF, and one analysis cannot tell
real from artificial DIF. Hagquist and Andrich (2017) recommend
resolving DIF items one at a time, starting with the item with the
largest DIF (for instance by splitting it into group-specific items),
and repeating the analysis before deciding on the next. Do not
overinterpret the flags from a single run.

## References

Andersen, E. B. (1973). A goodness of fit test for the Rasch model.
*Psychometrika, 38*(1), 123-140.
[doi:10.1007/BF02291180](https://doi.org/10.1007/BF02291180)

Hagquist, C., & Andrich, D. (2017). Recent advances in analysis of
differential item functioning in health research using the Rasch model.
*Health and Quality of Life Outcomes, 15*, 181.
[doi:10.1186/s12955-017-0755-0](https://doi.org/10.1186/s12955-017-0755-0)

Henninger, M., Radek, J., Sengewald, M.-A., & Strobl, C. (2025). Partial
credit trees meet the partial gamma coefficient for quantifying DIF and
DSF in polytomous items. *Behaviormetrika, 52*(2), 221-257.
[doi:10.1007/s41237-024-00252-3](https://doi.org/10.1007/s41237-024-00252-3)

Jin, Y., Myers, N. D., Ahn, S., & Penfield, R. D. (2013). A comparison
of uniform DIF effect size estimators under the MIMIC and Rasch models.
*Educational and Psychological Measurement, 73*(2), 339-358.
[doi:10.1177/0013164412462705](https://doi.org/10.1177/0013164412462705)

Johansson, M., Preuter, M., & Melin, J. (2026). Valid and reliable?
Basic criteria and modern methods for assessing psychometric quality.
OSF Preprints.
[doi:10.31219/osf.io/3htzc_v2](https://doi.org/10.31219/osf.io/3htzc_v2)

Kim, S., Ouellet, L. J., Mazza, J., & Spaulding, A. C. (2017). Rasch
analysis and differential item functioning of a social support measure
in jail inmates with HIV infection. *Evaluation & the Health
Professions, 40*(1), 33-60.
[doi:10.1177/0163278716644954](https://doi.org/10.1177/0163278716644954)

Pickard, A. S., Dalal, M. R., & Bushnell, D. M. (2006). A comparison of
depressive symptoms in stroke and primary care: Applying Rasch models to
evaluate the Center for Epidemiologic Studies-Depression Scale. *Value
in Health, 9*(1), 59-64.
[doi:10.1111/j.1524-4733.2006.00082.x](https://doi.org/10.1111/j.1524-4733.2006.00082.x)

Strobl, C., Kopf, J., & Zeileis, A. (2015). Rasch trees: A new method
for detecting differential item functioning in the Rasch model.
*Psychometrika, 80*(2), 289-316.
[doi:10.1007/s11336-013-9388-3](https://doi.org/10.1007/s11336-013-9388-3)

Van de Winckel, A., Kozlowski, A. J., Johnston, M. V., Weaver, J.,
Grampurohit, N., Terhorst, L., Juengst, S., Ehrlich-Jones, L.,
Heinemann, A. W., Melvin, J., Sood, P., & Mallinson, T. (2022).
Reporting guideline for RULER: Rasch reporting guideline for
rehabilitation research: Explanation and elaboration. *Archives of
Physical Medicine and Rehabilitation, 103*(7), 1487-1498.
[doi:10.1016/j.apmr.2022.03.019](https://doi.org/10.1016/j.apmr.2022.03.019)

## See also

[`RMdifGamma`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md),
with simulation-based cutoffs from
[`RMdifGammaCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md),
and
[`RMdifTree`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifTree.md),
for continuous DIF variables and interactions between them (Strobl et
al., 2015; Henninger et al., 2025).

## Examples

``` r
# \donttest{
if (requireNamespace("eRm", quietly = TRUE)) { # for the example data
  set.seed(1)
  data("pcmdat2", package = "eRm")
  grp <- factor(sample(c("A", "B"), nrow(pcmdat2), replace = TRUE))

  # Default: kable of per-group item locations
  RMdifLR(pcmdat2, dif_var = grp)

  # ggplot panel of item locations with 95% CIs
  RMdifLR(pcmdat2, dif_var = grp, output = "ggplot")

  # Threshold-level kable, sorted by MaxDiff
  RMdifLR(pcmdat2, dif_var = grp, level = "threshold", sort = TRUE)

  # Tidy data.frame for downstream use
  df <- RMdifLR(pcmdat2, dif_var = grp, output = "dataframe")
  attr(df, "lr_test")
  df[df$Flagged, ]
}
#> [1] Item    A       B       All     MaxDiff Flagged SE_A    SE_B    SE_All 
#> <0 rows> (or 0-length row.names)
# }
```
