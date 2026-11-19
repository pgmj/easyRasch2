# DIF analysis via Andersen's likelihood-ratio test

Tests whether the item parameters differ between the groups of an
external variable, with Andersen's likelihood-ratio (LR) test. The model
is fitted by conditional maximum likelihood (CML) in the full sample and
in each group. The test answers whether there is any DIF, not which
items or thresholds have it. For that, use
[`RMdifAnchor`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifAnchor.md).

## Usage

``` r
RMdifLR(
  data,
  dif_var,
  model = c("auto", "PCM", "RM"),
  output = c("kable", "dataframe"),
  ...
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

  One of `"auto"` (default), `"PCM"`, or `"RM"`. `"auto"` reports the
  Rasch model when the data are dichotomous (max response = 1) and the
  partial credit model otherwise (the two coincide for dichotomous
  data). `"RM"` errors on polytomous data.

- output:

  `"kable"` (default) or `"dataframe"`.

- ...:

  Not used. Arguments that earlier versions accepted for the per-group
  location table and plots (`level`, `cutoff`, `conf`, `sort`,
  `p_value`, `correction`, `alpha`, `plot_type`, `reference`) give an
  error that points to
  [`RMdifAnchor`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifAnchor.md).

## Value

A one-row data.frame, or a `knitr_kable` of it, with `Chisq` (the LR
statistic), `df`, `p`, `Model` (`"PCM"` or `"RM"`), `Items`, `Groups`,
`N` (respondents used), `N_excluded` (rows dropped for a missing
`dif_var` or no responses), `N_incomplete` (respondents used who skipped
some items) and one `n_<group>` column per group.

## Details

The LR statistic is twice the difference between the sum of the group
conditional log-likelihoods and the full-sample one, with df equal to
the extra free parameters of the group fits (Andersen, 1973). It equals
[`eRm::LRtest()`](https://rdrr.io/pkg/eRm/man/LRtest.html) on the same
data. Each fit is a CML partial credit model from
[`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)
(a dichotomous item is a 2-category PCM). Respondents with missing
responses are kept (CML uses the items they answered).

The test does not depend on how each group's scale is centred, but any
comparison of single items or thresholds between groups does. Centring
each group on its own mean makes real DIF in some items appear as
artificial DIF in the others (Andrich & Hagquist, 2012), so this
function no longer reports per-group locations.
[`RMdifAnchor`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifAnchor.md)
links the groups through anchor items and tests each item, threshold or
threshold spacing, with a correction for multiple testing. The
simulation study behind this change is reported at
<https://pgmj.github.io/rasch_book/chapters/dif-anchor-tests.html>.

Like other global tests, the LR test becomes significant for small
amounts of DIF in large samples, and has little power in small ones.

**Response categories.** Every category an item uses must be used in
every group, or the group fits have different thresholds and the test is
not defined. The same holds when a category is used in a group only by
respondents with a zero or perfect score there, since CML leaves them
out. The function then stops and names the items, categories and groups
([`eRm::LRtest()`](https://rdrr.io/pkg/eRm/man/LRtest.html) instead left
items with an unused category out of the test without saying so). With
small groups this is common for sparse categories: merge them with an
adjacent one, combine groups, or leave the item out.

## References

Andersen, E. B. (1973). A goodness of fit test for the Rasch model.
*Psychometrika, 38*(1), 123-140.
[doi:10.1007/BF02291180](https://doi.org/10.1007/BF02291180)

Andrich, D., & Hagquist, C. (2012). Real and artificial differential
item functioning. *Journal of Educational and Behavioral Statistics,
37*(3), 387-416.
[doi:10.3102/1076998611411913](https://doi.org/10.3102/1076998611411913)

## See also

[`RMdifAnchor`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifAnchor.md)
for anchored tests of single items, thresholds and threshold spacing,
[`RMdifGamma`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md)
and
[`RMdifTree`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifTree.md).

## Examples

``` r
# \donttest{
if (requireNamespace("eRm", quietly = TRUE)) { # for the example data
  set.seed(1)
  data("pcmdat2", package = "eRm")
  grp <- factor(sample(c("A", "B"), nrow(pcmdat2), replace = TRUE))

  RMdifLR(pcmdat2, dif_var = grp)
  RMdifLR(pcmdat2, dif_var = grp, output = "dataframe")
}
#> RMdifLR() now reports Andersen's LR test only, without the per-group location table and plots. For DIF in single items, thresholds and threshold spacing, use RMdifAnchor().
#> This message is displayed once per session.
#>      Chisq df         p Model Items Groups   N N_excluded N_incomplete n_A n_B
#> 1 4.063775  7 0.7724031   PCM     4      2 300          0            0 146 154
# }
```
