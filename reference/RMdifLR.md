# DIF analysis via Andersen's likelihood-ratio test

Splits a Rasch model by an external grouping variable using
[`eRm::LRtest()`](https://rdrr.io/pkg/eRm/man/LRtest.html) and reports
per-group item locations (or per-group threshold locations) together
with their standard errors. A single function replaces the four legacy
helpers (`RIdifTableLR`, `RIdifThreshTblLR`, `RIdifFigureLR`,
`RIdifThreshFigLR`) by exposing the two underlying axes – `level` (item
or threshold) and `output` (data.frame, kable, or ggplot) – as
arguments. The same data preparation pipeline feeds all six
combinations.

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

  One of `"auto"` (default), `"PCM"`, or `"RM"`. `"auto"` fits
  [`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html) when the data are
  dichotomous (max response = 1) and
  [`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html) otherwise. `"RM"`
  errors on polytomous data.

- level:

  One of `"item"` (default) or `"threshold"`. `"item"` reports each
  item's mean threshold location per group; `"threshold"` reports each
  individual threshold per group. For dichotomous (RM) data the two
  views are equivalent (one threshold per item).

- output:

  One of `"kable"` (default), `"dataframe"`, or `"ggplot"`. The
  data.frame view always carries a `Flagged` logical column; the kable
  view bolds flagged rows; the ggplot view shows confidence intervals.

- cutoff:

  Numeric or `NULL`. Threshold (in logits) for the `Flagged` column: a
  row is flagged when `MaxDiff > cutoff`, where `MaxDiff` is the
  difference between the largest and smallest per-group location for
  that item (or threshold). Set to `NULL` to suppress flagging. Default
  `0.5`.

- conf:

  Numeric in (0, 1). Confidence level used for the ggplot error bars.
  Default `0.95`.

- sort:

  Logical. `kable` output only: sort rows by `MaxDiff` (descending).
  Default `FALSE`.

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

For the data.frame and kable outputs, locations are reported on the
centred eRm parameterisation returned by
[`eRm::thresholds()`](https://rdrr.io/pkg/eRm/man/thresholds.html).
Per-group fits come from `eRm::LRtest(..., splitcr = dif_var)`; the
unsplit fit (`All` column) is the model fitted to the full dataset. The
Andersen LR statistic, df, and p-value reported as the `lr_test`
attribute / caption come directly from `LRtest()`'s return value.

`cell_spec()`-style HTML cell colouring used in the legacy easyRasch
package has been dropped in favour of a logical `Flagged` column (and
bold rendering in the kable output), so the kable renders correctly in
HTML, LaTeX, and pipe/markdown.

## Examples

``` r
# \donttest{
if (requireNamespace("eRm", quietly = TRUE)) {
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
