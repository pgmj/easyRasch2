# Raw-Score to Logit Score Transformation Table

For a given set of items, returns the score-to-theta lookup that maps
each possible raw sum score to a person-location estimate (in logits)
and its standard error. Useful when reporting a scale's measurement
properties or converting raw totals to interval-scaled scores for
downstream analysis.

## Usage

``` r
RMscoreSE(
  data,
  method = "WLE",
  output = "kable",
  ci_multiplier = 1.96,
  point_size = 3,
  error_width = 0.5,
  theta_range = c(-10, 10)
)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Missing values (`NA`) are
  allowed; the underlying model fit handles them.

- method:

  Character string. Either `"WLE"` (default) for Warm's Weighted
  Likelihood Estimator computed from a CML-fitted Rasch / Partial Credit
  Model via `eRm`, or `"EAP"` for Expected A Posteriori sum-score
  estimates from an MML-fitted model via `mirt`.

- output:

  Character string controlling the return value: `"kable"` (default) for
  a formatted
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table,
  `"dataframe"` for the underlying data.frame, or `"ggplot"` for a
  `ggplot2` figure showing each raw score's logit estimate with
  `ci_multiplier`-scaled error bars.

- ci_multiplier:

  Numeric. Multiplier applied to the standard error to draw error bars
  on the figure. Default `1.96` (≈95% CI under a Gaussian
  approximation). Ignored when `output != "ggplot"`.

- point_size:

  Numeric. Point size for the figure. Default `3`.

- error_width:

  Numeric. Cap width for error bars on the figure. Default `0.5`.

- theta_range:

  Numeric length 2. Theta search range used for boundary raw scores
  under WLE estimation. Default `c(-10, 10)`. Ignored when
  `method = "EAP"`.

## Value

- If `output = "kable"`: a `knitr_kable` object with columns "Ordinal
  sum score", "Logit score", and "Logit std.error", and a caption noting
  the estimation method.

- If `output = "dataframe"`: a data.frame with columns `raw_score`,
  `logit_score`, and `logit_se` (one row per possible raw sum score from
  0 to the theoretical maximum).

- If `output = "ggplot"`: a `ggplot` object — points at each
  (`logit_score`, `raw_score`) with horizontal error bars at ±
  `ci_multiplier × logit_se`.

## Details

The function automatically detects whether the data is dichotomous (max
score 1) or polytomous (max score \> 1) and selects the appropriate
Rasch / Partial Credit model.

**`method = "WLE"`** fits the model with
[`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html) or
[`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html) (CML) and applies
Warm's Weighted Likelihood correction for the bias inherent in MLE-based
person estimates at the score boundaries. Score-boundary estimates are
obtained via the `theta_range` search; if a boundary cannot be solved
within that range it is returned as `Inf` / `-Inf` with `NA` SE. This
path uses lightly patched copies of the iarm 0.4.x person-estimation
code so that boundary scores converge cleanly across the full theta
range.

**`method = "EAP"`** fits the model with
`mirt::mirt(..., itemtype = "Rasch")` (MML) and obtains sum-score-based
EAP estimates and posterior SDs via
`mirt::fscores(method = "EAPsum", full.scores = FALSE, full.scores.SE = TRUE)`.
EAP estimates are finite at all score boundaries (the prior shrinks them
inward), but they depend on the assumed normal prior on theta. Item
parameters from MML differ slightly from the CML values used by the WLE
path; for well-behaved data the difference is small.

## References

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika, 54*(3), 427-450.
[doi:10.1007/BF02294627](https://doi.org/10.1007/BF02294627)

Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability
in a microcomputer environment. *Applied Psychological Measurement,
6*(4), 431-444.
[doi:10.1177/014662168200600405](https://doi.org/10.1177/014662168200600405)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:3, 200 * 6, replace = TRUE), nrow = 200, ncol = 6)
)
colnames(sim_data) <- paste0("Item", 1:6)

# Default kable output, WLE
RMscoreSE(sim_data)

# Underlying data.frame
RMscoreSE(sim_data, output = "dataframe")

# ggplot figure
RMscoreSE(sim_data, output = "ggplot")

# EAP via mirt
RMscoreSE(sim_data, method = "EAP")
} # }
```
