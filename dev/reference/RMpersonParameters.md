# Person locations for a Rasch or partial credit model

Estimates a person location (theta) and its standard error of
measurement (SEM) for every respondent. Estimation is performed directly
on each person's observed response pattern, so partial missingness is
handled correctly: two respondents with the same sum score on different
subsets of items receive different estimates. This avoids the sum-score
lookup used by
[`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md),
which assumes complete data.

## Usage

``` r
RMpersonParameters(
  data,
  method = c("WLE", "EAP"),
  item_params = NULL,
  estimator = c("CML", "MML"),
  theta_range = c(-10, 10),
  prior_mean = 0,
  prior_sd = NULL,
  output = c("kable", "dataframe", "ggplot", "file"),
  filename = NULL
)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Missing values (`NA`) are
  allowed.

- method:

  Character. `"WLE"` (default) for Warm's weighted likelihood estimate,
  which has less bias than the MLE (Warm, 1989), or `"EAP"` for the
  expected a posteriori estimate under a normal prior (Bock & Mislevy,
  1982).

- item_params:

  Optional pre-specified item parameters, e.g. for anchoring or
  equating. Either a named list of Andrich-threshold vectors (one per
  item) or the long-format data.frame returned by
  [`RMitemParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemParameters.md).
  Items are matched to the columns of `data` by name, and `data` may
  hold extra columns. A name that is not a column of `data` stops with
  an error. When `NULL` (default) item parameters are estimated from
  `data`.

- estimator:

  Character. How item parameters are estimated when `item_params` is
  `NULL`: `"CML"` (default, via psychotools) or `"MML"` (via mirt).
  Ignored when `item_params` is supplied.

- theta_range:

  Numeric length 2. Search range for the WLE root and bounds for the EAP
  quadrature grid. Default `c(-10, 10)`. Only if the WLE root falls
  outside this range is the boundary returned.

- prior_mean:

  Numeric. Mean of the normal prior used by `method = "EAP"`. Default
  `0`, the natural choice when item parameters are centred at mean
  difficulty zero. Ignored for WLE.

- prior_sd:

  Numeric or `NULL`. Standard deviation of the normal prior used by
  `method = "EAP"`. When `NULL` (default) it is estimated from the data
  by marginal maximum likelihood, holding the item parameters fixed.
  Supply a number to fix it, for example `1` for an `N(0, 1)` prior, or
  a value carried over from a reference sample for equating. Ignored for
  WLE.

- output:

  Character. `"kable"` (default) for a formatted
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table,
  `"dataframe"` for the underlying data.frame, `"ggplot"` for a
  histogram of the estimated person locations, or `"file"` to write the
  data.frame to a CSV at `filename` (the data.frame is also returned
  invisibly).

- filename:

  Character. Path to the CSV file to write when `output = "file"`.
  Required in that case, ignored otherwise.

## Value

For `output = "dataframe"`, a data.frame with one row per respondent, in
input order (respondents with no responses are dropped with a message),
and the columns `theta` (the person location), `sem` (standard error of
measurement), `sum_score`, `n_answered` (number of responses) and
`extreme` (logical, the lowest or highest possible score on the items
answered). For `output = "kable"`, the same content as a `knitr_kable`
object. For `output = "ggplot"`, a histogram of `theta`.

With `method = "EAP"`, the prior used is attached to the data.frame as
`attr(result, "prior")`, a named numeric vector with `mean`, `sd` and
`estimated` (1 when `prior_sd` was estimated from the data, 0 when it
was supplied). The kable caption also states it. The CSV written by
`output = "file"` and the histogram do not include it, so keep the
data.frame, or the kable, if you need to report the prior or reuse an
estimated SD as `prior_sd`.

## Details

Item parameters are obtained once (by CML or MML, or taken from
`item_params`) and treated as fixed. Each person's theta is then
estimated from the items they answered.

**WLE** solves Warm's weighted-likelihood score equation
`l'(theta) + J(theta) / (2 I(theta)) = 0` by bracketed root finding,
where the bias-correction term `J / (2 I)` keeps the equation solvable
at the boundaries. The SEM is `1 / sqrt(I(theta))`. Unlike the MLE,
Warm's estimator therefore gives finite locations for the lowest and
highest possible scores, a step beyond the next score rather than +/-Inf
(Warm, 1989, Magis, 2015 for polytomous items), as in catR and TAM. Such
scores are still marked in the `extreme` column because they are
extrapolated and carry large standard errors. Only if the root lies
outside `theta_range` is the boundary returned with `NA` SEM.

**EAP** integrates the pattern likelihood against a normal prior over a
quadrature grid. The estimate is the posterior mean and the SEM the
posterior standard deviation. EAP is finite at extreme scores because
the prior shrinks them inward, at the cost of depending on the prior.
The prior used, including the estimate of `prior_sd` when it is left
`NULL`, is reported as described under Value.

## References

Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability
in a microcomputer environment. *Applied Psychological Measurement,
6*(4), 431-444.
[doi:10.1177/014662168200600405](https://doi.org/10.1177/014662168200600405)

Magis, D. (2015). A note on weighted likelihood and Jeffreys modal
estimation of proficiency levels in polytomous item response models.
*Psychometrika, 80*(1), 200-204.
[doi:10.1007/s11336-013-9378-5](https://doi.org/10.1007/s11336-013-9378-5)

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika, 54*(3), 427-450.
[doi:10.1007/BF02294627](https://doi.org/10.1007/BF02294627)

## See also

[`RMitemParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemParameters.md),
[`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md),
[`RMreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliability.md)

## Examples

``` r
# \donttest{
set.seed(1)
dat <- as.data.frame(
  matrix(sample(0:2, 200 * 6, replace = TRUE), nrow = 200, ncol = 6)
)
colnames(dat) <- paste0("Item", 1:6)
# Introduce some missingness
dat[cbind(sample(200, 30), sample(6, 30, replace = TRUE))] <- NA

# Default: WLE person locations
RMpersonParameters(dat, output = "dataframe") |> head()
#>        theta       sem sum_score n_answered extreme
#> 1 -0.2224510 0.5054892         5          6   FALSE
#> 2 -0.2224510 0.5054892         5          6   FALSE
#> 3 -1.0471828 0.6395816         2          6   FALSE
#> 4  0.2239135 0.5053530         7          6   FALSE
#> 5 -2.4725321 1.3253669         0          6    TRUE
#> 6 -1.0471828 0.6395816         2          6   FALSE

# EAP with a data-estimated prior SD
eap <- RMpersonParameters(dat, method = "EAP", output = "dataframe")
attr(eap, "prior") # mean, sd and whether the sd was estimated
#>      mean        sd estimated 
#> 0.0000000 0.1600579 1.0000000 

# EAP with a fixed N(0, 1) prior
RMpersonParameters(dat, method = "EAP", prior_sd = 1, output = "dataframe") |>
  head()
#>        theta       sem sum_score n_answered extreme
#> 1 -0.2162073 0.4680027         5          6   FALSE
#> 2 -0.2162073 0.4680027         5          6   FALSE
#> 3 -0.9457010 0.5321904         2          6   FALSE
#> 4  0.2169764 0.4679247         7          6   FALSE
#> 5 -1.6175238 0.6352806         0          6    TRUE
#> 6 -0.9457010 0.5321904         2          6   FALSE

# Write the person-location table to a CSV (also returned invisibly)
RMpersonParameters(dat, output = "file",
                   filename = tempfile(fileext = ".csv"))
#> Wrote 200 row(s) to '/tmp/RtmpSXG5Mp/file1b5135c11403.csv'.
# }
```
