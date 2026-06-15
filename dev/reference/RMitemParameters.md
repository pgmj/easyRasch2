# Item Parameters for a Rasch / Partial Credit Model

Estimates item difficulty (dichotomous) or item-category threshold
(polytomous) parameters and returns them in long or wide format, with
optional standard errors and Wald confidence intervals. Item parameters
are estimated by conditional maximum likelihood (CML, via eRm) by
default, with marginal maximum likelihood (MML, via mirt) available for
sparse data where CML can be unstable.

## Usage

``` r
RMitemParameters(
  data,
  estimator = c("CML", "MML"),
  format = c("long", "wide"),
  se = TRUE,
  ci_level = 0.95,
  center = TRUE,
  output = c("kable", "dataframe")
)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Missing values (`NA`) are
  allowed; both estimators handle them.

- estimator:

  Character. `"CML"` (default) estimates item parameters with
  [`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html) /
  [`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html); `"MML"` uses
  [`mirt::mirt()`](https://philchalmers.github.io/mirt/reference/mirt.html)
  with `itemtype = "Rasch"`. CML is preferred for Rasch measurement; MML
  can be more robust when data are sparse. When `estimator = "CML"` and
  sparse response categories are detected, a warning suggests switching
  to MML.

- format:

  Character. `"long"` (default) returns one row per item (dichotomous)
  or per item-threshold (polytomous); `"wide"` returns one row per item
  with threshold columns `t1`, `t2`, ... plus a mean `location` column.

- se:

  Logical. If `TRUE` (default), standard-error and confidence-interval
  columns are added.

- ci_level:

  Numeric in (0, 1). Confidence level for the Wald interval
  (`estimate +/- z * SE`). Default `0.95`.

- center:

  Logical. If `TRUE` (default), all thresholds are shifted so their
  grand mean is zero, the usual Rasch identification. CML estimates are
  already centred; the shift mainly affects the MML path, keeping the
  two estimators on a common scale.

- output:

  Character. `"kable"` (default) for a formatted
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table, or
  `"dataframe"` for the underlying data.frame.

## Value

For `output = "dataframe"`, a data.frame. In **long** format the columns
are `item`, `threshold` (integer; `1` for dichotomous items),
`location`, and – when `se = TRUE` – `se`, `ci_lower`, `ci_upper`. In
**wide** format the columns are `item`, the threshold locations (`t1`,
`t2`, ... or `location` for dichotomous items), and a mean `location`;
when `se = TRUE`, matching `se_t1`, `se_t2`, ... columns are appended.
For `output = "kable"`, the same content as a `knitr_kable` object.

## Details

Items are detected as dichotomous (maximum score 1) or polytomous
(maximum score \> 1), and the Rasch or Partial Credit model is chosen
accordingly. Thresholds are reported as Andrich thresholds (the person
locations at which adjacent response categories are equally probable) on
the logit difficulty scale, matching
[`RMtargeting()`](https://pgmj.github.io/easyRasch2/dev/reference/RMtargeting.md).

**Standard errors.** For the CML path, dichotomous difficulty SEs come
directly from eRm. Polytomous threshold SEs are obtained by the delta
method: the full item-category parameter covariance is reconstructed as
`W V W'` (with `W` the eRm design matrix and `V` the basic-parameter
covariance) and propagated through the linear threshold map. For the MML
path, SEs come from the mirt parameter covariance, propagated through
the same map. Confidence intervals are Wald intervals and are symmetric
on the logit scale.

## References

Andrich, D. (1978). A rating formulation for ordered response
categories. *Psychometrika, 43*(4), 561-573.
[doi:10.1007/BF02293814](https://doi.org/10.1007/BF02293814)

Mair, P., & Hatzinger, R. (2007). Extended Rasch modeling: The eRm
package for the application of IRT models in R. *Journal of Statistical
Software, 20*(9), 1-20.
[doi:10.18637/jss.v020.i09](https://doi.org/10.18637/jss.v020.i09)

## See also

[`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md),
[`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md),
[`RMtargeting()`](https://pgmj.github.io/easyRasch2/dev/reference/RMtargeting.md)

## Examples

``` r
# \donttest{
set.seed(1)
poly <- as.data.frame(
  matrix(sample(0:2, 250 * 5, replace = TRUE), nrow = 250, ncol = 5)
)
colnames(poly) <- paste0("Item", 1:5)

# Default: long-format kable with SE and 95% CI
RMitemParameters(poly)
#> 
#> 
#> Table: Item thresholds via CML (Andrich thresholds, logit scale). SE and 95% Wald CI shown.
#> 
#> |item  | threshold| location|    se| ci_lower| ci_upper|
#> |:-----|---------:|--------:|-----:|--------:|--------:|
#> |Item1 |         1|   -0.014| 0.142|   -0.292|    0.264|
#> |Item1 |         2|    0.086| 0.163|   -0.233|    0.406|
#> |Item2 |         1|    0.165| 0.145|   -0.119|    0.449|
#> |Item2 |         2|   -0.101| 0.166|   -0.426|    0.225|
#> |Item3 |         1|    0.009| 0.141|   -0.268|    0.285|
#> |Item3 |         2|    0.113| 0.164|   -0.209|    0.435|
#> |Item4 |         1|   -0.217| 0.144|   -0.500|    0.066|
#> |Item4 |         2|    0.054| 0.158|   -0.255|    0.363|
#> |Item5 |         1|    0.203| 0.149|   -0.090|    0.496|
#> |Item5 |         2|   -0.299| 0.166|   -0.625|    0.027|

# Wide format, point estimates only
RMitemParameters(poly, format = "wide", se = FALSE, output = "dataframe")
#>    item      t1      t2 location
#> 4 Item4 -0.2167  0.0543  -0.0812
#> 5 Item5  0.2028 -0.2989  -0.0481
#> 2 Item2  0.1650 -0.1007   0.0321
#> 1 Item1 -0.0138  0.0864   0.0363
#> 3 Item3  0.0086  0.1131   0.0609

# Dichotomous data
dich <- as.data.frame(
  matrix(sample(0:1, 250 * 6, replace = TRUE), nrow = 250, ncol = 6)
)
colnames(dich) <- paste0("Item", 1:6)
RMitemParameters(dich, output = "dataframe")
#>    item threshold location     se ci_lower ci_upper
#> 1 Item1         1  -0.0539 0.1162  -0.2817   0.1738
#> 2 Item2         1  -0.0215 0.1161  -0.2491   0.2060
#> 3 Item3         1  -0.0377 0.1162  -0.2654   0.1899
#> 4 Item4         1   0.0754 0.1160  -0.1520   0.3028
#> 5 Item5         1   0.0593 0.1160  -0.1681   0.2867
#> 6 Item6         1  -0.0215 0.1161  -0.2491   0.2060
# }
```
