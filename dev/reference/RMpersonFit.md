# Person-Fit Statistics for a Rasch / Partial Credit Model

Computes per-respondent person-fit statistics and resampling-based
p-values. Three statistics are reported: conditional **infit** and
**outfit** mean-squares (MSQ) – the magnitude (effect-size) measures
familiar from Rasch analysis – and the standardized log-likelihood
statistic **lz**. The conditional MSQ statistics use response
probabilities conditional on the total score, which require no person
estimate and are therefore unbiased (Kreiner & Christensen); they are
computed on each person's observed response pattern, so partial
missingness is handled directly.

## Usage

``` r
RMpersonFit(
  data,
  statistics = c("infit", "outfit", "lz"),
  estimator = c("CML", "MML"),
  theta_method = c("WLE", "EAP"),
  resample = c("theta", "conditional"),
  iterations = 500L,
  flag_alpha = 0.05,
  zstd = FALSE,
  parallel = FALSE,
  n_cores = NULL,
  seed = NULL,
  output = c("kable", "dataframe", "ggplot")
)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Missing values (`NA`) are
  allowed.

- statistics:

  Character vector. Which statistics to report; any of `"infit"`,
  `"outfit"`, `"lz"`. Defaults to all three.

- estimator:

  Character. How item parameters are estimated: `"CML"` (default, via
  eRm) or `"MML"` (via mirt).

- theta_method:

  Character. Person-location estimator used for the lz statistic:
  `"WLE"` (default) or `"EAP"`. Ignored if `"lz"` is not requested.

- resample:

  Character. Resampling scheme for the null distribution: `"theta"`
  (default) simulates each person's responses at their estimated
  location; `"conditional"` samples response patterns conditional on the
  person's total score (no person estimate required, fully consistent
  with the conditional MSQ statistics).

- iterations:

  Integer. Number of Monte-Carlo replications per person. `0` skips
  resampling and returns statistics only. Default `500`.

- flag_alpha:

  Numeric in (0, 1). Two-sided significance level for the `flagged`
  column. Default `0.05`.

- zstd:

  Logical. If `TRUE`, also report the Wilson-Hilferty standardized
  (ZSTD) versions of the MSQ statistics. These are provided only as a
  familiar cross-person comparability metric and are **not** used for
  inference (see Müller, 2020). Default `FALSE`.

- parallel, n_cores:

  Logical / integer. Parallelise the resampling across persons via mirai
  when available. Default sequential.

- seed:

  Optional integer for reproducible resampling.

- output:

  Character. `"kable"` (default), `"dataframe"`, or `"ggplot"` (a
  person-fit map of outfit MSQ against person location, flagged
  respondents highlighted).

## Value

For `output = "dataframe"`, a data.frame with one row per respondent
(input order): `id`, `n_answered`, `sum_score`, the requested statistics
(`infit_msq`, `outfit_msq`, `lz`), their resampled p-values (`p_infit`,
`p_outfit`, `p_lz`) when `iterations > 0`, `flagged`, and – if
`zstd = TRUE` – `infit_zstd`, `outfit_zstd`. Extreme scorers (minimum or
maximum possible given the items they answered) receive `NA` statistics.
For `output = "kable"` the same content as a `knitr_kable`; for
`output = "ggplot"` a person-fit map.

## Details

Statistical significance is assessed by Monte-Carlo resampling under the
fitted model rather than by an assumed asymptotic distribution. The
asymptotic null distributions of MSQ and lz are known to be unreliable –
the Wilson-Hilferty (ZSTD) transformation of MSQ in particular adds
nothing once conditional estimation is used (Müller, 2020) – so
resampling provides the valid reference (Sinharay, 2016).

**Conditional MSQ.** For person \\v\\ with answered items \\A_v\\ and
total score \\r_v\\, the conditional residual is \\z\_{vi} = (x\_{vi} -
E(X\_{vi} \mid R_v = r_v)) / \sqrt{\mathrm{Var}(X\_{vi} \mid R_v =
r_v)}\\, with conditional moments obtained from elementary symmetric
functions of the answered items' parameters. Outfit is the unweighted
mean of \\z\_{vi}^2\\; infit is the information-weighted mean. Because
the conditional moments use only item parameters, no (biased) person
estimate enters the residual.

**lz.** The standardized log-likelihood of the response pattern
evaluated at the estimated person location (Drasgow, Levine & Williams,
1985); small (negative) values indicate misfit. Its asymptotic
standard-normal null is unreliable when the person location is
estimated; the resampled p-value is therefore the basis for inference,
which makes the analytic lz\\ standardization (Snijders, 2001)
unnecessary here.

**Resampling.** For each person, `iterations` response patterns are
simulated under the fitted model – either at the estimated location
(`resample = "theta"`) or conditional on the observed total score
(`resample = "conditional"`) – and each statistic is recomputed. The
two-sided p-value is the proportion of replicates at least as extreme as
the observed value. This follows the resampling-based person-fit
approach (Sinharay, 2016) and the bootstrap recommendation of Müller
(2020).

## References

Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness
measurement with polychotomous item response models and standardized
indices. *British Journal of Mathematical and Statistical Psychology,
38*(1), 67-86.
[doi:10.1111/j.2044-8317.1985.tb00817.x](https://doi.org/10.1111/j.2044-8317.1985.tb00817.x)

Kreiner, S., & Christensen, K. B. (2011). Exact evaluation of bias in
Rasch model residuals. *Advances in Mathematics Research, 12*, 19-40.

Müller, M. (2020). Item fit statistics for Rasch analysis: can we trust
them? *Journal of Statistical Distributions and Applications, 7*(5).
[doi:10.1186/s40488-020-00108-7](https://doi.org/10.1186/s40488-020-00108-7)

Sinharay, S. (2016). Assessment of person fit using resampling-based
approaches. *Journal of Educational Measurement, 53*(1), 63-85.
[doi:10.1111/jedm.12101](https://doi.org/10.1111/jedm.12101)

## See also

[`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md),
[`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md),
[`RMitemParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemParameters.md)

## Examples

``` r
# \donttest{
set.seed(1)
dat <- as.data.frame(
  matrix(sample(0:2, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
)
colnames(dat) <- paste0("Item", 1:8)

# Conditional infit/outfit MSQ + lz with resampled p-values
RMpersonFit(dat, iterations = 200, output = "dataframe") |> head()
#>   id n_answered sum_score infit_msq outfit_msq      lz p_infit p_outfit  p_lz
#> 1  1          8         5    1.3900     1.3989  0.0619  0.2424   0.2727 0.525
#> 2  2          8         6    0.7271     0.7335 -0.0611  0.3300   0.3500 0.470
#> 3  3          8         2    0.8185     0.8472  0.2747  0.8108   0.8865 0.570
#> 4  4          8        10    0.7120     0.6872  0.4211  0.3500   0.3200 0.580
#> 5  5          8         4    1.5137     1.4542  0.3392  0.1443   0.1753 0.620
#> 6  6          8         4    1.5322     1.4837  0.2989  0.1218   0.1320 0.595
#>   flagged
#> 1   FALSE
#> 2   FALSE
#> 3   FALSE
#> 4   FALSE
#> 5   FALSE
#> 6   FALSE

# Compare the two resampling schemes
a <- RMpersonFit(dat, resample = "theta",       iterations = 200,
                 output = "dataframe")
b <- RMpersonFit(dat, resample = "conditional", iterations = 200,
                 output = "dataframe")
# }
```
