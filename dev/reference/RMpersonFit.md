# Person fit for a Rasch or partial credit model

Computes person-fit statistics for each respondent, with p-values from
simulations under the fitted model. Three statistics are reported: the
conditional infit and outfit mean-squares (MSQ), the size measures
familiar from Rasch analysis, and the standardized log-likelihood
statistic lz. The conditional MSQ statistics use response probabilities
given the total score, which need no person estimate and are therefore
unbiased (Kreiner & Christensen, 2011). They are computed on each
person's observed responses, so missing responses are handled directly.

## Usage

``` r
RMpersonFit(
  data,
  statistics = c("infit", "outfit", "lz"),
  estimator = c("CML", "MML"),
  theta_method = c("WLE", "EAP"),
  iterations = 500L,
  flag_alpha = 0.05,
  flag = c("both", "underfit"),
  zstd = FALSE,
  parallel = FALSE,
  n_cores = NULL,
  seed = NULL,
  output = c("kable", "dataframe", "ggplot", "list")
)
```

## Arguments

- data:

  A data.frame or matrix of item responses, items scored from 0. Missing
  values (`NA`) are allowed.

- statistics:

  Character vector. Which statistics to report, any of `"infit"`,
  `"outfit"` and `"lz"`. Default all three.

- estimator:

  Character. How the item parameters are estimated: `"CML"` (default) or
  `"MML"` (mirt).

- theta_method:

  Character. Person-location estimator for lz: `"WLE"` (default) or
  `"EAP"`. Ignored if `"lz"` is not requested.

- iterations:

  Integer. Number of simulated response patterns per person. `0` skips
  the simulations and returns the statistics only. Default `500`.

- flag_alpha:

  Numeric in (0, 1). Significance level for the `flagged` column and the
  plots. Default `0.05`.

- flag:

  Character. Which direction of misfit is flagged for the MSQ
  statistics. `"both"` (default) flags either direction with a two-sided
  p-value. `"underfit"` flags only underfit (MSQ above 1, noisy
  responding) with a one-sided p-value, which has more power for it and
  ignores overfit, usually harmless. `p_infit` and `p_outfit` follow the
  choice. lz is always one-sided. `"both"` is the default, since overfit
  can occasionally signal fabricated or copied response patterns.

- zstd:

  Logical. Also report the Wilson-Hilferty standardized (ZSTD) MSQ
  statistics, as a familiar metric for comparing persons. They are not
  used for inference (Müller, 2020). Default `FALSE`.

- parallel, n_cores:

  Logical / integer. Run the simulations in parallel over persons via
  mirai when available. Default sequential. When `n_cores` is `NULL`,
  `getOption("mc.cores")` is used, and two workers when that is unset
  too.

- seed:

  Optional integer for reproducible simulations. See
  [easyRasch2-reproducibility](https://pgmj.github.io/easyRasch2/dev/reference/easyRasch2-reproducibility.md)
  for what this guarantees and how it interacts with `parallel`.

- output:

  Character. `"kable"` (default), `"dataframe"`, `"ggplot"` or `"list"`.
  `"ggplot"` returns a named list of person-fit maps, one per requested
  statistic (`infit`, `outfit`, `lz`), each plotting the statistic
  against the person location. Respondents flagged by that statistic are
  highlighted, and the caption gives the sample size and the share
  flagged. For the MSQ maps the flagged count is split into underfit
  (MSQ above 1, noisy responding) and overfit (MSQ below 1, overly
  predictable responding).

## Value

For `output = "dataframe"`, a data.frame with one row per respondent, in
input order: `id`, `n_answered`, `sum_score`, the requested statistics
(`infit_msq`, `outfit_msq`, `lz`), their p-values (`p_infit`,
`p_outfit`, `p_lz`) when `iterations > 0`, `flagged`, and with
`zstd = TRUE` also `infit_zstd` and `outfit_zstd`. `flagged` is `TRUE`
when the uncorrected p-value of any requested statistic is below
`flag_alpha`. It is a screening flag per person, not corrected for the
number of respondents, so when the model holds about `flag_alpha` of the
respondents are flagged by chance. Each map from `output = "ggplot"`
instead highlights and counts the respondents flagged by its own
statistic. Respondents with the lowest or highest possible score on the
items they answered get `NA` statistics and are not assessed.

For `output = "kable"`, the same content as a `knitr_kable`. For
`output = "ggplot"`, the named list of maps described under `output`.
For `output = "list"`, all from one computation: `fit` (the data.frame),
`plots` (the maps) and `kable` (the table from `output = "kable"`).

## Details

The p-values come from simulations under the model rather than from an
assumed asymptotic distribution. The asymptotic distributions of MSQ and
lz are unreliable, and the Wilson-Hilferty (ZSTD) transformation of MSQ
adds nothing once conditional estimation is used (Müller, 2020), so
simulation gives the valid reference (Sinharay, 2016a).

**Conditional MSQ.** For person \\v\\ with answered items \\A_v\\ and
total score \\r_v\\, the conditional residual is \\z\_{vi} = (x\_{vi} -
E(X\_{vi} \mid R_v = r_v)) / \sqrt{\mathrm{Var}(X\_{vi} \mid R_v =
r_v)}\\, with the conditional moments from elementary symmetric
functions of the answered items' parameters. Outfit is the unweighted
mean of \\z\_{vi}^2\\, and infit the information-weighted mean. The
conditional moments use only item parameters, so no biased person
estimate enters the residual.

**Direction.** Values near 1 indicate fit. MSQ above 1 (underfit) means
the responses are noisier than the model expects, the direction that
matters for validity, as in careless or aberrant responding. MSQ below 1
(overfit) means the responses are overly predictable (Guttman-like),
usually harmless for the scores, though a suspiciously perfect pattern
can occasionally signal copied or fabricated data. Use
`flag = "underfit"` to flag only the former. lz is one-sided: low
(negative) values flag the aberrant direction.

**lz** is the standardized log-likelihood of the response pattern at the
estimated person location (Drasgow, Levine & Williams, 1985). Small
(negative) values indicate misfit. Since the location is estimated, the
variance of lz is below 1 and its asymptotic standard normal
distribution does not hold, which makes a naive test conservative
(Snijders, 2001, Sinharay, 2016a). The analytic lz\* standardization
corrects for this, for dichotomous items (Snijders, 2001) and for
polytomous and mixed-format items (Sinharay, 2016b, 2026). easyRasch2
instead simulates the reference distribution with the person location
re-estimated for every simulated pattern (see Simulations), the
simulation counterpart of lz\*, so the p-value is calibrated. de la
Torre and Deng (2008) also replaced the asymptotic distribution of lz
with a simulated reference distribution.

**Simulations.** Each statistic is compared with the simulations that
match it. The conditional MSQ statistics use patterns simulated given
the person's total score, the exact null under the Rasch model, which
needs no person estimate. lz, which is defined at the estimated
location, uses patterns simulated at that location and then scored with
the location re-estimated from each simulated pattern, so the
simulations include the same effect of estimating the location as the
observed lz (in a Rasch or partial credit model the location is a
function of the total score, so this is a cheap lookup). Both are
computed per person over the items answered. The p-value is
`(1 + count) / (B + 1)`, with count the number of the B simulated values
at least as extreme as the observed one (two-sided: twice the smaller
tail). This follows the simulation-based approach to person fit of
Sinharay (2016a) and the recommendation of Müller (2020).

## References

de la Torre, J., & Deng, W. (2008). Improving person-fit assessment by
correcting the ability estimate and its reference distribution. *Journal
of Educational Measurement, 45*(2), 159-177.
[doi:10.1111/j.1745-3984.2008.00058.x](https://doi.org/10.1111/j.1745-3984.2008.00058.x)

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

Sinharay, S. (2016a). Assessment of person fit using resampling-based
approaches. *Journal of Educational Measurement, 53*(1), 63-85.
[doi:10.1111/jedm.12101](https://doi.org/10.1111/jedm.12101)

Sinharay, S. (2016b). Asymptotically correct standardization of
person-fit statistics beyond dichotomous items. *Psychometrika, 81*(4),
992-1013.
[doi:10.1007/s11336-015-9465-x](https://doi.org/10.1007/s11336-015-9465-x)

Sinharay, S. (2026). Refining the asymptotically correct standardization
of person-fit statistics for mixed-format tests. *British Journal of
Mathematical and Statistical Psychology*.
[doi:10.1111/bmsp.70049](https://doi.org/10.1111/bmsp.70049)

Snijders, T. A. B. (2001). Asymptotic null distribution of person fit
statistics with estimated person parameter. *Psychometrika, 66*(3),
331-342. [doi:10.1007/BF02294440](https://doi.org/10.1007/BF02294440)

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

# Conditional infit and outfit MSQ and lz, with simulated p-values
RMpersonFit(dat, iterations = 200, output = "dataframe") |> head()
#>   id n_answered sum_score infit_msq outfit_msq          lz   p_infit  p_outfit
#> 1  1          8         5 1.3900179  1.3988496  0.06189744 0.2786070 0.2885572
#> 2  2          8         6 0.7271090  0.7334992 -0.06110551 0.4079602 0.4477612
#> 3  3          8         2 0.8184723  0.8471554  0.27470409 0.7960199 0.7960199
#> 4  4          8        10 0.7119447  0.6871710  0.42107989 0.2587065 0.1791045
#> 5  5          8         4 1.5137117  1.4541483  0.33920680 0.2288557 0.2288557
#> 6  6          8         4 1.5321777  1.4836487  0.29890695 0.2189055 0.2189055
#>        p_lz flagged
#> 1 0.3930348   FALSE
#> 2 0.2835821   FALSE
#> 3 0.4228856   FALSE
#> 4 0.7860697   FALSE
#> 5 0.7064677   FALSE
#> 6 0.6119403   FALSE

# Person-fit maps: a named list with one plot per statistic
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plots <- RMpersonFit(dat, iterations = 200, output = "ggplot")
  plots$infit    # infit map (each plot's caption reports the % flagged)
}


# Flag only underfit (noisy responding), the validity-relevant direction
RMpersonFit(dat, iterations = 200, flag = "underfit",
            output = "dataframe") |> head()
#>   id n_answered sum_score infit_msq outfit_msq          lz    p_infit
#> 1  1          8         5 1.3900179  1.3988496  0.06189744 0.14925373
#> 2  2          8         6 0.7271090  0.7334992 -0.06110551 0.74129353
#> 3  3          8         2 0.8184723  0.8471554  0.27470409 0.49253731
#> 4  4          8        10 0.7119447  0.6871710  0.42107989 0.86069652
#> 5  5          8         4 1.5137117  1.4541483  0.33920680 0.10945274
#> 6  6          8         4 1.5321777  1.4836487  0.29890695 0.08457711
#>     p_outfit      p_lz flagged
#> 1 0.13432836 0.3482587   FALSE
#> 2 0.73631841 0.2537313   FALSE
#> 3 0.49253731 0.4626866   FALSE
#> 4 0.91044776 0.8059701   FALSE
#> 5 0.10945274 0.7213930   FALSE
#> 6 0.08457711 0.6119403   FALSE
# }
```
