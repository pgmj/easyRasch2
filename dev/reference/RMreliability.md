# Reliability coefficients for a Rasch model

Computes four reliability coefficients for a Rasch or partial credit
model: Cronbach's alpha, the person separation index (PSI, the WLE-based
separation reliability), the marginal reliability (the mean of the
conditional reliability curve over the estimated latent distribution),
and the relative measurement uncertainty (RMU,
[`RMUreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMUreliability.md))
from plausible values drawn with
[`mirt::fscores()`](https://philchalmers.github.io/mirt/reference/fscores.html).

## Usage

``` r
RMreliability(
  data,
  conf_int = 0.95,
  draws = 1000,
  rmu_iter = 50,
  estim = "WLE",
  boot = FALSE,
  boot_iter = 200,
  parallel = TRUE,
  n_cores = NULL,
  seed = NULL,
  verbose = FALSE,
  theta_range = c(-10, 10),
  output = "kable"
)
```

## Arguments

- data:

  A data.frame or matrix of item responses, items scored from 0. See
  Details for how missing values are handled.

- conf_int:

  Numeric in (0, 1). HDCI width of the bootstrap intervals and of RMU.
  Default `0.95`.

- draws:

  Integer. Number of plausible values drawn from the `mirt` model for
  RMU. Default `1000`. More draws give a more stable RMU, at a cost that
  grows about linearly.

- rmu_iter:

  Integer. How many times
  [`RMUreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMUreliability.md)
  is repeated on the same draws, each with a new random split of the
  columns. The estimates are averaged over the repetitions. Default
  `50`.

- estim:

  Character. Estimator used by
  [`mirt::fscores()`](https://philchalmers.github.io/mirt/reference/fscores.html)
  as the starting point for the plausible values: `"WLE"` (default),
  `"EAP"`, `"MAP"` or `"ML"`. The plausible values themselves are drawn
  by Metropolis-Hastings. PSI and the marginal reliability do not use
  it.

- boot:

  Logical. Bootstrap confidence intervals for alpha, PSI and the
  marginal reliability. Default `FALSE`.

- boot_iter:

  Integer. Number of bootstrap samples with `boot = TRUE`. Default
  `200`.

- parallel:

  Logical. Run the bootstrap in parallel via mirai when it is installed.
  Default `TRUE`.

- n_cores:

  Integer or `NULL`. Number of parallel workers. When `NULL`,
  `getOption("mc.cores")` is used. If that is unset too, the bootstrap
  runs sequentially.

- seed:

  Integer or `NULL`. Random seed for reproducibility. See
  [easyRasch2-reproducibility](https://pgmj.github.io/easyRasch2/dev/reference/easyRasch2-reproducibility.md)
  for what this guarantees and how it interacts with `parallel`.

- verbose:

  Logical. Print progress messages and a progress bar for the bootstrap.
  Default `FALSE`.

- theta_range:

  Numeric vector of length 2. Theta limits passed to
  [`mirt::fscores()`](https://philchalmers.github.io/mirt/reference/fscores.html).
  Default `c(-10, 10)`.

- output:

  Character. `"kable"` (default) for a
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table, or
  `"dataframe"`.

## Value

- `output = "kable"`: a `knitr_kable` table with one row per
  coefficient.

- `output = "dataframe"`: a data.frame with the columns `metric`
  (`"Cronbach's alpha"`, `"PSI"`, `"Marginal (curve mean)"` and
  `"RMU (<estim>)"`), `estimate`, `lower`, `upper` and `notes`.

## Details

With `boot = TRUE`, alpha, PSI and the marginal reliability get
confidence intervals from a nonparametric bootstrap that resamples
respondents and recomputes the coefficients without refitting a `mirt`
model. The RMU interval is the HDCI of the correlations across
plausible-value draws, averaged over `rmu_iter` random splits of the
draws.

**Data.** Respondents with no responses are dropped. Cronbach's alpha
uses complete cases. PSI and the marginal reliability use every
respondent's answered items, and PSI leaves out respondents with the
lowest or highest possible score (see below). At least one complete case
is needed.

**Marginal reliability** is the mean of the conditional reliability
curve over the latent density, \\\int \sigma^2/(\sigma^2 +
1/I(\theta))\\ g(\theta)\\d\theta\\. The test information \\I(\theta)\\
is summed from the CML item parameters, and \\g\\ is the normal latent
density \\N(\mu, \sigma^2)\\ with both parameters estimated by marginal
maximum likelihood with the items fixed. Integrating over the estimated
distribution, rather than the \\N(0,1)\\ that
[`mirt::marginal_rxx()`](https://philchalmers.github.io/mirt/reference/marginal_rxx.html)
assumes, keeps it correct on the Rasch logit scale, where \\\sigma\\ is
often well above 1 and an \\N(0,1)\\ assumption underestimates
reliability.

Up to version 1.3.0 \\\mu\\ was held at 0. Centring the item thresholds
fixes the item mean, not where the respondents are, so a normal density
pinned to 0 could reach an off-target sample only by widening.
\\\sigma\\ then absorbed the mistargeting, and marginal reliability rose
as targeting got worse. See the note under Value in
[`RMreliabilityCurve()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliabilityCurve.md)
for the size of the effect.

The marginal reliability is the model-based counterpart of the
sample-based PSI. PSI divides by the observed spread of the WLE
estimates, the marginal reliability uses the fitted latent density. A
large gap between them flags an off-target or non-normal sample, and
since version 1.3.1 both fall when targeting gets worse.

**How it relates to published coefficients.** The reported value is none
of the coefficients named in the literature.

*It is the mean of a curve, not a ratio of averages.* The classical
latent-scale reliability is one variance ratio, \\\sigma^2/(\sigma^2 +
\overline{SEM^2})\\ (Milanzi et al., 2015, section 3.2, with a constant
error variance). Here the ratio is averaged over the latent density
instead. Since \\\sigma^2/(\sigma^2 + x)\\ is convex in \\x\\, the mean
of the curve is the larger of the two (Jensen's inequality). The choice
was made on accuracy: averaging the curve tracked the exact
expected-sum-score reliability more closely than the ratio of averages
in an internal replication (mean absolute error 0.015 against 0.018 for
binary data, 0.010 against 0.015 for polytomous).

*It is a latent-scale quantity checked against a manifest criterion.*
Milanzi et al. (2015) argue that interest usually lies in the
reliability of observed scores, and that latent coefficients run
consistently higher. This one was therefore validated against the exact
expected-sum-score reliability, \\Var(\mu)/(Var(\mu) + Var(\epsilon))\\,
and tracked it to within about 0.02 in the conditions tested. It remains
a latent-scale coefficient that approximates a manifest one.

*Milanzi et al. recommend something this package does not implement.*
They favour Taylor-series manifest reliability measures, which
approximate \\Var(\mu)\\ and \\Var(\epsilon)\\ on the observed-score
scale directly. In the same internal replication those were the most
accurate (0.012). They are cited here for their criticism of the
subtractive coefficient, which this package acted on, not for their
remedy.

Up to version 1.2.0 this row used the subtractive coefficient \\1 -
\overline{1/I(\theta)}/\sigma^2\\ (Green et al., 1984, Milanzi et al.'s
equation 12, attributed there to Lord, 1980). It can fall below zero
when the average error variance exceeds the trait variance, and much of
the gap between PSI and that coefficient came from the different
formulas rather than the sample. It is still returned as the
`marginal_green` attribute of
[`RMreliabilityCurve()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliabilityCurve.md)
for comparison.

**PSI** is the WLE-based separation reliability, \\1 - \overline{SEM^2}
/ \mathrm{Var}(\hat\theta)\\, from the CML item thresholds
(`psychotools`) and the WLE person locations and their SEMs. Respondents
with the lowest or highest possible raw score are left out, since their
boundary estimates would inflate the person variance and overstate
reliability. Earlier versions used
[`eRm::SepRel()`](https://rdrr.io/pkg/eRm/man/SepRel.html) with MLE, and
the values can differ, most for scales with many extreme scores, such as
dichotomous items.

**RMU** follows Bignardi, Kievit and Bürkner (2026), here with plausible
values from `mirt` rather than fully Bayesian posterior draws (see
Mislevy, 1991, for plausible values).

Bootstrap samples that cannot be fitted are dropped, with a message
giving how many and why. The usual reason is a sparse lowest or middle
response category that is missing from the sample.

## References

Bignardi, G., Kievit, R., & Bürkner, P. C. (2026). A general method for
estimating reliability using Bayesian Measurement Uncertainty.
*PsyArXiv*. <https://osf.io/preprints/psyarxiv/h54k8>

Green, B. F., Bock, R. D., Humphreys, L. G., Linn, R. L., & Reckase, M.
D. (1984). Technical Guidelines for Assessing Computerized Adaptive
Tests. *Journal of Educational Measurement, 21*(4), 347–360.
[doi:10.1111/j.1745-3984.1984.tb01039.x](https://doi.org/10.1111/j.1745-3984.1984.tb01039.x)

Milanzi, E., Molenberghs, G., Alonso, A., Verbeke, G., & De Boeck, P.
(2015). Reliability measures in item response theory: Manifest versus
latent correlation functions. *British Journal of Mathematical and
Statistical Psychology, 68*(1), 43-64.
[doi:10.1111/bmsp.12033](https://doi.org/10.1111/bmsp.12033)

Mislevy, R. J. (1991). Randomization-Based Inference about Latent
Variables from Complex Samples. *Psychometrika, 56*(2), 177-196.
[doi:10.1007/BF02294457](https://doi.org/10.1007/BF02294457)

## See also

[`RMreliabilityCurve()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliabilityCurve.md)
for precision along the scale,
[`RMUreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMUreliability.md),
[`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md).

## Examples

``` r
# \donttest{
if (requireNamespace("ggdist", quietly = TRUE) &&
    requireNamespace("eRm", quietly = TRUE)) {
  set.seed(1)
  RMreliability(eRm::raschdat1[, 1:20], draws = 1000)

  # Bootstrap intervals for alpha, PSI and the marginal reliability
  # (use 200 or more bootstrap samples in real analyses)
  RMreliability(eRm::raschdat1[, 1:20], draws = 1000,
                boot = TRUE, boot_iter = 25, parallel = FALSE, seed = 42)
}
#> 
#> 
#> Table: Reliability for 20 items, n = 100 respondents. PSI is the WLE-based separation reliability and excludes min/max scoring respondents. Marginal is the latent-density-weighted mean of the conditional reliability curve, sigma^2 / (sigma^2 + SEM(theta)^2).
#> 
#> |Metric                | Estimate| Lower (95% HDCI)| Upper (95% HDCI)|Notes                       |
#> |:---------------------|--------:|----------------:|----------------:|:---------------------------|
#> |Cronbach's alpha      |    0.754|            0.701|            0.813|25 bootstrap resamples      |
#> |PSI                   |    0.725|            0.673|            0.768|25 bootstrap resamples      |
#> |Marginal (curve mean) |    0.748|            0.699|            0.798|25 bootstrap resamples      |
#> |RMU (WLE)             |    0.753|            0.684|            0.817|1000 PVs, 50 RMU iterations |
# }
```
