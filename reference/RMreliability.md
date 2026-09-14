# Reliability metrics for a Rasch model

Computes three reliability indices for a Rasch / partial credit model:
the Person Separation Index (PSI) – WLE-based separation reliability –
the marginal reliability (native, CML test information integrated over
the estimated normal latent density), and Relative Measurement
Uncertainty (RMU) via
[`RMUreliability()`](https://pgmj.github.io/easyRasch2/reference/RMUreliability.md)
applied to plausible values from
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

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers).

- conf_int:

  Numeric in (0, 1). HDCI width for both bootstrap CIs and RMU. Default
  `0.95`.

- draws:

  Integer. Number of plausible-value draws drawn from the mirt model for
  the RMU calculation. Default `1000`. More gives a more stable RMU;
  computational cost is mostly linear.

- rmu_iter:

  Integer. Number of times
  [`RMUreliability()`](https://pgmj.github.io/easyRasch2/reference/RMUreliability.md)
  is repeated on the same set of plausible-value draws (each repetition
  uses a fresh random column split). Estimates are averaged across
  repetitions to stabilise against split-induced variability. Default
  `50`.

- estim:

  Character. Theta estimator used by
  [`mirt::fscores()`](https://philchalmers.github.io/mirt/reference/fscores.html)
  for the RMU plausible-value seed. One of `"WLE"` (default), `"EAP"`,
  `"MAP"`, `"ML"`. Plausible draws themselves are produced by
  Metropolis-Hastings. (PSI and marginal reliability are computed
  natively and do not use this.)

- boot:

  Logical. If `TRUE`, run a non-parametric bootstrap to obtain CIs for
  PSI and Marginal reliability. Default `FALSE`.

- boot_iter:

  Integer. Number of bootstrap iterations when `boot = TRUE`. Default
  `200`.

- parallel:

  Logical. Use parallel processing via `mirai` for the bootstrap if
  available. Default `TRUE`.

- n_cores:

  Integer or `NULL`. Number of parallel workers. When `NULL`,
  `getOption("mc.cores")` is checked first; if neither is set,
  bootstrapping falls back to sequential.

- seed:

  Integer or `NULL`. Master random seed for reproducibility. See
  [easyRasch2-reproducibility](https://pgmj.github.io/easyRasch2/reference/easyRasch2-reproducibility.md)
  for what this guarantees and how it interacts with `parallel`.

- verbose:

  Logical. Print progress messages and a progress bar for the bootstrap.
  Default `FALSE`.

- theta_range:

  Numeric length-2 vector. Theta limits passed to
  [`mirt::fscores()`](https://philchalmers.github.io/mirt/reference/fscores.html).
  Default `c(-10, 10)`.

- output:

  Character. `"kable"` (default) for a formatted
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table, or
  `"dataframe"` for the underlying data.frame.

## Value

- If `output = "kable"`: a `knitr_kable` object with one row per metric.

- If `output = "dataframe"`: a data.frame with columns `metric`,
  `estimate`, `lower`, `upper`, `notes`.

## Details

Confidence intervals for **PSI** and **Marginal** reliability are
obtained by non-parametric bootstrap (resampling respondents; all three
indices are recomputed natively per resample, no model is refitted by
`mirt`). The RMU interval is the HDCI of correlations across
plausible-value draws, averaged over `rmu_iter` random splits of the
draws.

Marginal reliability is the latent-density-weighted mean of the
conditional reliability curve, \\\int \sigma^2/(\sigma^2 +
1/I(\theta))\\ g(\theta)\\d\theta\\, where the test information
\\I(\theta)\\ is summed from the CML item parameters and \\g\\ is the
estimated normal latent density \\N(\mu, \sigma^2)\\, **both**
parameters coming from marginal ML with the items held fixed.
Integrating over the estimated latent distribution, rather than the
\\N(0,1)\\ assumed by
[`mirt::marginal_rxx()`](https://philchalmers.github.io/mirt/reference/marginal_rxx.html),
keeps it correct on the Rasch logit scale, where \\\sigma\\ is typically
well above 1 and the \\N(0,1)\\ assumption underestimates reliability.

Through version 1.3.0 \\\mu\\ was held at 0 rather than estimated. That
is an assumption, not a consequence of centring the item thresholds,
which fixes the *item* mean and says nothing about where the respondents
sit. A normal density pinned to 0 can only reach an off-target sample by
widening, so \\\sigma\\ absorbed the mistargeting and marginal
reliability **rose** as targeting worsened. See the note under Value in
[`RMreliabilityCurve()`](https://pgmj.github.io/easyRasch2/reference/RMreliabilityCurve.md)
for the size of the effect.

It is the model-based complement to the sample-based PSI, and the two
are now the same coefficient by two routes: PSI divides by the observed
spread of the WLE estimates, marginal reliability by the fitted latent
density. A large gap between them therefore does flag an off-target or
non-normal sample, and since 1.3.1 both move in the same direction when
targeting worsens rather than apart.

**How this relates to the published coefficients.** Three things are
worth being explicit about, because the reported value is not any of the
coefficients the literature names.

*It is the mean of a curve, not a ratio of averages.* The classical
latent-scale reliability is a single variance ratio,
\\\sigma^2/(\sigma^2 + \overline{SEM^2})\\ (Milanzi et al., 2015,
section 3.2, where the error variance is a constant). What is reported
here averages the ratio over the latent density instead. Since
\\\sigma^2/(\sigma^2 + x)\\ is convex in \\x\\, the mean of the curve is
the **larger** of the two by Jensen's inequality. The choice was made on
accuracy, not convention: averaging the curve tracked the exact
expected-sum-score reliability more closely than the ratio of averages
in both arms of an internal replication (mean absolute error 0.015
against 0.018 for binary data, 0.010 against 0.015 for polytomous).

*It is a latent-scale quantity checked against a manifest criterion.*
Milanzi et al. (2015) argue that interest usually lies in the
reliability of observed scores rather than latent ones, and that latent
coefficients run consistently higher. This one was therefore validated
against the exact expected-sum-score reliability, \\Var(\mu)/(Var(\mu) +
Var(\epsilon))\\, rather than assumed to match it, and tracked it to
within about 0.02 across the conditions tested. It remains a
latent-scale coefficient that approximates a manifest one, not a
manifest coefficient.

*Milanzi et al. recommend something this package does not implement.*
Their conclusion favours Taylor-series manifest reliability measures,
which approximate \\Var(\mu)\\ and \\Var(\epsilon)\\ on the
observed-score scale directly rather than working from test information.
In the same internal replication those were the most accurate of the
estimators compared (0.012). They are cited here for their criticism of
the subtractive coefficient, which this package acted on, and not for
their remedy, which it has not adopted.

The superseded subtractive coefficient is Milanzi et al.'s equation
(12), attributed there to Lord (1980), and is still available as the
`marginal_green` attribute of
[`RMreliabilityCurve()`](https://pgmj.github.io/easyRasch2/reference/RMreliabilityCurve.md).
Through version 1.2.0 this row used Green's subtractive \\1 -
\overline{1/I(\theta)}/\sigma^2\\ instead, under which much of the
PSI-to-marginal gap was an artefact of the differing formulas rather
than a property of the sample. See `dev/TODO-reliability-form.md`.

PSI is the WLE-based separation reliability, \\1 - \overline{SEM^2} /
\mathrm{Var}(\hat\theta)\\, computed from CML item thresholds
(`psychotools`) and Warm's WLE person locations / analytic SEMs.
Respondents with extreme (min/max) raw scores are excluded – their
boundary estimates would inflate the person variance and overstate
reliability. (Earlier versions used
[`eRm::SepRel()`](https://rdrr.io/pkg/eRm/man/SepRel.html) with MLE; the
values can differ, most noticeably for scales with many extreme scorers,
e.g. dichotomous items.)

RMU is from Bignardi, Kievit, & Bürkner (2025), modified here to use
mirt plausible values rather than fully Bayesian posterior draws (see
Mislevy, 1991, for the plausible-values framework).

Marginal reliability here is the subtractive Green/Lord coefficient.
Milanzi et al. (2015) show that this form can fall below zero when the
average error variance exceeds the trait variance, which happens with
few items or a narrow sample, and it is floored at 0 above.
[`RMreliabilityCurve()`](https://pgmj.github.io/easyRasch2/reference/RMreliabilityCurve.md)
reports the same quantity in the bounded ratio form alongside this one,
so the two can be compared directly. See `dev/TODO-reliability-form.md`
for the open question of which form this row should use.

Bootstrap iterations that fail to converge are silently dropped.

## References

Bignardi, G., Kievit, R., & Bürkner, P. C. (2025). A general method for
estimating reliability using Bayesian Measurement Uncertainty.
*PsyArXiv*.
[doi:10.31234/osf.io/h54k8_v1](https://doi.org/10.31234/osf.io/h54k8_v1)

Green, B. F., Bock, R. D., Humphreys, L. G., Linn, R. L., & Reckase, M.
D. (1984). Technical Guidelines for Assessing Computerized Adaptive
Tests. *Journal of Educational Measurement, 21*(4), 347–360.
[doi:10.1111/j.1745-3984.1984.tb01039.x](https://doi.org/10.1111/j.1745-3984.1984.tb01039.x)

Mislevy, R. J. (1991). Randomization-Based Inference about Latent
Variables from Complex Samples. *Psychometrika, 56*(2), 177-196.
[doi:10.1007/BF02294457](https://doi.org/10.1007/BF02294457)

Adams, R. J. (2005). Reliability as a measurement design effect.
*Studies in Educational Evaluation, 31*(2), 162-172.
[doi:10.1016/j.stueduc.2005.05.008](https://doi.org/10.1016/j.stueduc.2005.05.008)

Milanzi, E., Molenberghs, G., Alonso, A., Verbeke, G., & De Boeck, P.
(2015). Reliability measures in item response theory: Manifest versus
latent correlation functions. *British Journal of Mathematical and
Statistical Psychology, 68*(1), 43-64.
[doi:10.1111/bmsp.12033](https://doi.org/10.1111/bmsp.12033)

## See also

[`RMUreliability()`](https://pgmj.github.io/easyRasch2/reference/RMUreliability.md),
[`RMreliabilityCurve()`](https://pgmj.github.io/easyRasch2/reference/RMreliabilityCurve.md)

## Examples

``` r
# \donttest{
if (requireNamespace("ggdist", quietly = TRUE) &&
    requireNamespace("eRm", quietly = TRUE)) {
  set.seed(1)
  RMreliability(eRm::raschdat1[, 1:20], draws = 1000)

  # Bootstrap CI for PSI and Marginal
  # (use more bootstrap iterations, e.g. 200+, in real analyses)
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
