# Reliability metrics for a Rasch model

Computes three reliability indices for a Rasch / partial credit model:
the Person Separation Index (PSI) – WLE-based separation reliability –
the marginal reliability (native, CML test information integrated over
the estimated normal latent density), and Relative Measurement
Uncertainty (RMU) via
[`RMUreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMUreliability.md)
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
  [`RMUreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMUreliability.md)
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
  [easyRasch2-reproducibility](https://pgmj.github.io/easyRasch2/dev/reference/easyRasch2-reproducibility.md)
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
estimated normal latent density \\N(0, \sigma^2)\\ (\\\sigma\\ from
marginal ML). Integrating over the estimated latent variance, rather
than the \\N(0,1)\\ assumed by
[`mirt::marginal_rxx()`](https://philchalmers.github.io/mirt/reference/marginal_rxx.html),
keeps it correct on the Rasch logit scale, where \\\sigma\\ is typically
well above 1 and the \\N(0,1)\\ assumption underestimates reliability.

It is the model-based complement to the sample-based PSI, and the two
are now the same coefficient by two routes: PSI divides by the observed
spread of the WLE estimates, marginal reliability by the fitted latent
density. A large gap between them therefore does flag an off-target or
non-normal sample. Through version 1.2.0 this row used Green's
subtractive \\1 - \overline{1/I(\theta)}/\sigma^2\\ instead, under which
much of the PSI-to-marginal gap was an artefact of the differing
formulas rather than a property of the sample. See
`dev/TODO-reliability-form.md`.

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
[`RMreliabilityCurve()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliabilityCurve.md)
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

[`RMUreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMUreliability.md),
[`RMreliabilityCurve()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliabilityCurve.md)

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
#> Table: Reliability for 20 items, n = 100 respondents. PSI is the WLE-based separation reliability and excludes min/max scoring respondents.
#> 
#> |Metric                | Estimate| Lower (95% HDCI)| Upper (95% HDCI)|Notes                       |
#> |:---------------------|--------:|----------------:|----------------:|:---------------------------|
#> |Cronbach's alpha      |    0.754|            0.701|            0.813|25 bootstrap resamples      |
#> |PSI                   |    0.725|            0.673|            0.768|25 bootstrap resamples      |
#> |Marginal (ratio form) |    0.769|            0.717|            0.817|25 bootstrap resamples      |
#> |RMU (WLE)             |    0.753|            0.684|            0.817|1000 PVs, 50 RMU iterations |
# }
```
