# Occasion-to-occasion SD from a test-retest study

Estimates the per-occasion SD, in logits, of a respondent's
occasion-specific deviation: the variation between two administrations
that is neither measurement error nor change in the construct. It is the
`retest_sd` argument of
[`RMpersonChange()`](https://pgmj.github.io/easyRasch2/reference/RMpersonChange.md)
when `null = "retest"`.

## Usage

``` r
RMretestSD(
  data_t1,
  data_t2,
  anchor = "stack",
  item_params = NULL,
  method = "WLE",
  estimator = "CML",
  sim_iter = 500,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  boot = TRUE,
  boot_iter = 500,
  conf_int = 0.95,
  seed = NULL,
  theta_range = c(-10, 10),
  output = "kable"
)
```

## Arguments

- data_t1, data_t2:

  Data.frames or matrices of item responses from a test-retest study,
  same items in the same order, row `i` being the same respondent at
  both administrations.

- anchor, item_params, method, estimator, theta_range:

  As in
  [`RMpersonChange()`](https://pgmj.github.io/easyRasch2/reference/RMpersonChange.md).
  Use the same settings there as here, since the estimate is defined
  relative to the measurement-error variance the model implies.

- sim_iter:

  Integer. Iterations used to simulate the measurement-error component.
  Default `500`.

- parallel:

  Logical. Use `mirai` for the simulation if available. Default `TRUE`.

- n_cores:

  Integer or `NULL`. Parallel workers. When `NULL`,
  `getOption("mc.cores")` is checked first.

- verbose:

  Logical. Progress bar for the simulation. Default `FALSE`.

- boot:

  Logical. Bootstrap a confidence interval by resampling respondents.
  Default `TRUE`.

- boot_iter:

  Integer. Bootstrap iterations. Default `500`.

- conf_int:

  Numeric in (0, 1). HDCI width. Default `0.95`.

- seed:

  Integer or `NULL`. Random seed.

- output:

  Character. `"kable"` (default) or `"dataframe"`.

## Value

A one-row data.frame (or its `knitr_kable`) with `variance`, `sd`,
`lower`, `upper`, `var_change` (the observed variance of the change) and
`var_error` (the simulated measurement-error variance of the change),
plus `n`.

## Details

Writing \\\hat\theta\_{it} = \theta_i + u\_{it} + e\_{it}\\ with
\\u\_{it} \sim N(0, \sigma^2\_{retest})\\ the occasion deviation and
\\e\_{it}\\ the measurement error, the observed change has variance
\\2\sigma^2\_{retest} + SE\_{i1}^2 + SE\_{i2}^2\\. So

\$\$\hat\sigma^2\_{retest} = \frac{\mathrm{Var}(d) -
\mathrm{Var}\_{err}(d)}{2}\$\$

**How the error term is obtained.** The obvious choice for
\\\mathrm{Var}\_{err}(d)\\ is \\\overline{SE_1^2 + SE_2^2}\\ from the
information-based standard errors, and that is what Caronni et al.
(2026) and the classical formulations use. It is badly biased on short
scales, because the information standard error is asymptotic and
overstates the real sampling variance of \\\hat\theta\\ when items are
few, so the subtraction removes too much. Simulation over a range of
scale lengths, at a true occasion variance of 0.09 and n = 2000:

|           |                |               |
|-----------|----------------|---------------|
| **items** | **asymptotic** | **simulated** |
| 6         | -0.023         | 0.072         |
| 12        | 0.053          | 0.093         |
| 20        | 0.086          | 0.099         |

Across eight replications at n = 1000 the simulated version is close to
unbiased on a six-item scale (mean 0.002 against a truth of 0, and 0.081
against a truth of 0.09), with a replication SD near 0.015.

This function therefore obtains \\\mathrm{Var}\_{err}(d)\\ by
simulation: respondents are held at their observed locations, two
response vectors are drawn with no occasion variance, and the variance
of the resulting change is averaged over `sim_iter` iterations. The
observed missingness pattern is carried into the simulated responses.
`var_error` in the output is that simulated quantity, not the mean of
the squared standard errors.

The retest interval has to be long enough for recall to fade and short
enough that no real change occurs. That is an assumption about the
design and the data cannot check it.

**The estimate can be negative.** It is a variance obtained by
subtraction, so when the observed change is tighter than the measurement
model says it should be, the difference goes below zero. A negative
value is reported as such rather than floored, because flooring hides
the signal. It points at standard errors overstated by misfit or local
dependence, at regression toward a common value over the interval, or at
sampling noise. `sd` is `NA` whenever `variance` is negative.

It is a variance estimated from `n` differences, so its own sampling
error is roughly \\\sigma^2\sqrt{2/(n-1)}\\. A few dozen respondents
will not pin it down, and the bootstrap interval is the honest way to
see that. The bootstrap resamples respondents and recomputes
\\\mathrm{Var}(d)\\ only, holding the simulated error term at its
full-sample value, since that term is model-implied and its Monte Carlo
error is small beside the sampling error of the change variance.

**Converting a published test-retest coefficient.** A published
test-retest ICC or Pearson correlation can be converted, but it is not
itself the quantity wanted here, and neither is a published test-retest
SEM (\\SD\sqrt{1-ICC}\\), which already contains measurement error.
Since \\\hat\theta = \theta + u + e\\ gives \\r = \mathrm{Var}(\theta) /
\mathrm{Var}(\hat\theta)\\,

\$\$\hat\sigma^2\_{retest} = (1 - r)\\\mathrm{Var}(\hat\theta) -
\mathrm{Var}\_{err}(d)/2\$\$

with \\\mathrm{Var}(\hat\theta)\\ the observed variance of the person
estimates and \\\mathrm{Var}\_{err}(d)\\ the simulated error term this
function reports as `var_error`. Checked against known truth at n =
3000, true occasion variance 0.09:

|  |  |  |  |  |
|----|----|----|----|----|
| **items** | **r (logit)** | **r (sum score)** | **from logit r** | **from sum-score r** |
| 6 | .808 | .829 | 0.094 | 0.046 |
| 12 | .873 | .890 | 0.084 | 0.048 |
| 20 | .903 | .913 | 0.088 | 0.067 |

Three conditions have to hold, and the first is usually violated:

- **The coefficient must be on the logit metric.** The score-to-logit
  map is nonlinear, so a coefficient computed on sum scores runs
  consistently higher than the same coefficient on logits. After
  subtraction the residual is small, so that modest gap does real
  damage: on a six-item scale, substituting the sum-score value roughly
  halves the recovered occasion variance. The error has a known
  direction, understating occasion noise and leaving
  [`RMpersonChange()`](https://pgmj.github.io/easyRasch2/reference/RMpersonChange.md)
  too permissive. Nearly all published coefficients are on sum scores.

- **The published sample must resemble yours.** A retest coefficient
  rises with the heterogeneity of the people measured, so it is a
  property of that sample rather than of the instrument. Combining
  someone else's `r` with your \\\mathrm{Var}(\hat\theta)\\ assumes the
  structure transfers.

- **Pearson and an agreement ICC are not interchangeable.** Pearson
  ignores a systematic mean shift between occasions and so reports
  stability that practice or drift has removed. ICC(A,1) penalises it
  and is the closer choice here.

The conversion also needs the simulated error term, which this function
produces only when given two occasions of raw data. Converting from a
published coefficient and a single administration is not implemented.
When that is the situation you are in, `retest_sd_tip` from
[`RMpersonChange()`](https://pgmj.github.io/easyRasch2/reference/RMpersonChange.md)
answers the practical question from the other end: it reports how much
occasion noise each flagged result would tolerate, which you can weigh
against a published coefficient without converting anything.

## See also

[`RMpersonChange()`](https://pgmj.github.io/easyRasch2/reference/RMpersonChange.md)

## Examples

``` r
# \donttest{
set.seed(1)
thr <- lapply(seq(-1.2, 1.2, length.out = 6), function(b) b + c(-0.8, 0, 0.8))
theta <- rnorm(200, 0, 1.4)
# a true occasion SD of 0.3 logits
r1 <- as.data.frame(easyRasch2:::sim_partial_score(thr, theta + rnorm(200, 0, 0.3)))
r2 <- as.data.frame(easyRasch2:::sim_partial_score(thr, theta + rnorm(200, 0, 0.3)))
colnames(r1) <- colnames(r2) <- paste0("I", 1:6)
RMretestSD(r1, r2, sim_iter = 100, parallel = FALSE, boot = FALSE)
#> 
#> 
#> Table: Per-occasion SD of occasion-specific deviation, in logits, for use as `retest_sd` in RMpersonChange(null = "retest"). Estimated as (Var(change) - mean error variance) / 2. n = 200 respondents measured twice.
#> 
#> |   n| Var(change)| Mean error var| Occasion var| Occasion SD| Lower (95% HDCI)| Upper (95% HDCI)|
#> |---:|-----------:|--------------:|------------:|-----------:|----------------:|----------------:|
#> | 200|      0.7728|         0.6471|       0.0629|       0.251|               NA|               NA|
# }
```
