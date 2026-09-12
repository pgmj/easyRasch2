# Single-subject change between two occasions

Tests, for each respondent, whether their person location moved between
two occasions by more than measurement error allows. The statistic is
the Rasch change index, \\RCI = (\hat\theta_2 -
\hat\theta_1)/SE\_{diff}\\, referred by default to a simulated rather
than a normal null.

## Usage

``` r
RMpersonChange(
  data_t1,
  data_t2,
  id = NULL,
  anchor = "stack",
  item_params = NULL,
  method = "WLE",
  estimator = "CML",
  null = "measurement",
  retest_sd = NULL,
  critical = "exact",
  alpha = 0.05,
  direction = "two.sided",
  conditional_crit = FALSE,
  sim_iter = 1000,
  parallel = TRUE,
  n_cores = NULL,
  seed = NULL,
  theta_range = c(-10, 10),
  verbose = FALSE,
  output = "dataframe"
)
```

## Arguments

- data_t1, data_t2:

  Data.frames or matrices of item responses at the two occasions. Same
  items, in the same order, one row per respondent, with row `i` of each
  being the same person. Items must be scored from 0.

- id:

  Optional vector of respondent identifiers, length `nrow(data_t1)`.
  Defaults to row numbers.

- anchor:

  Character. Where the item calibration comes from: `"stack"` (default,
  both occasions calibrated together), `"t1"`, or `"t2"`. Ignored when
  `item_params` is supplied.

- item_params:

  Optional pre-specified item parameters, either a named list of
  Andrich-threshold vectors or the long-format data.frame from
  [`RMitemParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemParameters.md),
  as in
  [`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md).
  **Required for a single respondent**, since item parameters cannot be
  estimated from one person. Supplied thresholds are used as given and
  are not re-centred.

- method:

  Character. Person-location estimator, `"WLE"` (default) or `"EAP"`.

- estimator:

  Character. How item parameters are estimated when `item_params` is
  `NULL`: `"CML"` (default) or `"MML"`.

- null:

  Character. Which null hypothesis is tested. `"measurement"` (default)
  treats the response process as the only source of variation.
  `"retest"` additionally treats occasion-to-occasion fluctuation as
  noise and requires `retest_sd`. See Details, and note that the two
  answer different questions.

- retest_sd:

  Numeric, or `NULL` (default). The **per-occasion** SD, in logits, of a
  respondent's occasion-specific deviation. Required when
  `null = "retest"` and not accepted otherwise. Estimate it with
  [`RMretestSD()`](https://pgmj.github.io/easyRasch2/dev/reference/RMretestSD.md).

- critical:

  How the critical values are obtained. `"exact"` (default) enumerates
  the null distribution, `"simulate"` estimates it by Monte Carlo, and a
  positive number is used as a symmetric cutoff on the RCI referred to a
  normal null (for example `1.96`).

- alpha:

  Numeric in (0, 1). Error rate for a single respondent. Default `0.05`.

- direction:

  Character. `"two.sided"` (default), or `"increase"` / `"decrease"` for
  a one-sided test. These name the direction of \\\theta\\; whether an
  increase is an improvement depends on how the scale is oriented.

- conditional_crit:

  Logical. If `TRUE`, critical values are simulated separately for each
  respondent at their own location rather than pooled across the sample.
  Default `FALSE`. Ignored when `critical` is numeric.

- sim_iter:

  Integer. Simulation iterations when `critical = "simulate"`. Default
  `1000`.

- parallel:

  Logical. Use `mirai` for the simulation if available. Default `TRUE`.

- n_cores:

  Integer or `NULL`. Parallel workers. When `NULL`,
  `getOption("mc.cores")` is checked first; if neither is set, the
  simulation runs sequentially.

- seed:

  Integer or `NULL`. Random seed. See
  [easyRasch2-reproducibility](https://pgmj.github.io/easyRasch2/dev/reference/easyRasch2-reproducibility.md).

- theta_range:

  Numeric length 2. Search range for the WLE root and bounds for the EAP
  grid. Default `c(-10, 10)`.

- verbose:

  Logical. Print a progress bar for the simulation. Default `FALSE`.

- output:

  Character. `"dataframe"` (default), `"kable"`, or `"ggplot"`.

## Value

- If `output = "dataframe"`: one row per respondent, with columns `id`,
  `sum_t1`, `sum_t2`, `theta_t1`, `se_t1`, `theta_t2`, `se_t2`,
  `extreme_t1`, `extreme_t2`, `change`, `se_diff`, `rci`, `p_value`,
  `crit_lower`, `crit_upper`, `change_class`, and `retest_sd_tip`.
  Attributes `null`, `retest_sd`, `anchor`, `alpha`, `direction`,
  `critical` and `sim_iter` record the analysis.

- If `output = "kable"`: the same content as a `knitr_kable`.

- If `output = "ggplot"`: occasion 1 against occasion 2, with the
  no-change band and points coloured by `change_class`.

## Details

**Which null is tested.** Following Zumbo (2026), the estimand is fixed
by what counts as systematic and what counts as residual, and that
choice has to be made before an estimator is picked. Two are available
here:

- `null = "measurement"`, the default. The respondent's \\\theta\\ is
  fixed and identical at both occasions and only the response process is
  random, so \\SE\_{diff}^2 = SE_1^2 + SE_2^2\\. This asks whether the
  two estimates differ by more than responding alone would produce.

- `null = "retest"`. Occasion-to-occasion fluctuation that is not change
  in the construct (state, mood, recall, practice) also counts as noise,
  so \\SE\_{diff}^2 = SE_1^2 + SE_2^2 + 2\sigma^2\_{retest}\\. The
  factor of 2 is there because `retest_sd` is a per-occasion SD and both
  occasions carry one.

The default is the narrower of the two, so it flags more change than a
retest-based null would. It is the only one computable from two response
vectors alone, which is why it is the default, but a significant result
under it is a necessary condition for change rather than evidence of it.
`retest_sd_tip` reports, for each flagged respondent, how large
\\\sigma\_{retest}\\ would have to be to overturn their result.

Caronni et al. (2026) compute \\1.96\sqrt{SE_A^2 + SE_B^2}\\ from the
Rasch score-to-measure map and call it a minimal detectable change. That
is the `null = "measurement"` quantity, whereas the classical MDC is
built from a test-retest ICC and is the `null = "retest"` quantity. The
two are not interchangeable, which is why the null in force is named in
every caption here.

Regression to the mean is present under both nulls and is addressed by
neither. A regression-based reliable change index (Maassen, 2004) is out
of scope for a single-subject function, since it needs population
parameters the individual case does not supply.

**Why 1.96 is the wrong cutoff.** The sum score is a sufficient
statistic, so \\\hat\theta\\ takes only as many values as there are
scores, and \\SE(\hat\theta)\\ is a deterministic function of the score
rather than a constant. The denominator of the RCI is therefore not
independent of its numerator: a respondent who moves to an extreme score
produces a large change and a large standard error together, and the
ratio is damped. Both tails of the null are pulled in, so the critical
value sits below 1.96, and the shortfall grows as the scale shortens.
Enumerated values for one six-item family: 1.62 at four items, 1.70 at
six, 1.79 at ten, 1.85 at twenty, 1.92 at forty. Referring the RCI to a
normal null is conservative, and materially so on short scales.

**How the null is obtained.** `critical = "exact"` enumerates it.
Because the score is sufficient, \\\hat\theta(r)\\ and \\SE(r)\\ are
deterministic lookups and \\P(r \mid \theta)\\ follows from the
Lord-Wingersky recursion, so the null is a discrete distribution over
score pairs with no Monte Carlo error. Respondents are grouped by their
pair of answered-item sets, so incomplete data is enumerated within each
pattern. Under `null = "retest"` the occasion deviations are integrated
out by quadrature. `critical = "simulate"` estimates the same
distribution by drawing response vectors instead, which is slower and
noisier but makes no use of sufficiency.

Quantiles are taken from the **signed** RCI, so `crit_lower` and
`crit_upper` are reported separately. Under the null the two occasions
are exchangeable, which makes the RCI symmetric about zero however
skewed \\\hat\theta\\ itself is, and the two bounds then agree in
magnitude. They part company when exchangeability fails, most commonly
when the occasions differ in which items were answered.

The pooled critical value is a property of the sample as well as the
items, because it averages over where the respondents sit. Targeting
matters more than the shape of the distribution: skewing the sample
barely moves it, while pushing the same sample off target until 40
percent of respondents score at an extreme moved it from 1.70 to 1.52 in
one check.

**What the RCI measures, and what it does not.** The RCI answers whether
a change is distinguishable from the null. It is not a measure of how
much someone changed, and the two orderings genuinely differ. Its value
mixes two ingredients that it cannot separate: how far the respondent
moved, and how precisely each of their two positions was pinned down. On
a six-item scale scored 0 to 3 the standard error runs from 0.49 in the
middle to 1.45 at the boundary, a factor of three within one instrument.
A respondent going from the minimum to the maximum score moves 7.24
logits and scores RCI 3.53, while one going from 3 to 15 moves 3.11
logits and scores 3.59. The first traversed the whole scale, the second
less than half of it, and the second has the larger statistic. Both
statements are correct about their own question.

In a group comparison every case shares one standard error, so ordering
by the test statistic matches ordering by effect. Here the scale factor
changes from respondent to respondent, and the output is read one person
at a time, which is where the misreading does damage. Read `change` for
magnitude, in logits, and `rci` only for the decision. Do not rank
respondents by `rci`, and do not compare RCIs across instruments or
across regions of one scale as though they were a common currency.

Simulating fixes the reference distribution, not the estimand. A
simulated critical value under `null = "measurement"` is a well
calibrated answer to the measurement-error question, not to the retest
question.

**Calibration and its precondition.** `anchor = "stack"` calibrates on
both occasions row-bound together. It uses all the data and puts both
occasions on one metric by construction, at the cost of assuming the
items behave the same way on both occasions. When they do not, the
compromise calibration pulls the occasions toward each other and shrinks
apparent change. Test that precondition first with the package's DIF
tools, using occasion as the grouping variable on the stacked data.
`anchor = "t1"` measures follow-up on the baseline metric and keeps
post-treatment data out of it, at the cost of halving the calibration
sample. `"t2"` is the mirror.

Note that stacking puts each respondent in the data twice, so the item
parameters' own standard errors are optimistic. That does not propagate
into the person standard errors, which treat the thresholds as fixed.

**`retest_sd_tip`.** For a respondent whose change is flagged, this is
the per-occasion retest SD that would bring their RCI back to the
critical value, \\\sigma^2\_{tip} = (\mathrm{change}^2/c^2 - SE_1^2 -
SE_2^2)/2\\. It is `NA` for respondents not flagged, for whom the
question does not arise. The critical value is held fixed at the one in
force, which is exact when `critical` is numeric and an approximation
when it is simulated, since a simulated critical value drifts toward
\\\pm 1.96\\ as the added normal component smooths the discreteness.

## References

Caronni, A., et al. (2026). Improving single-subject change assessment:
deriving the minimal detectable change of questionnaires' ordinal scores
from the Rasch analysis measures.

Jacobson, N. S., & Truax, P. (1991). Clinical significance: A
statistical approach to defining meaningful change in psychotherapy
research. *Journal of Consulting and Clinical Psychology, 59*(1), 12-19.
[doi:10.1037/0022-006X.59.1.12](https://doi.org/10.1037/0022-006X.59.1.12)

Maassen, G. H. (2004). The standard error in the Jacobson and Truax
Reliable Change Index. *Journal of Clinical and Experimental
Neuropsychology, 26*(5), 643-657.
[doi:10.1080/13803390409609791](https://doi.org/10.1080/13803390409609791)

Zumbo, B. D. (2026). Conditional standard error of measurement as an
estimand of individual score precision. *Psychometrika*.
[doi:10.1017/psy.2026.10141](https://doi.org/10.1017/psy.2026.10141)

## See also

[`RMretestSD()`](https://pgmj.github.io/easyRasch2/dev/reference/RMretestSD.md),
[`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md),
[`RMreliabilityCurve()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliabilityCurve.md),
[`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md)

## Examples

``` r
# \donttest{
set.seed(1)
thr <- lapply(seq(-1.2, 1.2, length.out = 6), function(b) b + c(-0.8, 0, 0.8))
theta <- rnorm(80, 0, 1.4)
t1 <- as.data.frame(easyRasch2:::sim_partial_score(thr, theta))
t2 <- as.data.frame(easyRasch2:::sim_partial_score(thr, theta + 0.8))
colnames(t1) <- colnames(t2) <- paste0("I", 1:6)

RMpersonChange(t1, t2)
#>    id sum_t1 sum_t2    theta_t1     se_t1    theta_t2     se_t2 extreme_t1
#> 1   1      5      8 -1.05265787 0.5342416 -0.30246523 0.5039273      FALSE
#> 2   2      9     12 -0.05692244 0.5055699  0.72136345 0.5361040      FALSE
#> 3   3      6     10 -0.79407316 0.5164849  0.19368859 0.5113013      FALSE
#> 4   4     12     17  0.72136345 0.5361040  2.76643368 0.9204042      FALSE
#> 5   5     10     15  0.19368859 0.5113013  1.68336335 0.6470375      FALSE
#> 6   6      5     10 -1.05265787 0.5342416  0.19368859 0.5113013      FALSE
#> 7   7     12     14  0.72136345 0.5361040  1.32134359 0.5929126      FALSE
#> 8   8     15     13  1.68336335 0.6470375  1.00745673 0.5585944      FALSE
#> 9   9     12     14  0.72136345 0.5361040  1.32134359 0.5929126      FALSE
#> 10 10      6     10 -0.79407316 0.5164849  0.19368859 0.5113013      FALSE
#> 11 11     14     16  1.32134359 0.5929126  2.13324117 0.7384993      FALSE
#> 12 12     10     15  0.19368859 0.5113013  1.68336335 0.6470375      FALSE
#> 13 13      7      7 -0.54647281 0.5071109 -0.54647281 0.5071109      FALSE
#> 14 14      0      2 -3.76951909 1.4711205 -2.05231802 0.6992678       TRUE
#> 15 15     17     16  2.76643368 0.9204042  2.13324117 0.7384993      FALSE
#> 16 16      7     13 -0.54647281 0.5071109  1.00745673 0.5585944      FALSE
#> 17 17      8     14 -0.30246523 0.5039273  1.32134359 0.5929126      FALSE
#> 18 18     16     16  2.13324117 0.7384993  2.13324117 0.7384993      FALSE
#> 19 19     12     17  0.72136345 0.5361040  2.76643368 0.9204042      FALSE
#> 20 20      7     15 -0.54647281 0.5071109  1.68336335 0.6470375      FALSE
#> 21 21     15     17  1.68336335 0.6470375  2.76643368 0.9204042      FALSE
#> 22 22     17     15  2.76643368 0.9204042  1.68336335 0.6470375      FALSE
#> 23 23      9     10 -0.05692244 0.5055699  0.19368859 0.5113013      FALSE
#> 24 24      1      1 -2.61813570 0.8735700 -2.61813570 0.8735700      FALSE
#> 25 25     11     17  0.45204624 0.5211334  2.76643368 0.9204042      FALSE
#> 26 26     12     14  0.72136345 0.5361040  1.32134359 0.5929126      FALSE
#> 27 27      6     13 -0.79407316 0.5164849  1.00745673 0.5585944      FALSE
#> 28 28      3      3 -1.65451471 0.6135234 -1.65451471 0.6135234      FALSE
#> 29 29      7     14 -0.54647281 0.5071109  1.32134359 0.5929126      FALSE
#> 30 30     13     15  1.00745673 0.5585944  1.68336335 0.6470375      FALSE
#> 31 31     16     18  2.13324117 0.7384993  4.00101566 1.5358175      FALSE
#> 32 32      9     10 -0.05692244 0.5055699  0.19368859 0.5113013      FALSE
#> 33 33     12     14  0.72136345 0.5361040  1.32134359 0.5929126      FALSE
#> 34 34      7     12 -0.54647281 0.5071109  0.72136345 0.5361040      FALSE
#> 35 35      3      2 -1.65451471 0.6135234 -2.05231802 0.6992678      FALSE
#> 36 36      8     14 -0.30246523 0.5039273  1.32134359 0.5929126      FALSE
#> 37 37      5     11 -1.05265787 0.5342416  0.45204624 0.5211334      FALSE
#> 38 38      8     10 -0.30246523 0.5039273  0.19368859 0.5113013      FALSE
#> 39 39     16     16  2.13324117 0.7384993  2.13324117 0.7384993      FALSE
#> 40 40     12     16  0.72136345 0.5361040  2.13324117 0.7384993      FALSE
#> 41 41     11     10  0.45204624 0.5211334  0.19368859 0.5113013      FALSE
#> 42 42      8      8 -0.30246523 0.5039273 -0.30246523 0.5039273      FALSE
#> 43 43     13     15  1.00745673 0.5585944  1.68336335 0.6470375      FALSE
#> 44 44     16     15  2.13324117 0.7384993  1.68336335 0.6470375      FALSE
#> 45 45      4      8 -1.33332113 0.5641476 -0.30246523 0.5039273      FALSE
#> 46 46      6      8 -0.79407316 0.5164849 -0.30246523 0.5039273      FALSE
#> 47 47     11     15  0.45204624 0.5211334  1.68336335 0.6470375      FALSE
#> 48 48     14     16  1.32134359 0.5929126  2.13324117 0.7384993      FALSE
#> 49 49     13     12  1.00745673 0.5585944  0.72136345 0.5361040      FALSE
#> 50 50     14     15  1.32134359 0.5929126  1.68336335 0.6470375      FALSE
#> 51 51     11     11  0.45204624 0.5211334  0.45204624 0.5211334      FALSE
#> 52 52      3     11 -1.65451471 0.6135234  0.45204624 0.5211334      FALSE
#> 53 53      7     16 -0.54647281 0.5071109  2.13324117 0.7384993      FALSE
#> 54 54      3     12 -1.65451471 0.6135234  0.72136345 0.5361040      FALSE
#> 55 55     17     17  2.76643368 0.9204042  2.76643368 0.9204042      FALSE
#> 56 56     18     17  4.00101566 1.5358175  2.76643368 0.9204042       TRUE
#> 57 57      6      9 -0.79407316 0.5164849 -0.05692244 0.5055699      FALSE
#> 58 58      7      8 -0.54647281 0.5071109 -0.30246523 0.5039273      FALSE
#> 59 59     13     14  1.00745673 0.5585944  1.32134359 0.5929126      FALSE
#> 60 60      7     11 -0.54647281 0.5071109  0.45204624 0.5211334      FALSE
#> 61 61     17     18  2.76643368 0.9204042  4.00101566 1.5358175      FALSE
#> 62 62     10     12  0.19368859 0.5113013  0.72136345 0.5361040      FALSE
#> 63 63     13     18  1.00745673 0.5585944  4.00101566 1.5358175      FALSE
#> 64 64      9     16 -0.05692244 0.5055699  2.13324117 0.7384993      FALSE
#> 65 65      1      4 -2.61813570 0.8735700 -1.33332113 0.5641476      FALSE
#> 66 66      7     16 -0.54647281 0.5071109  2.13324117 0.7384993      FALSE
#> 67 67      1      1 -2.61813570 0.8735700 -2.61813570 0.8735700      FALSE
#> 68 68     16     18  2.13324117 0.7384993  4.00101566 1.5358175      FALSE
#> 69 69     11     13  0.45204624 0.5211334  1.00745673 0.5585944      FALSE
#> 70 70     16     18  2.13324117 0.7384993  4.00101566 1.5358175      FALSE
#> 71 71     17     17  2.76643368 0.9204042  2.76643368 0.9204042      FALSE
#> 72 72      6      7 -0.79407316 0.5164849 -0.54647281 0.5071109      FALSE
#> 73 73     12     14  0.72136345 0.5361040  1.32134359 0.5929126      FALSE
#> 74 74      4      9 -1.33332113 0.5641476 -0.05692244 0.5055699      FALSE
#> 75 75      1      5 -2.61813570 0.8735700 -1.05265787 0.5342416      FALSE
#> 76 76     11     15  0.45204624 0.5211334  1.68336335 0.6470375      FALSE
#> 77 77      9     12 -0.05692244 0.5055699  0.72136345 0.5361040      FALSE
#> 78 78      9     12 -0.05692244 0.5055699  0.72136345 0.5361040      FALSE
#> 79 79     10     14  0.19368859 0.5113013  1.32134359 0.5929126      FALSE
#> 80 80      6      7 -0.79407316 0.5164849 -0.54647281 0.5071109      FALSE
#>    extreme_t2     change   se_diff        rci      p_value crit_lower
#> 1       FALSE  0.7501926 0.7344092  1.0214914 0.3060315291    -1.7412
#> 2       FALSE  0.7782859 0.7368910  1.0561750 0.2630420758    -1.7412
#> 3       FALSE  0.9877617 0.7267638  1.3591234 0.1544650281    -1.7412
#> 4       FALSE  2.0450702 1.0651532  1.9199775 0.0316457520    -1.7412
#> 5       FALSE  1.4896748 0.8246736  1.8063810 0.0412876851    -1.7412
#> 6       FALSE  1.2463465 0.7394885  1.6854171 0.0682679217    -1.7412
#> 7       FALSE  0.5999801 0.7993452  0.7505895 0.4020234001    -1.7412
#> 8       FALSE -0.6759066 0.8548013 -0.7907178 0.3809688244    -1.7412
#> 9       FALSE  0.5999801 0.7993452  0.7505895 0.4020234001    -1.7412
#> 10      FALSE  0.9877617 0.7267638  1.3591234 0.1544650281    -1.7412
#> 11      FALSE  0.8118976 0.9470620  0.8572803 0.3567022977    -1.7412
#> 12      FALSE  1.4896748 0.8246736  1.8063810 0.0412876851    -1.7412
#> 13      FALSE  0.0000000 0.7171631  0.0000000 1.0000000000    -1.7412
#> 14      FALSE  1.7172011 1.6288557  1.0542377 0.2676546678    -1.7412
#> 15      FALSE -0.6331925 1.1800530 -0.5365797 0.5744960963    -1.7412
#> 16      FALSE  1.5539295 0.7544463  2.0596955 0.0205922019    -1.7412
#> 17      FALSE  1.6238088 0.7781311  2.0868061 0.0181212067    -1.7412
#> 18      FALSE  0.0000000 1.0443957  0.0000000 1.0000000000    -1.7412
#> 19      FALSE  2.0450702 1.0651532  1.9199775 0.0316457520    -1.7412
#> 20      FALSE  2.2298362 0.8220821  2.7124251 0.0015726735    -1.7412
#> 21      FALSE  1.0830703 1.1250785  0.9626621 0.3338917140    -1.7412
#> 22      FALSE -1.0830703 1.1250785 -0.9626621 0.3338917140    -1.7412
#> 23      FALSE  0.2506110 0.7190480  0.3485317 0.7714279194    -1.7412
#> 24      FALSE  0.0000000 1.2354146  0.0000000 1.0000000000    -1.7412
#> 25      FALSE  2.3143874 1.0576975  2.1881373 0.0118718001    -1.7412
#> 26      FALSE  0.5999801 0.7993452  0.7505895 0.4020234001    -1.7412
#> 27      FALSE  1.8015299 0.7607788  2.3680076 0.0075684545    -1.7412
#> 28      FALSE  0.0000000 0.8676531  0.0000000 1.0000000000    -1.7412
#> 29      FALSE  1.8678164 0.7801966  2.3940329 0.0065276664    -1.7412
#> 30      FALSE  0.6759066 0.8548013  0.7907178 0.3809688244    -1.7412
#> 31       TRUE  1.8677745 1.7041469  1.0960173 0.2315434611    -1.7412
#> 32      FALSE  0.2506110 0.7190480  0.3485317 0.7714279194    -1.7412
#> 33      FALSE  0.5999801 0.7993452  0.7505895 0.4020234001    -1.7412
#> 34      FALSE  1.2678363 0.7379491  1.7180538 0.0580826487    -1.7412
#> 35      FALSE -0.3978033 0.9302615 -0.4276252 0.6146606174    -1.7412
#> 36      FALSE  1.6238088 0.7781311  2.0868061 0.0181212067    -1.7412
#> 37      FALSE  1.5047041 0.7463204  2.0161636 0.0259764515    -1.7412
#> 38      FALSE  0.4961538 0.7178940  0.6911241 0.4882863699    -1.7412
#> 39      FALSE  0.0000000 1.0443957  0.0000000 1.0000000000    -1.7412
#> 40      FALSE  1.4118777 0.9125725  1.5471403 0.0857589874    -1.7412
#> 41      FALSE -0.2583577 0.7300747 -0.3538784 0.7494164884    -1.7412
#> 42      FALSE  0.0000000 0.7126608  0.0000000 1.0000000000    -1.7412
#> 43      FALSE  0.6759066 0.8548013  0.7907178 0.3809688244    -1.7412
#> 44      FALSE -0.4498778 0.9818547 -0.4581918 0.6070605206    -1.7412
#> 45      FALSE  1.0308559 0.7564425  1.3627684 0.1469363378    -1.7412
#> 46      FALSE  0.4916079 0.7215950  0.6812796 0.5423847296    -1.7412
#> 47      FALSE  1.2313171 0.8308054  1.4820764 0.0972227850    -1.7412
#> 48      FALSE  0.8118976 0.9470620  0.8572803 0.3567022977    -1.7412
#> 49      FALSE -0.2860933 0.7742320 -0.3695188 0.6956156963    -1.7412
#> 50      FALSE  0.3620198 0.8776120  0.4125055 0.6379025959    -1.7412
#> 51      FALSE  0.0000000 0.7369940  0.0000000 1.0000000000    -1.7412
#> 52      FALSE  2.1065610 0.8049789  2.6169145 0.0031752904    -1.7412
#> 53      FALSE  2.6797140 0.8958474  2.9912616 0.0004378706    -1.7412
#> 54      FALSE  2.3758782 0.8147505  2.9160805 0.0008587425    -1.7412
#> 55      FALSE  0.0000000 1.3016482  0.0000000 1.0000000000    -1.7412
#> 56      FALSE -1.2345820 1.7904970 -0.6895192 0.5048526219    -1.7412
#> 57      FALSE  0.7371507 0.7227431  1.0199346 0.3167815724    -1.7412
#> 58      FALSE  0.2440076 0.7149155  0.3413097 0.8303248175    -1.7412
#> 59      FALSE  0.3138869 0.8145999  0.3853264 0.6695244154    -1.7412
#> 60      FALSE  0.9985191 0.7271461  1.3732027 0.1428504543    -1.7412
#> 61       TRUE  1.2345820 1.7904970  0.6895192 0.5048526219    -1.7412
#> 62      FALSE  0.5276749 0.7408350  0.7122704 0.4468168495    -1.7412
#> 63       TRUE  2.9935589 1.6342470  1.8317665 0.0367236311    -1.7412
#> 64      FALSE  2.1901636 0.8949760  2.4471757 0.0041633503    -1.7412
#> 65      FALSE  1.2848146 1.0398976  1.2355202 0.1760409997    -1.7412
#> 66      FALSE  2.6797140 0.8958474  2.9912616 0.0004378706    -1.7412
#> 67      FALSE  0.0000000 1.2354146  0.0000000 1.0000000000    -1.7412
#> 68       TRUE  1.8677745 1.7041469  1.0960173 0.2315434611    -1.7412
#> 69      FALSE  0.5554105 0.7639423  0.7270320 0.4269539073    -1.7412
#> 70       TRUE  1.8677745 1.7041469  1.0960173 0.2315434611    -1.7412
#> 71      FALSE  0.0000000 1.3016482  0.0000000 1.0000000000    -1.7412
#> 72      FALSE  0.2476003 0.7238218  0.3420736 0.8127906708    -1.7412
#> 73      FALSE  0.5999801 0.7993452  0.7505895 0.4020234001    -1.7412
#> 74      FALSE  1.2763987 0.7575378  1.6849307 0.0712117584    -1.7412
#> 75      FALSE  1.5654778 1.0239818  1.5288141 0.0887072204    -1.7412
#> 76      FALSE  1.2313171 0.8308054  1.4820764 0.0972227850    -1.7412
#> 77      FALSE  0.7782859 0.7368910  1.0561750 0.2630420758    -1.7412
#> 78      FALSE  0.7782859 0.7368910  1.0561750 0.2630420758    -1.7412
#> 79      FALSE  1.1276550 0.7829268  1.4403070 0.1078702298    -1.7412
#> 80      FALSE  0.2476003 0.7238218  0.3420736 0.8127906708    -1.7412
#>    crit_upper  change_class retest_sd_tip
#> 1      1.7412 none detected            NA
#> 2      1.7412 none detected            NA
#> 3      1.7412 none detected            NA
#> 4      1.7412      increase     0.3499579
#> 5      1.7412      increase     0.1610447
#> 6      1.7412 none detected            NA
#> 7      1.7412 none detected            NA
#> 8      1.7412 none detected            NA
#> 9      1.7412 none detected            NA
#> 10     1.7412 none detected            NA
#> 11     1.7412 none detected            NA
#> 12     1.7412      increase     0.1610447
#> 13     1.7412 none detected            NA
#> 14     1.7412 none detected            NA
#> 15     1.7412 none detected            NA
#> 16     1.7412      increase     0.3371007
#> 17     1.7412      increase     0.3634683
#> 18     1.7412 none detected            NA
#> 19     1.7412      increase     0.3499579
#> 20     1.7412      increase     0.6943341
#> 21     1.7412 none detected            NA
#> 22     1.7412 none detected            NA
#> 23     1.7412 none detected            NA
#> 24     1.7412 none detected            NA
#> 25     1.7412      increase     0.5692212
#> 26     1.7412 none detected            NA
#> 27     1.7412      increase     0.4958395
#> 28     1.7412 none detected            NA
#> 29     1.7412      increase     0.5205850
#> 30     1.7412 none detected            NA
#> 31     1.7412 none detected            NA
#> 32     1.7412 none detected            NA
#> 33     1.7412 none detected            NA
#> 34     1.7412 none detected            NA
#> 35     1.7412 none detected            NA
#> 36     1.7412      increase     0.3634683
#> 37     1.7412      increase     0.3080643
#> 38     1.7412 none detected            NA
#> 39     1.7412 none detected            NA
#> 40     1.7412 none detected            NA
#> 41     1.7412 none detected            NA
#> 42     1.7412 none detected            NA
#> 43     1.7412 none detected            NA
#> 44     1.7412 none detected            NA
#> 45     1.7412 none detected            NA
#> 46     1.7412 none detected            NA
#> 47     1.7412 none detected            NA
#> 48     1.7412 none detected            NA
#> 49     1.7412 none detected            NA
#> 50     1.7412 none detected            NA
#> 51     1.7412 none detected            NA
#> 52     1.7412      increase     0.6386334
#> 53     1.7412      increase     0.8848707
#> 54     1.7412      increase     0.7739702
#> 55     1.7412 none detected            NA
#> 56     1.7412 none detected            NA
#> 57     1.7412 none detected            NA
#> 58     1.7412 none detected            NA
#> 59     1.7412 none detected            NA
#> 60     1.7412 none detected            NA
#> 61     1.7412 none detected            NA
#> 62     1.7412 none detected            NA
#> 63     1.7412      increase     0.3775317
#> 64     1.7412      increase     0.6249793
#> 65     1.7412 none detected            NA
#> 66     1.7412      increase     0.8848707
#> 67     1.7412 none detected            NA
#> 68     1.7412 none detected            NA
#> 69     1.7412 none detected            NA
#> 70     1.7412 none detected            NA
#> 71     1.7412 none detected            NA
#> 72     1.7412 none detected            NA
#> 73     1.7412 none detected            NA
#> 74     1.7412 none detected            NA
#> 75     1.7412 none detected            NA
#> 76     1.7412 none detected            NA
#> 77     1.7412 none detected            NA
#> 78     1.7412 none detected            NA
#> 79     1.7412 none detected            NA
#> 80     1.7412 none detected            NA

# A single respondent needs an external calibration
names(thr) <- paste0("I", 1:6)
RMpersonChange(t1[1, ], t2[1, ], item_params = thr)
#>   id sum_t1 sum_t2   theta_t1     se_t1  theta_t2     se_t2 extreme_t1
#> 1  1      5      8 -0.9554345 0.5288866 -0.228083 0.4890212      FALSE
#>   extreme_t2    change   se_diff     rci   p_value crit_lower crit_upper
#> 1      FALSE 0.7273516 0.7203213 1.00976 0.3107644  -1.791785   1.791785
#>    change_class retest_sd_tip
#> 1 none detected            NA

# Monte Carlo alternative, for comparison
RMpersonChange(t1, t2, critical = "simulate", sim_iter = 200,
               parallel = FALSE, seed = 1)
#>    id sum_t1 sum_t2    theta_t1     se_t1    theta_t2     se_t2 extreme_t1
#> 1   1      5      8 -1.05265787 0.5342416 -0.30246523 0.5039273      FALSE
#> 2   2      9     12 -0.05692244 0.5055699  0.72136345 0.5361040      FALSE
#> 3   3      6     10 -0.79407316 0.5164849  0.19368859 0.5113013      FALSE
#> 4   4     12     17  0.72136345 0.5361040  2.76643368 0.9204042      FALSE
#> 5   5     10     15  0.19368859 0.5113013  1.68336335 0.6470375      FALSE
#> 6   6      5     10 -1.05265787 0.5342416  0.19368859 0.5113013      FALSE
#> 7   7     12     14  0.72136345 0.5361040  1.32134359 0.5929126      FALSE
#> 8   8     15     13  1.68336335 0.6470375  1.00745673 0.5585944      FALSE
#> 9   9     12     14  0.72136345 0.5361040  1.32134359 0.5929126      FALSE
#> 10 10      6     10 -0.79407316 0.5164849  0.19368859 0.5113013      FALSE
#> 11 11     14     16  1.32134359 0.5929126  2.13324117 0.7384993      FALSE
#> 12 12     10     15  0.19368859 0.5113013  1.68336335 0.6470375      FALSE
#> 13 13      7      7 -0.54647281 0.5071109 -0.54647281 0.5071109      FALSE
#> 14 14      0      2 -3.76951909 1.4711205 -2.05231802 0.6992678       TRUE
#> 15 15     17     16  2.76643368 0.9204042  2.13324117 0.7384993      FALSE
#> 16 16      7     13 -0.54647281 0.5071109  1.00745673 0.5585944      FALSE
#> 17 17      8     14 -0.30246523 0.5039273  1.32134359 0.5929126      FALSE
#> 18 18     16     16  2.13324117 0.7384993  2.13324117 0.7384993      FALSE
#> 19 19     12     17  0.72136345 0.5361040  2.76643368 0.9204042      FALSE
#> 20 20      7     15 -0.54647281 0.5071109  1.68336335 0.6470375      FALSE
#> 21 21     15     17  1.68336335 0.6470375  2.76643368 0.9204042      FALSE
#> 22 22     17     15  2.76643368 0.9204042  1.68336335 0.6470375      FALSE
#> 23 23      9     10 -0.05692244 0.5055699  0.19368859 0.5113013      FALSE
#> 24 24      1      1 -2.61813570 0.8735700 -2.61813570 0.8735700      FALSE
#> 25 25     11     17  0.45204624 0.5211334  2.76643368 0.9204042      FALSE
#> 26 26     12     14  0.72136345 0.5361040  1.32134359 0.5929126      FALSE
#> 27 27      6     13 -0.79407316 0.5164849  1.00745673 0.5585944      FALSE
#> 28 28      3      3 -1.65451471 0.6135234 -1.65451471 0.6135234      FALSE
#> 29 29      7     14 -0.54647281 0.5071109  1.32134359 0.5929126      FALSE
#> 30 30     13     15  1.00745673 0.5585944  1.68336335 0.6470375      FALSE
#> 31 31     16     18  2.13324117 0.7384993  4.00101566 1.5358175      FALSE
#> 32 32      9     10 -0.05692244 0.5055699  0.19368859 0.5113013      FALSE
#> 33 33     12     14  0.72136345 0.5361040  1.32134359 0.5929126      FALSE
#> 34 34      7     12 -0.54647281 0.5071109  0.72136345 0.5361040      FALSE
#> 35 35      3      2 -1.65451471 0.6135234 -2.05231802 0.6992678      FALSE
#> 36 36      8     14 -0.30246523 0.5039273  1.32134359 0.5929126      FALSE
#> 37 37      5     11 -1.05265787 0.5342416  0.45204624 0.5211334      FALSE
#> 38 38      8     10 -0.30246523 0.5039273  0.19368859 0.5113013      FALSE
#> 39 39     16     16  2.13324117 0.7384993  2.13324117 0.7384993      FALSE
#> 40 40     12     16  0.72136345 0.5361040  2.13324117 0.7384993      FALSE
#> 41 41     11     10  0.45204624 0.5211334  0.19368859 0.5113013      FALSE
#> 42 42      8      8 -0.30246523 0.5039273 -0.30246523 0.5039273      FALSE
#> 43 43     13     15  1.00745673 0.5585944  1.68336335 0.6470375      FALSE
#> 44 44     16     15  2.13324117 0.7384993  1.68336335 0.6470375      FALSE
#> 45 45      4      8 -1.33332113 0.5641476 -0.30246523 0.5039273      FALSE
#> 46 46      6      8 -0.79407316 0.5164849 -0.30246523 0.5039273      FALSE
#> 47 47     11     15  0.45204624 0.5211334  1.68336335 0.6470375      FALSE
#> 48 48     14     16  1.32134359 0.5929126  2.13324117 0.7384993      FALSE
#> 49 49     13     12  1.00745673 0.5585944  0.72136345 0.5361040      FALSE
#> 50 50     14     15  1.32134359 0.5929126  1.68336335 0.6470375      FALSE
#> 51 51     11     11  0.45204624 0.5211334  0.45204624 0.5211334      FALSE
#> 52 52      3     11 -1.65451471 0.6135234  0.45204624 0.5211334      FALSE
#> 53 53      7     16 -0.54647281 0.5071109  2.13324117 0.7384993      FALSE
#> 54 54      3     12 -1.65451471 0.6135234  0.72136345 0.5361040      FALSE
#> 55 55     17     17  2.76643368 0.9204042  2.76643368 0.9204042      FALSE
#> 56 56     18     17  4.00101566 1.5358175  2.76643368 0.9204042       TRUE
#> 57 57      6      9 -0.79407316 0.5164849 -0.05692244 0.5055699      FALSE
#> 58 58      7      8 -0.54647281 0.5071109 -0.30246523 0.5039273      FALSE
#> 59 59     13     14  1.00745673 0.5585944  1.32134359 0.5929126      FALSE
#> 60 60      7     11 -0.54647281 0.5071109  0.45204624 0.5211334      FALSE
#> 61 61     17     18  2.76643368 0.9204042  4.00101566 1.5358175      FALSE
#> 62 62     10     12  0.19368859 0.5113013  0.72136345 0.5361040      FALSE
#> 63 63     13     18  1.00745673 0.5585944  4.00101566 1.5358175      FALSE
#> 64 64      9     16 -0.05692244 0.5055699  2.13324117 0.7384993      FALSE
#> 65 65      1      4 -2.61813570 0.8735700 -1.33332113 0.5641476      FALSE
#> 66 66      7     16 -0.54647281 0.5071109  2.13324117 0.7384993      FALSE
#> 67 67      1      1 -2.61813570 0.8735700 -2.61813570 0.8735700      FALSE
#> 68 68     16     18  2.13324117 0.7384993  4.00101566 1.5358175      FALSE
#> 69 69     11     13  0.45204624 0.5211334  1.00745673 0.5585944      FALSE
#> 70 70     16     18  2.13324117 0.7384993  4.00101566 1.5358175      FALSE
#> 71 71     17     17  2.76643368 0.9204042  2.76643368 0.9204042      FALSE
#> 72 72      6      7 -0.79407316 0.5164849 -0.54647281 0.5071109      FALSE
#> 73 73     12     14  0.72136345 0.5361040  1.32134359 0.5929126      FALSE
#> 74 74      4      9 -1.33332113 0.5641476 -0.05692244 0.5055699      FALSE
#> 75 75      1      5 -2.61813570 0.8735700 -1.05265787 0.5342416      FALSE
#> 76 76     11     15  0.45204624 0.5211334  1.68336335 0.6470375      FALSE
#> 77 77      9     12 -0.05692244 0.5055699  0.72136345 0.5361040      FALSE
#> 78 78      9     12 -0.05692244 0.5055699  0.72136345 0.5361040      FALSE
#> 79 79     10     14  0.19368859 0.5113013  1.32134359 0.5929126      FALSE
#> 80 80      6      7 -0.79407316 0.5164849 -0.54647281 0.5071109      FALSE
#>    extreme_t2     change   se_diff        rci      p_value crit_lower
#> 1       FALSE  0.7501926 0.7344092  1.0214914 0.3078557590    -1.7412
#> 2       FALSE  0.7782859 0.7368910  1.0561750 0.2647959503    -1.7412
#> 3       FALSE  0.9877617 0.7267638  1.3591234 0.1576776451    -1.7412
#> 4       FALSE  2.0450702 1.0651532  1.9199775 0.0318105118    -1.7412
#> 5       FALSE  1.4896748 0.8246736  1.8063810 0.0428098244    -1.7412
#> 6       FALSE  1.2463465 0.7394885  1.6854171 0.0682457346    -1.7412
#> 7       FALSE  0.5999801 0.7993452  0.7505895 0.4055371539    -1.7412
#> 8       FALSE -0.6759066 0.8548013 -0.7907178 0.3848509468    -1.7412
#> 9       FALSE  0.5999801 0.7993452  0.7505895 0.4055371539    -1.7412
#> 10      FALSE  0.9877617 0.7267638  1.3591234 0.1576776451    -1.7412
#> 11      FALSE  0.8118976 0.9470620  0.8572803 0.3611024311    -1.7412
#> 12      FALSE  1.4896748 0.8246736  1.8063810 0.0428098244    -1.7412
#> 13      FALSE  0.0000000 0.7171631  0.0000000 1.0000000000    -1.7412
#> 14      FALSE  1.7172011 1.6288557  1.0542377 0.2693581651    -1.7412
#> 15      FALSE -0.6331925 1.1800530 -0.5365797 0.5769014437    -1.7412
#> 16      FALSE  1.5539295 0.7544463  2.0596955 0.0205612149    -1.7412
#> 17      FALSE  1.6238088 0.7781311  2.0868061 0.0177488907    -1.7412
#> 18      FALSE  0.0000000 1.0443957  0.0000000 1.0000000000    -1.7412
#> 19      FALSE  2.0450702 1.0651532  1.9199775 0.0318105118    -1.7412
#> 20      FALSE  2.2298362 0.8220821  2.7124251 0.0019373789    -1.7412
#> 21      FALSE  1.0830703 1.1250785  0.9626621 0.3381663646    -1.7412
#> 22      FALSE -1.0830703 1.1250785 -0.9626621 0.3381663646    -1.7412
#> 23      FALSE  0.2506110 0.7190480  0.3485317 0.7727642022    -1.7412
#> 24      FALSE  0.0000000 1.2354146  0.0000000 1.0000000000    -1.7412
#> 25      FALSE  2.3143874 1.0576975  2.1881373 0.0108743204    -1.7412
#> 26      FALSE  0.5999801 0.7993452  0.7505895 0.4055371539    -1.7412
#> 27      FALSE  1.8015299 0.7607788  2.3680076 0.0069370664    -1.7412
#> 28      FALSE  0.0000000 0.8676531  0.0000000 1.0000000000    -1.7412
#> 29      FALSE  1.8678164 0.7801966  2.3940329 0.0060621211    -1.7412
#> 30      FALSE  0.6759066 0.8548013  0.7907178 0.3848509468    -1.7412
#> 31       TRUE  1.8677745 1.7041469  1.0960173 0.2341728642    -1.7412
#> 32      FALSE  0.2506110 0.7190480  0.3485317 0.7727642022    -1.7412
#> 33      FALSE  0.5999801 0.7993452  0.7505895 0.4055371539    -1.7412
#> 34      FALSE  1.2678363 0.7379491  1.7180538 0.0593712893    -1.7412
#> 35      FALSE -0.3978033 0.9302615 -0.4276252 0.6182738579    -1.7412
#> 36      FALSE  1.6238088 0.7781311  2.0868061 0.0177488907    -1.7412
#> 37      FALSE  1.5047041 0.7463204  2.0161636 0.0266858321    -1.7412
#> 38      FALSE  0.4961538 0.7178940  0.6911241 0.4925317168    -1.7412
#> 39      FALSE  0.0000000 1.0443957  0.0000000 1.0000000000    -1.7412
#> 40      FALSE  1.4118777 0.9125725  1.5471403 0.0863071058    -1.7412
#> 41      FALSE -0.2583577 0.7300747 -0.3538784 0.7506405850    -1.7412
#> 42      FALSE  0.0000000 0.7126608  0.0000000 1.0000000000    -1.7412
#> 43      FALSE  0.6759066 0.8548013  0.7907178 0.3848509468    -1.7412
#> 44      FALSE -0.4498778 0.9818547 -0.4581918 0.6089619399    -1.7412
#> 45      FALSE  1.0308559 0.7564425  1.3627684 0.1497406412    -1.7412
#> 46      FALSE  0.4916079 0.7215950  0.6812796 0.5455909006    -1.7412
#> 47      FALSE  1.2313171 0.8308054  1.4820764 0.0989938129    -1.7412
#> 48      FALSE  0.8118976 0.9470620  0.8572803 0.3611024311    -1.7412
#> 49      FALSE -0.2860933 0.7742320 -0.3695188 0.6975189051    -1.7412
#> 50      FALSE  0.3620198 0.8776120  0.4125055 0.6419598775    -1.7412
#> 51      FALSE  0.0000000 0.7369940  0.0000000 1.0000000000    -1.7412
#> 52      FALSE  2.1065610 0.8049789  2.6169145 0.0032497969    -1.7412
#> 53      FALSE  2.6797140 0.8958474  2.9912616 0.0005624648    -1.7412
#> 54      FALSE  2.3758782 0.8147505  2.9160805 0.0008124492    -1.7412
#> 55      FALSE  0.0000000 1.3016482  0.0000000 1.0000000000    -1.7412
#> 56      FALSE -1.2345820 1.7904970 -0.6895192 0.5097806387    -1.7412
#> 57      FALSE  0.7371507 0.7227431  1.0199346 0.3190425598    -1.7412
#> 58      FALSE  0.2440076 0.7149155  0.3413097 0.8304480970    -1.7412
#> 59      FALSE  0.3138869 0.8145999  0.3853264 0.6736453972    -1.7412
#> 60      FALSE  0.9985191 0.7271461  1.3732027 0.1453659146    -1.7412
#> 61       TRUE  1.2345820 1.7904970  0.6895192 0.5097806387    -1.7412
#> 62      FALSE  0.5276749 0.7408350  0.7122704 0.4515342791    -1.7412
#> 63       TRUE  2.9935589 1.6342470  1.8317665 0.0376226486    -1.7412
#> 64      FALSE  2.1901636 0.8949760  2.4471757 0.0039372539    -1.7412
#> 65      FALSE  1.2848146 1.0398976  1.2355202 0.1793012937    -1.7412
#> 66      FALSE  2.6797140 0.8958474  2.9912616 0.0005624648    -1.7412
#> 67      FALSE  0.0000000 1.2354146  0.0000000 1.0000000000    -1.7412
#> 68       TRUE  1.8677745 1.7041469  1.0960173 0.2341728642    -1.7412
#> 69      FALSE  0.5554105 0.7639423  0.7270320 0.4312855447    -1.7412
#> 70       TRUE  1.8677745 1.7041469  1.0960173 0.2341728642    -1.7412
#> 71      FALSE  0.0000000 1.3016482  0.0000000 1.0000000000    -1.7412
#> 72      FALSE  0.2476003 0.7238218  0.3420736 0.8125742141    -1.7412
#> 73      FALSE  0.5999801 0.7993452  0.7505895 0.4055371539    -1.7412
#> 74      FALSE  1.2763987 0.7575378  1.6849307 0.0716830198    -1.7412
#> 75      FALSE  1.5654778 1.0239818  1.5288141 0.0898693832    -1.7412
#> 76      FALSE  1.2313171 0.8308054  1.4820764 0.0989938129    -1.7412
#> 77      FALSE  0.7782859 0.7368910  1.0561750 0.2647959503    -1.7412
#> 78      FALSE  0.7782859 0.7368910  1.0561750 0.2647959503    -1.7412
#> 79      FALSE  1.1276550 0.7829268  1.4403070 0.1103056059    -1.7412
#> 80      FALSE  0.2476003 0.7238218  0.3420736 0.8125742141    -1.7412
#>    crit_upper  change_class retest_sd_tip
#> 1    1.743296 none detected            NA
#> 2    1.743296 none detected            NA
#> 3    1.743296 none detected            NA
#> 4    1.743296      increase     0.3475814
#> 5    1.743296      increase     0.1582903
#> 6    1.743296 none detected            NA
#> 7    1.743296 none detected            NA
#> 8    1.743296 none detected            NA
#> 9    1.743296 none detected            NA
#> 10   1.743296 none detected            NA
#> 11   1.743296 none detected            NA
#> 12   1.743296      increase     0.1582903
#> 13   1.743296 none detected            NA
#> 14   1.743296 none detected            NA
#> 15   1.743296 none detected            NA
#> 16   1.743296      increase     0.3356781
#> 17   1.743296      increase     0.3620277
#> 18   1.743296 none detected            NA
#> 19   1.743296      increase     0.3475814
#> 20   1.743296      increase     0.6929135
#> 21   1.743296 none detected            NA
#> 22   1.743296 none detected            NA
#> 23   1.743296 none detected            NA
#> 24   1.743296 none detected            NA
#> 25   1.743296      increase     0.5673533
#> 26   1.743296 none detected            NA
#> 27   1.743296      increase     0.4945406
#> 28   1.743296 none detected            NA
#> 29   1.743296      increase     0.5192552
#> 30   1.743296 none detected            NA
#> 31   1.743296 none detected            NA
#> 32   1.743296 none detected            NA
#> 33   1.743296 none detected            NA
#> 34   1.743296 none detected            NA
#> 35   1.743296 none detected            NA
#> 36   1.743296      increase     0.3620277
#> 37   1.743296      increase     0.3066043
#> 38   1.743296 none detected            NA
#> 39   1.743296 none detected            NA
#> 40   1.743296 none detected            NA
#> 41   1.743296 none detected            NA
#> 42   1.743296 none detected            NA
#> 43   1.743296 none detected            NA
#> 44   1.743296 none detected            NA
#> 45   1.743296 none detected            NA
#> 46   1.743296 none detected            NA
#> 47   1.743296 none detected            NA
#> 48   1.743296 none detected            NA
#> 49   1.743296 none detected            NA
#> 50   1.743296 none detected            NA
#> 51   1.743296 none detected            NA
#> 52   1.743296      increase     0.6372548
#> 53   1.743296      increase     0.8832610
#> 54   1.743296      increase     0.7725234
#> 55   1.743296 none detected            NA
#> 56   1.743296 none detected            NA
#> 57   1.743296 none detected            NA
#> 58   1.743296 none detected            NA
#> 59   1.743296 none detected            NA
#> 60   1.743296 none detected            NA
#> 61   1.743296 none detected            NA
#> 62   1.743296 none detected            NA
#> 63   1.743296      increase     0.3727979
#> 64   1.743296      increase     0.6234564
#> 65   1.743296 none detected            NA
#> 66   1.743296      increase     0.8832610
#> 67   1.743296 none detected            NA
#> 68   1.743296 none detected            NA
#> 69   1.743296 none detected            NA
#> 70   1.743296 none detected            NA
#> 71   1.743296 none detected            NA
#> 72   1.743296 none detected            NA
#> 73   1.743296 none detected            NA
#> 74   1.743296 none detected            NA
#> 75   1.743296 none detected            NA
#> 76   1.743296 none detected            NA
#> 77   1.743296 none detected            NA
#> 78   1.743296 none detected            NA
#> 79   1.743296 none detected            NA
#> 80   1.743296 none detected            NA
# }
```
