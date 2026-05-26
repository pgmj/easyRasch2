# Martin-Lof Test of Unidimensionality

Likelihood-ratio test of unidimensionality against an *a priori*
specified multidimensional alternative, generalised to polytomous Rasch
/ partial credit models (Christensen, Bjorner, Kreiner, & Petersen,
2002). The p-value is obtained by parametric-bootstrap (Monte Carlo)
sampling under the unidimensional null, following Christensen & Kreiner
(2007), because the asymptotic chi-square approximation is biased toward
conservatism for realistic sample sizes – especially with polytomous
items, where the degrees of freedom can be very large.

## Usage

``` r
RMdimMartinLof(
  data,
  partition,
  iterations = 1000L,
  stopping = c("none", "sequential"),
  h = 50L,
  alpha = 0.05,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL
)
```

## Arguments

- data:

  A data.frame or matrix of item responses (0-based, non-negative
  integers). Rows with any `NA` are dropped.

- partition:

  The hypothesised partition of items into subscales. One of:

  - a **list** of column-name or column-index vectors, e.g.
    `list(c("I1","I2","I3"), c("I4","I5","I6"))`;

  - a **vector** of length `ncol(data)` indicating each item's subscale
    (factor, character, or integer), e.g. `c(1,1,1,2,2,2)`. Each
    subscale must contain at least two items. Subscales must not
    overlap; items not assigned to any subscale are dropped with a
    warning.

- iterations:

  Integer. Maximum number of Monte Carlo iterations (default `1000`).

- stopping:

  Character. `"none"` (default) runs all iterations. `"sequential"` uses
  Besag & Clifford's (1991) sequential rule: stop as soon as `h`
  simulated statistics have exceeded the observed value. The sequential
  strategy substantially reduces compute time when H0 holds but cannot
  be parallelised.

- h:

  Integer. Sequential-stopping count threshold (default `50`). Ignored
  when `stopping = "none"`.

- alpha:

  Numeric in (0, 1). Nominal significance level used only for the
  `rejected` flag in the result; default `0.05`.

- parallel:

  Logical. Use parallel processing via `mirai` (default `TRUE`). Ignored
  when `stopping = "sequential"`.

- n_cores:

  Integer or `NULL`. Number of parallel workers. When `NULL`,
  `getOption("mc.cores")` is checked first; if neither is set, falls
  back to sequential with a warning.

- verbose:

  Logical. Show a progress bar (default `FALSE`).

- seed:

  Integer or `NULL`. Random seed for reproducibility.

## Value

A list with components:

- `T_obs`:

  Observed Martin-Lof likelihood-ratio statistic.

- `p_value`:

  Monte Carlo p-value with `(n_exceed + 1) / (n + 1)` correction.

- `actual_iterations`:

  Number of successful MC iterations completed.

- `rejected`:

  Logical: is `p_value < alpha`?

- `partition`:

  Normalised partition (list of integer indices).

- `n_subscales`:

  Number of subscales.

- `is_polytomous`:

  Whether a PCM was fitted.

- `sample_n`:

  Number of complete cases analysed.

- `n_items`:

  Number of items.

- `stopping`:

  The stopping strategy used.

- `h`:

  The sequential-stopping count, or `NA` for `stopping = "none"`.

- `T_rep`:

  Numeric vector of successful MC test statistics.

- `wle_scores`:

  data.frame with one row per person and one column per subscale
  (`subscale_1_wle`, ..., `subscale_D_wle`), giving Warm's Weighted
  Likelihood Estimate of theta from a CML fit on each subscale alone.
  Persons whose subscore equals the minimum or maximum on a subscale
  produce non-finite WLEs (`Inf` / `-Inf`) and are excluded from
  `wle_correlation` pairwise.

- `wle_correlation`:

  data.frame of pairwise Pearson correlations between subscale WLEs,
  with columns `subscale_a`, `subscale_b`, `r`, `ci_lower`, `ci_upper`
  (95% CI from
  [`stats::cor.test`](https://rdrr.io/r/stats/cor.test.html)),
  `p_value`, and `n` (number of persons with finite WLEs on both
  subscales). One row per pair; for D = 2, a single row. Useful as an
  effect-size companion to `p_value` – a rejected test with `r` near 1
  indicates a small effect; `r` clearly below 1 indicates substantive
  multidimensionality.

## Details

This is **not** a routine screening tool. The test requires an *a
priori* partition of items into subscales; using it post-hoc on, e.g.,
the partition suggested by
[`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/reference/RMdimResidualPCA.md)'s
PC1 sign would inflate the Type-I error rate. Both source papers state
this explicitly.

**Test statistic.** With items partitioned into D subscales, total score
\\t\\ and subscores \\(t_1, \ldots, t_D)\\ (Christensen et al. 2002, eq.
22): \$\$T = 2\Bigl\[\sum\_{t_1, \ldots, t_D} n\_{t_1, \ldots,
t_D}\log(n\_{t_1, \ldots, t_D}/N) - \sum_t n_t\log(n_t/N) -
\ell_C(\hat{\epsilon}) + \sum_d \ell_C(\hat{\epsilon}^{(d)})\Bigr\]\$\$
where \\\ell_C\\ is the conditional log-likelihood and the
\\\hat{\epsilon}^{(d)}\\ are CML estimates on the d-th subscale alone.
CML fits use
[`psychotools::raschmodel()`](https://rdrr.io/pkg/psychotools/man/raschmodel.html)
(RM) or
[`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)
(PCM) for speed.

**Monte Carlo sampling under H0.** Following Christensen & Kreiner
(2007): (a) sample N total scores from the empirical score distribution
\\n_t/N\\; (b) for each sampled score, sample an item-response vector
from the conditional distribution \\p(x \mid t, \hat{\epsilon})\\ given
by eq. 4 of the paper. For dichotomous data the fast algorithm of
Christensen & Kreiner (2007, p. 23) is used (sample without replacement
weighted by item easinesses). For polytomous data the recursive
\\\gamma\\-function approach is used, with each item's response sampled
conditional on the remaining items' joint score distribution (computed
via
[`psychotools::elementary_symmetric_functions()`](https://rdrr.io/pkg/psychotools/man/elementary_symmetric_functions.html)).

Iterations that fail (e.g., simulated dataset has an empty category for
an item) are silently dropped.

Item parameters are estimated once on the observed data and held fixed
across MC iterations. Christensen & Kreiner (2007) use the extended
likelihood function (Tjur, 1982) with the empirical score distribution
as a non-parametric estimate of the latent distribution, so no
distributional assumption about \\\theta\\ is needed.

## References

Christensen, K. B., Bjorner, J. B., Kreiner, S., & Petersen, J. H.
(2002). Testing unidimensionality in polytomous Rasch models.
*Psychometrika, 67*(4), 563-574.
[doi:10.1007/BF02295132](https://doi.org/10.1007/BF02295132)

Christensen, K. B., & Kreiner, S. (2007). A Monte Carlo approach to
unidimensionality testing in polytomous Rasch models. *Applied
Psychological Measurement, 31*(1), 20-30.
[doi:10.1177/0146621605286204](https://doi.org/10.1177/0146621605286204)

Besag, J., & Clifford, P. (1991). Sequential Monte Carlo p-values.
*Biometrika, 78*(2), 301-304.

## See also

[`RMdimResidualPCA`](https://pgmj.github.io/easyRasch2/reference/RMdimResidualPCA.md),
[`RMdimResidualPCACutoff`](https://pgmj.github.io/easyRasch2/reference/RMdimResidualPCACutoff.md)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
# Build 2-dimensional polytomous data: 4 items per subscale, 5 categories
n     <- 400
theta1 <- rnorm(n)
theta2 <- 0.6 * theta1 + sqrt(1 - 0.6^2) * rnorm(n)
make_pcm <- function(theta, n_items, taus) {
  sapply(seq_len(n_items), function(j) {
    # ... toy simulation here
    sample(0:4, n, replace = TRUE)
  })
}
dat <- cbind(make_pcm(theta1, 4, NULL), make_pcm(theta2, 4, NULL))
colnames(dat) <- paste0("I", 1:8)

RMdimMartinLof(dat,
            partition = list(c("I1","I2","I3","I4"),
                             c("I5","I6","I7","I8")),
            iterations = 200, parallel = FALSE, seed = 1)

# Sequential stopping: stop as soon as h = 50 simulated statistics exceed
# the observed one (cuts compute time under H0).
RMdimMartinLof(dat,
            partition = c(1,1,1,1,2,2,2,2),
            iterations = 1000, stopping = "sequential", h = 50,
            seed = 1)
} # }
```
