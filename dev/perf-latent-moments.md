# Deferred: speed up `.estimate_prior_moments()`

**Applied 2026-09-30** in the development version after 1.3.1, with one
addition the patch lacked: a log-space fallback for rows whose products all
underflow to 0. Without it such respondents became `-Inf` and were dropped by
the `is.finite()` filter, which makes a narrow, far-off candidate prior look
better. It does not happen on ordinary scales (log-likelihood ranges of a few
hundred), but on a 150-item scale 259 of 300 respondents underflowed at a
prior of N(-8, 0.05^2). With the fallback the objective matches the log-space
one to 3e-9 relative across 10 test datasets, the estimated moments match to
1e-6, and `RMreliability()` output is identical. Measured on a 300 x 10 PCM
dataset: `.latent_moments()` 0.050 s to 0.020 s per call, a 100-iteration
`RMreliabilityCurve(boot = TRUE)` 5.5 s to 2.5 s. Regression test in
`test-reliability_curve.R`.

The original patch is kept below for reference.

## Why

The 1.3.1 fix (estimating the latent mean alongside the SD) made
`RMreliability(boot = TRUE)` about twice as slow. Two causes:

1. The SD-only search was 1-D golden section over roughly 12 objective
   evaluations. The bounded 2-D search needs roughly 35, and L-BFGS-B adds
   finite-difference gradients on top.
2. Every evaluation re-exponentiated the whole person-by-node likelihood
   matrix, and `.latent_moments()` runs two quadrature passes, so that matrix
   work happens twice per call.

## What the patch does

Exponentiate once per call instead of once per evaluation. With
`E = exp(loglik - offset)` precomputed, each evaluation is the matrix-vector
product `E %*% w`, `w = dnorm(grid, mu, sigma)`, rather than a `sweep()`, a
full `exp()` and a row-wise `apply()`. The offset is a fixed upper bound on
each row's maximum, so log-sum-exp stays stable without recomputing row maxima
per candidate. The first quadrature pass, which only has to locate the mean,
drops to a coarser grid.

## Measured

300 respondents, 10 items, 1000 bootstrap iterations, sequential:

| | time |
|---|---|
| 1.3.0 | ~42 s |
| 1.3.1 as released | ~79 s |
| with this patch | ~40 s |

Results unchanged: `dev/latent_mean_check.R` reproduces to three decimals, and
`test-reliability.R` and `test-reliability_curve.R` pass.

## Before applying

Re-verify rather than trusting the numbers above. The component-level timings
I first took to attribute the cost did not reconcile with the end-to-end
figures, so only the end-to-end A/B should be relied on.
