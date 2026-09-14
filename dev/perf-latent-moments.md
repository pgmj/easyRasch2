# Deferred: speed up `.estimate_prior_moments()`

**Not applied.** Patch in `dev/perf-latent-moments.patch`, against
`R/utils-theta.R` as released in 1.3.1. Apply with:

```
git apply dev/perf-latent-moments.patch
```

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
