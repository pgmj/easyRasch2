# easyRasch2 1.0.0

## Resubmission

This is a resubmission addressing the check-time NOTE (tests ran 11 min):
the simulation-heavy tests are now skipped on CRAN via
`testthat::skip_on_cran()` and run only when `NOT_CRAN` is set — locally
and on GitHub Actions CI (three platforms, every commit). CRAN still runs
a fast smoke suite of 500+ assertions covering input validation, error
paths, and output structure for every exported function. Test time on
local machine dropped from ~134 s to ~12 s; the full `--as-cran` check now
completes well under 10 minutes.

## Submission

This is an update (0.8.0 -> 1.0.0). The headline changes are a unified
CML/WLE estimation engine (moving `eRm` from Imports to Suggests),
new person fit and person/item-parameter functions, bootstrap
p-values with multiplicity correction for the simulation-based
diagnostics, and consistent sample-size reporting across all table and
figure captions. See NEWS.md for the full list.

One function from 0.8.0 was renamed in this release; the deprecated
alias is retained (see NEWS.md). There are no CRAN reverse dependencies
affected.

## Test environments

* Local: macOS v26.5.1, R v4.5.3 — `R CMD check --as-cran`
* R CMD check: macos-latest, windows-latest, ubuntu-latest

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for the reviewer

* **Vignette** is precompiled (chunks pre-executed locally via
  `vignettes/precompile.R`, source retained in `easyRasch2.Rmd.orig`)
  because it demonstrates Monte-Carlo simulations; rebuilding the
  shipped vignette on CRAN only formats static output and takes
  seconds.
* **Parallel processing** (optional, via `mirai` in Suggests) never
  starts more than one worker unless the user explicitly sets
  `n_cores` or `options(mc.cores)`. Examples, tests, and the vignette
  all run sequentially, respecting CRAN's two-core limit.

## Downstream dependencies

There are no reverse dependencies on CRAN.
