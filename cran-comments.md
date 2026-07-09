# easyRasch2 1.1.0

## Submission

This is an update (1.0.0 -> 1.1.0). It is submitted relatively soon
after 1.0.0 because it fixes a crash: response data containing a
respondent with no responses at all (an all-NA row) caused a segfault
(dichotomous data) or an opaque error (polytomous data) in several
functions. Such rows are now dropped up front with an informative
message, consistently across the package.

The other headline change: `output = "dataframe"` now returns unrounded
values across all functions (rounding was inconsistently applied and is
a presentation concern); the formatted kable output is unchanged.
Flag labels and sorting are now computed on exact rather than rounded
values. There are no changes to function signatures and no functions
were added, removed, or renamed. See NEWS.md for the full list.

There are no CRAN reverse dependencies affected.

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
