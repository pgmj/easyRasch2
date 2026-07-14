# easyRasch2 1.1.0

## Submission

This is an update (1.0.0 -> 1.1.0). It is submitted relatively soon
after 1.0.0 because it fixes three bugs, two of them serious:

* **Crash on all-NA respondents**: response data containing a
  respondent with no responses at all caused a segfault (dichotomous
  data) or an opaque error (polytomous data) in several functions.
  Such rows are now dropped up front with an informative message,
  consistently across the package.
* **Correctness of the Martin-Löf Monte-Carlo test**:
  `RMdimMartinLof()` passed item parameters to
  `psychotools::elementary_symmetric_functions()` with the wrong sign
  convention, so the simulated null distribution did not follow the
  fitted model and the p-value depended on the arrangement of the data
  columns; `RMdimMartinLofResiduals()`'s expected counts were affected
  by the same error. The corrected sampler reproduces the exact
  conditional pattern distribution (verified against brute-force
  enumeration in the tests). **Martin-Löf p-values and residuals
  therefore change relative to 1.0.0**; this is documented prominently
  in NEWS.md. Results are now also invariant to column order and
  reproducible under a seed regardless of column arrangement.
* `RMreliability()`'s RMU estimate was not reproducible under a seed.

Beyond the bug fixes, the release adds resampling-based p-values with
multiple-comparison correction to the simulation-based fit functions,
and `output = "dataframe"` now returns unrounded values across all
functions (rounding was inconsistently applied and is a presentation
concern; the formatted kable output is unchanged, and flag labels and
sorting are now computed on exact values).

All API changes are backwards compatible: no functions were added,
removed, or renamed; a number of functions gained optional arguments
or additional `output` values, with defaults preserving prior
behaviour (except where the bug fixes above correct results). See
NEWS.md for the full list.

There are no CRAN reverse dependencies affected.

## Test environments

* Local: macOS v26.5.1, R v4.6.1
* R CMD check: macos-latest, windows-latest, ubuntu-latest
* check_win_devel() and check_win_release()

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
