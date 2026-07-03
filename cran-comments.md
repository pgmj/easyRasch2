# easyRasch2 1.0.0

## Submission

This is an update (0.8.0 -> 1.0.0). The headline changes are a unified
CML/WLE estimation engine (moving `eRm` from Imports to Suggests),
new dimensionality and person/item-parameter functions, bootstrap
p-values with multiplicity correction for the simulation-based
diagnostics, and consistent sample-size reporting across all table and
figure captions. See NEWS.md for the full list.

Some functions from 0.8.0 were renamed in this release; deprecated
aliases are retained where feasible (see NEWS.md). There are no CRAN
reverse dependencies affected.

## Test environments

* Local: macOS [VERSION], R [X.Y.Z] — `R CMD check --as-cran`
* win-builder: R-release and R-devel
* macOS builder (`devtools::check_mac_release()`)
* R-hub v2: linux, windows, macos

## R CMD check results

0 errors | 0 warnings | 0 notes

[Replace with actual results; if the DESCRIPTION spell-check NOTE
appears: CFA, Rasch, infit, unidimensionality, Yen's, Q3 are correctly
spelled technical terms, proper nouns, or citation fragments.]

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
