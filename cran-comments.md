# easyRasch2 1.2.0

## Submission

This is a minor update (1.1.1 -> 1.2.0). The minor version is bumped because
the default decision rule changes in three exported functions.

`RMitemInfit()`, `RMlocdepQ3()` and `RMlocdepGamma()` now flag items and item
pairs on a multiplicity-corrected bootstrap p-value rather than against the
simulated interval, and the three cutoff functions move to a common
`hdci_width = 0.95` and `iterations = 400`. Flagging against an interval tests
every item or pair at once, so its width sets a family-wise error rate
implicitly; the corrected p-value targets the stated level directly. The change
follows a simulation study, Johansson (2026),
<https://doi.org/10.31234/osf.io/7pqz4_v2>. Existing scripts keep the old
behaviour by passing `p_value = FALSE`.

The release also fixes a sampling bug in `RMdimMartinLof()` that affected
dichotomous data since 1.0.0, and makes `parallel = TRUE` and
`parallel = FALSE` return identical results for the same `seed`.

There are no CRAN reverse dependencies. The `easyRasch2jmv` module for jamovi
depends on this package but is not distributed through CRAN; it has been
updated alongside this release.

## Test environments

* Local: macOS v26.6.1, R v4.6.1
* R CMD check: macos-latest, windows-latest, ubuntu-latest
* check_win_devel()

## R CMD check results

0 errors | 0 warnings | 0 notes
