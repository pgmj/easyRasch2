# easyRasch2 1.3.1

## Submission

Apologies for submitting one day after 1.3.0. The reliability functionality
added in that release contains a bug that returns incorrect values, and it
seems better to correct it at once than to leave it in place.

`RMreliability()` and `RMreliabilityCurve()` integrate the conditional
reliability over an estimated normal latent density. Only its SD was estimated,
with the mean held at zero, so the SD absorbed any mistargeting and marginal
reliability rose as a sample became less well targeted, when it should fall. It
was overstated by up to .18 in the case checked. Well-targeted samples are
essentially unaffected. Every existing test used well-targeted data, which is
why none of them caught it, and tests for the off-target behaviour are added.

No interface changes. There are no CRAN reverse dependencies.

## Test environments

* Local: macOS v26.6.1, R v4.6.1
* R CMD check: macos-latest, windows-latest, ubuntu-latest
* check_win_devel()

## R CMD check results

0 errors | 0 warnings | 0 notes
