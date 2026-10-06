# easyRasch2 1.4.0

## Submission

This release adds a parametric bootstrap test for the item-restscore
statistic (`RMitemRestscoreCutoff()`, `RMitemRestscorePlot()`), fixes two
simulation nulls that were too narrow (`RMdifGammaCutoff()`, and
`RMlocdepQ3Cutoff()` with missing data), and corrects a multiplicity
adjustment that was labelled Benjamini-Hochberg but was Bonferroni. Details are
in NEWS.md. It also adds a simulation-based alternative to the asymptotic 
item-restscore test, which is miscalibrated under a fitting Rasch model.

There are two small interface changes. `cutoff` is now the second argument of
`RMitemRestscore()`, and `RMdifGammaCutoff()` has new defaults. Both are listed
under "Breaking changes" in NEWS.md.

There are no CRAN reverse dependencies.

I am aware of the number of recent updates (7 in the past 6 months) and
apologise for the frequency. This release corrects results that the current
CRAN version gets wrong: a DIF null distribution that flags items too often,
a local-dependence null that is too narrow with missing data, and p-values
labelled Benjamini-Hochberg that were Bonferroni-adjusted. I would rather not
leave these in place until a later release.

## Test environments

* Local: macOS v26.6.2, R v4.6.1
* R CMD check: macos-latest, windows-latest, ubuntu-latest
* check_win_devel()

## R CMD check results

0 errors | 0 warnings | 0 notes
