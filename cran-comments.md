# easyRasch2 1.5.0
## Submission

This release adds three functions: `RMitemCatFit()` and `RMitemCatFitPlot()`
test and show whether each item's response categories work as the partial
credit model predicts, `RMdifAnchor()` gives anchored Wald tests of DIF, and
`RMdifImpact()` estimates the practical size of DIF. It also fixes bootstrap
p-values that were too small when simulated values were missing, and makes
the simulation cutoff functions stop with an error instead of returning
cutoffs from a few successful simulations. Details are in NEWS.md.

Some changes alter documented output. The flag column is now named `Flagged`
in all functions, `RMdifLR()` reports only Andersen's LR test, and
`RMdifTree()` has new defaults. These are marked in bold in NEWS.md, and the
minor version is bumped for this reason.

There are no CRAN reverse dependencies.

I am aware that this package has had frequent updates (eight releases in the
past six months). It is seven weeks since 1.4.0, and the next release will
wait longer.

## Test environments

* Local: macOS v26.7.3, R v4.6.1
* GitHub Actions (R CMD check): macos-latest, windows-latest, ubuntu-latest
* check_win_devel()

## R CMD check results

0 errors | 0 warnings | 0 notes
