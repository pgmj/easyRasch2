# easyRasch2 1.3.0

## Submission

This is a minor update (1.2.0 -> 1.3.0). The minor version is bumped because
the value of one reported statistic changes and one figure changes its default
panel.

`RMreliability()`'s marginal reliability is now the latent-density-weighted
mean of the conditional reliability curve, sigma^2 / (sigma^2 + SEM(theta)^2),
rather than Green's subtractive coefficient. The subtractive form falls below
zero whenever the average error variance exceeds the trait variance, which
happens with few items or a narrow sample, and was floored at zero. Values
therefore move upward, more so on short scales. The row is renamed
"Marginal (curve mean)" so the change is visible in the output itself, and the
previous value is still available as the `marginal_green` attribute of
`RMreliabilityCurve()`.

`RMtargeting()` now draws response-category bands in its bottom panel. Estimates
are unchanged and `panel = "thresholds"` restores the previous panel.

The release adds three exported functions: `RMreliabilityCurve()` for
conditional measurement precision across the scale, and `RMpersonChange()` with
`RMretestSD()` for assessing change in individual respondents between two
occasions.

There are no CRAN reverse dependencies. The `easyRasch2jmv` module for jamovi
depends on this package but is not distributed through CRAN.

## Test environments

* Local: macOS v26.6.1, R v4.6.1
* R CMD check: macos-latest, windows-latest, ubuntu-latest
* check_win_devel()

## R CMD check results

0 errors | 0 warnings | 0 notes
