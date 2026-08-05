# easyRasch2 1.1.1

## Submission

This is a minor update (1.1.0 -> 1.1.1). It primarily fixes a simple but troublesome bug in RMtargeting().

There are no CRAN reverse dependencies affected.

## Test environments

* Local: macOS v26.5.2, R v4.6.1
* R CMD check: macos-latest, windows-latest, ubuntu-latest
* check_win_devel() and check_mac_release()

## R CMD check results

0 errors | 0 warnings | 0 notes

All good except macos-latest and check_mac_release that produced a new error, related to a function (`RMreliability()`) that did not change in this release at all, so I suspect this is not an issue with my package.

