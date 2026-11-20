# B10 Cross-cutting: help-page review

2026-11-20.

- Titles: three were still in Title Case (RMitemCatFit, RMitemCatFitPlot,
  phq9), now sentence case. The package title comes from DESCRIPTION, where
  CRAN asks for Title Case, so it stays.
- "marginal p-value" remains only in internal `@noRd` helpers. All user
  pages say "uncorrected".
- `n_cores` wording is the same on every page now. The behaviour differs
  (F5, FINDINGS.md).
- `_pkgdown.yml` covers every export except the deprecated alias
  `RMitemInfitCutoffPlot`, which is documented on the RMitemInfitPlot page.
  Nothing to change. Titles come from the Rd files.
- README: every export listed. The RMitemRestscoreBoot bullet now says what
  the function is (flag rates of the asymptotic test, a stability check,
  not a calibrated test).
- Vignette (.Rmd.orig): "threshold reversal" now "threshold disorder"
  (line 362), as on the RMtargeting page. Precompile still pending on the
  host.
- `devtools::check_man()`: no issues. `tools::checkRd()` on all Rd files:
  only the usual non-ASCII notes (the package declares UTF-8).
- All findings collected in FINDINGS.md.
