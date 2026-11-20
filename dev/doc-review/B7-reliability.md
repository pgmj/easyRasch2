# B7 Reliability: help-page review

2026-11-20. RMreliability, RMUreliability, RMreliabilityCurve. Roxygen only.
Rd regenerated, checkRd clean (UTF-8 note only), Rd2txt scanned,
test-documentation passes, the examples run (`run.donttest = TRUE`).

## Changes

| Page | Changes |
|---|---|
| RMreliability | **Said "three reliability indices"**, but the table has four rows: Cronbach's alpha was missing from the description, the Value section and the bootstrap sentence (alpha gets a bootstrap interval too). **Contradiction removed**: a late paragraph said "Marginal reliability here is the subtractive Green/Lord coefficient ... floored at 0" and pointed to an open question in `dev/TODO-reliability-form.md`, while the rest of the page and the code (row "Marginal (curve mean)") use the ratio form. That paragraph is gone, and the history of the subtractive coefficient is in one place. Two references to `dev/` files, which users cannot see, removed. New Data paragraph: empty rows dropped, alpha on complete cases, PSI and marginal on answered items, PSI without extreme scores. "where the respondents sit" (style list) and the announcing sentence "Three things are worth being explicit about" rewritten. Green et al. (1984) now cited. Title "Reliability coefficients for a Rasch model". References alphabetical. See also adds RMpersonParameters. |
| RMUreliability | Title sentence case. "Bayesian-style" and "silently" rephrased. Says it takes no `seed` and follows the session's random stream, as the reproducibility page states. |
| RMreliabilityCurve | `dev/milanzi_check.R` path removed. Green et al. (1984) cited (was listed only). `parallel` and `n_cores` wording aligned. One semicolon. See also adds RMpersonParameters. |

## Findings (not changed)

- **F30** RMreliabilityCurve returns its summary quantities (`sigma`,
  `marginal_ratio`, `marginal_green`, `latent_mean`, `benchmark_percent`, ...)
  only as attributes of the data.frame output (attributes rule). The kable
  shows them, but a list output would make them accessible.
- **F31** RMreliability lists Adams (2005) without citing it. Cite (for
  what?) or drop.
- **F32** RMUreliability has no example.
