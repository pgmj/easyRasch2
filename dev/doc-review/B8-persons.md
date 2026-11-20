# B8 Persons: help-page review

2026-11-20. RMpersonFit, RMpersonChange, RMretestSD. Roxygen only. Rd
regenerated, checkRd clean (UTF-8 note only), Rd2txt scanned,
test-documentation passes, all three examples run (`run.donttest = TRUE`).

## Changes

| Page | Changes |
|---|---|
| RMpersonFit | Rewritten without the many "--" asides and semicolons. Title sentence case. **Two references had the same author and year** ("Sinharay (2016)" for both the resampling paper and the Psychometrika standardization paper), so the text could not say which was meant. Now Sinharay (2016a) and (2016b), cited where each applies. The p-value is now given as `(1 + count) / (B + 1)` (two-sided: twice the smaller tail), as the code computes; the page said "the proportion of the replicates at least as extreme". "marginal (uncorrected)" now "uncorrected". "resampling" now "simulation" where it means simulating from the model. de la Torre & Deng (2008) is still listed but not cited (F33). |
| RMpersonChange | Jacobson & Truax (1991) was listed but not cited, now cited for the reliable change index. "where the respondents sit" and "genuinely differ" (style list) rephrased. "percent" as "%". `parallel` and `n_cores` wording aligned. One overlong line rewrapped. One semicolon. The page is otherwise long but coherent, and was left as it is. |
| RMretestSD | Cited Caronni et al. (2026) but had **no references section**, added. `seed` now links the reproducibility page. `n_cores` said only that `mc.cores` is "checked first", now says what happens when it is unset. |

## Findings (not changed)

- **F33** RMpersonFit lists de la Torre & Deng (2008) without citing it.
  It is about correcting the ability estimate and its reference
  distribution, which is what the lz simulations do, so it could be cited
  next to Snijders and Sinharay. Cite or drop?
- **F34** RMpersonChange records the analysis settings (`null`, `retest_sd`,
  `anchor`, `alpha`, `direction`, `critical`, `sim_iter`) only as
  attributes of the data.frame (attributes rule).
- **F35** The RMpersonChange and RMretestSD examples call the internal
  `easyRasch2:::sim_partial_score()`. CRAN accepts `:::` to one's own
  package in examples, but users cannot copy the code as it is. An
  exported simulation helper, or `phq9`-based examples, would avoid it.
