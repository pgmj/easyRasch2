# Code review: R/item_cat_fit.R + R/item_cat_fit_plot.R

Run 2026-11-19, `/code-review high`, on the uncommitted category-conditional fit code (`RMitemCatFit()`, `RMitemCatFitPlot()`), its tests, and the doc changes. Reproduction runs were in the reviewer's session and are not kept.

Status key as in utils-theta.md.

## Findings

| # | Where | Finding | Status | Reaches users? |
|---|---|---|---|---|
| CF1 | item_cat_fit.R:214 `.catfit_prepare()` | There is no extreme-only category check: it calls only `.check_complete_case_categories()`, not `.extreme_only_categories()`, unlike cfa_cutoff, dif_impact, item_parameters and utils-dif. On pcmdat2, with category 2 of I1 used only by 40 perfect scorers, the I1 thresholds are (-2.34, 13.65). `RMitemCatFit()` runs every simulation and then stops with a generic "None of the simulated datasets could be refitted" message that names no item. `type = "adjacent"`/`"residuals"` returns output built on the diverged threshold with no warning | fixed | yes, silent nonsense in plots |
| CF2 | item_cat_fit_plot.R:227, 231 | Duplicate `item_labels` give a cryptic error, "factor level [2] is duplicated". `RMitemCatFit()` accepts them | fixed | rare |
| CF3 | item_cat_fit_plot.R:163 | The `score_groups` check does not test for finite whole numbers. `Inf` gives a base-R seq() error. `2.5` gives group codes 1, 2 and 2.5, so 3 groups | fixed | rare |
| CF4 | zzz_reproducibility.R:81 | The new functions are missing from the easyRasch2-reproducibility help page. They set Mersenne-Twister per iteration, so the session RNGkind changes (L'Ecuyer-CMRG becomes Mersenne-Twister) | fixed | docs |
| CF5 | item_cat_fit.R:17 | The docs say p-values are out of `iterations + 1`. The real denominator is `iterations_used + 1`, after refit failures are dropped | fixed | docs |
| CF6 | README.md:144 | The new bullet is truncated: "four plots of where an item departs" | fixed | docs |
| CF7 | test-item_cat_fit.R:147 | `tapply(r, g, range)` puts en-dash labels into names, which gives 4 "unable to translate" warnings under a C locale | fixed | check machines with a non-UTF-8 locale |

## Status

All fixed 2026-11-19 (uncommitted). 139 expectations pass, including the parallel test, under a C locale with no warnings.

- CF1: `.catfit_prepare()` now calls `.extreme_only_categories()` on the complete cases before the CML fit and stops with an `easyRasch2_category_error` that names the item and category, the same pattern as `RMdifImpact()`. Both functions are covered, since they share the helper. Test: category 4 of I1 used only by perfect scorers, for both functions.
- CF2: `.catfit_item_labels()` refuses duplicate labels, so both functions do. Test added.
- CF3: `score_groups` must be finite and a whole number. Tests for 2.5 and Inf.
- CF4: both functions added to the list in zzz_reproducibility.R. Rd regenerated.
- CF5: `iterations` doc now says the denominator is `iterations_used + 1`.
- CF6: README bullet completed: "four plots that show where an item departs from the model."
- CF7: the test builds the score ranges with an unnamed vapply instead of tapply.

None of the fixes changes a result for valid data, so the rasch_book chapter needs no re-render.
