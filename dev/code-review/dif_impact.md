# Code review: R/dif_impact.R + R/utils-dif-impact.R

Run 2026-11-19, `/code-review high`, on both files as new uncommitted code (after the approved review fixes). Reproduction scripts t1.R–t4.R were in the session scratchpad and are not kept.

Status key as in utils-theta.md. Fixes 2026-11-19, same day (see Status below).

## Findings

| # | Where | Finding | Status | Reaches users? |
|---|---|---|---|---|
| DI1 | dif_impact.R:517 `.dimp_estimate()` | No per-group category check on split items. If an item does not vary in one group, the shift diverges: DIF = -41.8, SE = 187575, impact SE ~0, CI [0.318, 0.318], and no warning. If focal Q2 is all 0, the function stops with the generic "small or sparse groups" error | confirmed (run) | yes, silent nonsense |
| DI2 | dif_impact.R:666 | Target rates leave out latent mean/SD uncertainty (no `Vms`, no `Jms` chain) when `target$mean`/`sd` take their defaults. Same estimate as the sample, but SE 0.0023 against 0.0079 (3.4x narrower). Contradicts the Details text | confirmed (run) | yes, CIs too narrow |
| DI3 | dif_impact.R:630 | With a theta-scale cut, the target "keep" rates use the sample-pooled keep rule. The projection uses target-refitted pooled thresholds (`pj$thr_keep`), so the two disagree | confirmed (reading) | yes, inconsistent output |
| DI4 | utils-dif-impact.R:523 | The fixed theta grid [-10, 10] (`.dimp_th_grid`, `.dimp_wle_by_score` L264, `.dimp_mml_grid`) gives NA `theta_c_sum`, so rates come back NA with no message. WLEs are clipped at +/-10 | confirmed (run, extreme thresholds) | rare, silent NA |
| DI5 | utils-dif-impact.R:361 `.dimp_mml_fit()` | No diagonal/PD check on `solve(H)`, unlike `.dimp_fit()`. Can give a "NaNs produced" warning followed by a generic stop, or NA `SE_Delta`/`Impact_in_SE` in the projection with no warning | plausible (reading) | rare |
| DI6 | dif_impact.R:321 | The category check runs on complete cases but its message reads as if it were about the full data. It reports an "empty category" that the user can see in their data | plausible (reading) | yes, confusing error |

## Status

Ordinary data give the same output to the last digit as before the fixes (two pcmdat2 runs, theta and sum cut, with a target: every table, max abs difference 0), so the Study 3c results and the earlier checks stand. test-dif_impact.R: 170 expectations (was 155), all pass.

| # | Status | What changed |
|---|---|---|
| DI1 | fixed | New `.dimp_check_shifts()`, called after the group-size check. For each split item and group it compares the group's item total with the lowest and highest totals the persons' raw scores allow, max(0, r - (R - m_j)) and min(m_j, r). At either end the shift has no finite estimate (exponential family: the sufficient statistic is on the boundary of its support). The function now stops with an `easyRasch2_category_error` that names the item, the group and the end. An item that does not vary in a group is the common case. Two tests. |
| DI2 | fixed | Target means and SDs left at their defaults now follow `ms` in the delta method: `Jms` is chained in for the item parameters, and the direct `Vms` term is added. With both at their defaults, the target rates and SEs now equal the sample ones (0.0217, against 0.0070 before). User-given target moments are fixed, so nothing changes for them. One test, and one sentence in Details. The projected impact has the same omission, but the effect is negligible there: on pcmdat2 the latent moments add SE 0.001 to 0.002 against 0.066 to 0.070 from the shifts, because keep and split see the same population. Left as is. |
| DI3 | documented (user chose option 1) | Code unchanged. Details now say that the target rates are those of the rule fitted to these data (as when a scoring table from one study is used in another population), that the projected impact refits the pooled model, and that only a theta cut is affected. Studies 3b and 3c validated the target rates with this fixed rule. |
| DI4 | fixed | New `.dimp_span()`: -10 to 10, widened in steps of 5 when a threshold lies beyond +-5, so the range reaches at least 5 logits past the most extreme threshold. The WLE bisection, the inverse true score (span +-5), the MML grid (step 0.05) and the rates grid (step 0.02, stored in the rate parts and reused by the perturbed calls) all use it. The `.dimp_range()` root search widens to mu +- 10 SD. Thresholds within +-5 keep the old grids exactly. One test, with thresholds at +-12: the true-status point is found, and the WLEs are not clipped. |
| DI5 | fixed | `.dimp_mml_fit()` returns `vcov1 = NULL` when `solve(H)` fails, is not finite, or has a non-positive diagonal. `.dimp_estimate()` now checks keep as well as split, and stops naming the cause. The projection warns when its latent-mean SE is NA. Not tested, because no small data set triggers it. |
| DI6 | fixed | `RMdifImpact()` calls the shared `.check_complete_case_categories()` (new argument `who`, whose default keeps the "This test" wording for the other callers). A category lost with the incomplete rows is now reported as "RMdifImpact() uses complete cases only (156 of 300 respondents), and in them: ...". The extreme-only message says when it refers to the complete cases. One test. |

Seen while testing, not fixed (shared helper, pre-existing): if the highest category of an item occurs only in incomplete rows, the complete cases lose it silently. `.check_response_categories()` tabulates up to the observed maximum, so the item then has fewer categories with no message. The same applies to the other callers of `.check_complete_case_categories()`: with I1's top category only in incomplete rows of pcmdat2, `RMitemRestscore()` and `RMitemInfit()` run on I1 as a dichotomous item, with only the "Model was refitted with only complete cases" message. Fixed 2026-11-19 (user approved): `.check_complete_case_categories()` compares each item's maximum in the complete cases with its maximum in the full data, and stops with an `easyRasch2_category_error` naming the categories that occur only in respondents with missing responses. This affects `RMitemInfit()`, `RMitemInfitPlot()`, `RMitemRestscore()`, `RMitemRestscoreBoot()` and `RMdifImpact()`. There is a NEWS bullet and a test in test-category-check.R. The test files of all affected functions pass (some tests skipped locally).

