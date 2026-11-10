# Code review: R/utils-theta.R

Started 2026-11-10. Source: `/code-review high R/utils-theta.R` (user-run, unanchored), each finding then re-checked by me in R in the container (single datasets, small n). Reproduction scripts are not kept in the repo. The snippets below are enough to rerun each check.

Status key: **confirmed** (reproduced, wrong output or error), **plausible** (real in code, small or rare effect), **no issue** (not a problem in practice).

No fixes applied. Fixes wait for the user's go-ahead per finding.

## Summary

| # | Line | Finding | Status | Reaches users? |
|---|------|---------|--------|----------------|
| 1 | 113 | WLE clamp picks the wrong boundary | confirmed | yes, rare |
| 2 | 443 | Empty interior category drops a threshold | confirmed | yes, error or silent |
| 3 | 213 | EAP gives all-NA rows the prior mean | confirmed | yes, RMpersonChange EAP |
| 4 | 291 | Coarse first-pass grid distorts the first fit | plausible | only if the second pass fails |
| 5 | 370 | Second pass ignores the first-pass start | plausible | no, speed only |
| 6 | 402 | `sem_average` is tail-dominated at large sigma | plausible | RMreliabilityCurve caption |
| 7 | 136 | `log(exp())` in `.logp_tables` | no issue | no |
| 8 | 77 | `info_at()` duplicates `.test_information()` | plausible | no, maintenance |
| 9 | 180 | `.estimate_prior_sd` duplicates the MML objective | plausible | no, maintenance |
| 10 | 507 | Zero-row input fails in `.estimate_thetas` | confirmed | no known caller path |
| 2b | 443 | Empty lowest category: psychotools rescales the item | confirmed | yes, error or silent |
| 11 | 443 | Constant item dropped, thresholds misaligned with columns | confirmed | yes, error or silent (Q3 cutoff) |
| 12 | 442 | All-NA rows crash R in `raschmodel()` | confirmed | no public path; RMpersonChange/RMretestSD give a misleading error |
| 13 | 52 | `.theta_wle` doc says Newton-Raphson, code uses `uniroot()` | confirmed | no, doc only |

## 1. WLE clamp picks the wrong boundary (confirmed)

When the score function has the same sign at both ends of `theta_range`, line 113 picks `lo` or `hi` by comparing the raw score with half the maximum score. The sign of the score function decides which side the root is on, and the raw-score rule can disagree with it.

- `.theta_wle(0, list(7.5), c(-6, 6))` returns -6. The WLE root is +6.40.
- `.theta_wle(c(1, 0), list(9, 9), c(-6, 6))` returns -6. The root is near +9.

This happens when a respondent answered only items far from their location, or when items sit outside the range. It is rare with complete data and centred thresholds, but partial missingness and the `c(-6, 6)` bootstrap pools make it possible.

Proposed fix: `theta <- if (g_lo > 0) hi else lo` (with a guard for non-finite `g_lo`/`g_hi`).

## 2. Empty interior category drops a threshold (confirmed)

`psychotools::pcmodel()` keeps null categories by dropping their threshold (it warns "There are items with null categories"). For an item scored only 0 and 2, `.fit_cml_thresholds()` returns one threshold. Nothing downstream compares threshold counts with observed categories.

Reproduced with n = 300, 4 items 0-2, item 1 using only 0 and 2:

- `.latent_moments()` fails with "subscript out of bounds" in `.grid_loglik` (line 163).
- **`RMreliability()` fails on the public path with that same bare error message.**
- `.estimate_thetas(method = "WLE")` silently puts 8 respondents at +6, because their raw score exceeds the model maximum. None reach +6 when item 1 is removed.
- By reading only (not run): `.rasch_std_residuals()` would form residuals against E in [0, 1] for an observed 2, so Q3 would be silently wrong.

`RMpersonParameters()` uses `.rasch_fit_cml()` and warns about sparse categories. `RMpersonChange()` already guards the threshold count with `.pc_check_categories()`. So the problem is known locally but not handled in the shared helper. Bootstrap and simulation resamples are the most likely place to hit it.

Options to discuss: a check in `.fit_cml_thresholds()` that stops with a clear message (callers inside bootstrap loops already treat errors as failed iterations), or recoding empty interior categories. Recoding changes the model and is a design decision, not a bug fix.

## 3. EAP gives all-NA rows the prior mean (confirmed)

An all-NA row has a log-likelihood of zero at every node, so its EAP is the prior mean and its SEM is the prior SD. The `!is.finite(eap)` check never fires. This contradicts the doc line "Persons with no answered items return NA", and WLE does return NA.

- `.estimate_thetas(rbind(c(1,0,1), c(NA,NA,NA), c(0,1,1)), list(-1,0,1), "EAP")` gives row 2 a theta of about 0 and an SEM equal to the prior SD.
- **`RMpersonChange(method = "EAP")` with one respondent all-NA at t2:** theta_t2 = 0 (the prior mean), change = 0.65, class "none detected". With WLE the same row is NA throughout.
- `RMpersonParameters`, `RMreliability` and `RMreliabilityCurve` call `.drop_empty_respondents()` first, so they are not affected.

Side observation outside this file: `sum_t2` shows 0 for that all-NA row (it uses `rowSums(na.rm = TRUE)`). That belongs to the person_change.R review.

Proposed fix: in `.theta_eap()` (or `.estimate_thetas()`), set rows with `rowSums(!is.na(data)) == 0` to NA.

## 4. Coarse first-pass grid (plausible)

The first-pass grid has 41 nodes over [-10, 10] (spacing 0.5), but the SD search goes down to 0.05, and `.estimate_prior_moments()` does not normalise the prior weights. With a small true SD the first pass snaps the SD to 0.05 and the mean to a node.

18 datasets (5 items, n = 500, true SD from 0.1 to 1, mean 0 to 3) compared the first pass, the two-pass result and a single fine grid (801 nodes). The first pass was off in 8 of 18. The two-pass result matched the fine grid to within 0.01 in all 18. The distortion only reaches output when the second pass fails and `.latent_moments()` returns `first`.

Proposed fix (cheap): normalise `w` so that `sum(w) * dgrid == 1`, as `.marginal_summaries()` does, or put a larger SD floor on the first pass.

## 5. Second pass starts at c(0, 0) (plausible, speed only)

Both passes start L-BFGS-B at mean 0, SD 1. An off-target sample (true mean 4, SD 0.4) reached the same answer as the fine grid (3.871/0.377 against 3.871/0.377). So this has no effect on correctness. It costs extra evaluations in bootstrap loops. A low-priority fix: start the second pass at `c(first$mean, log(first$sd))`.

## 6. `sem_average` at large sigma (plausible)

`1/I(theta)` grows exponentially in the tails, so the density-weighted mean of the error variance is dominated by the tails when sigma is large. With 10 dichotomous items from -2 to 2:

| sigma | sem_average | ratio | green |
|-------|-------------|-------|-------|
| 1 | 0.784 | 0.622 | 0.385 |
| 2 | 1.085 | 0.814 | 0.705 |
| 3 | 3.047 | 0.826 | -0.032 |
| 5 | 149.7 | 0.747 | -896 |

At sigma = 5, the truncation point matters (6 sigma gives 149, 8 sigma gives 162, 10 sigma gives 162). The node count does not (161 nodes give 149.7, 321 give 149.4). So the review's claim that the result depends on `n_nodes` is wrong, but its point about the truncation span holds. `ratio` (the reported reliability) is bounded and stable. `sem_average` is shown in the RMreliabilityCurve output (reliability_curve.R:841) and is unbounded for well-spread or off-target samples. This is a property of the quantity rather than a coding error. Whether to document it, cap it, or drop it from the caption is a decision for the user.

## 7. `log(exp())` in `.logp_tables` (no issue)

On the widest grid used (-10 to 10), a 10-category item with thresholds from -5 to 5 has a minimum log-probability of -90. Underflow (below about -745) cannot occur at realistic settings. A log-scale `.pcm_cat_probs` would be cleaner but is not needed.

## 8, 9. Duplicated computations (plausible, maintenance)

- `info_at()` and `score()` in `.theta_wle()` repeat the PCM moment loop in `.test_information()`.
- `.estimate_prior_sd()` repeats the MML objective of `.estimate_prior_moments()` with the mean fixed. The two differ in how they treat non-finite rows, but `.estimate_prior_sd()` already uses log-sum-exp, so its rows only go non-finite through finding 7's underflow, which does not occur. No behavioural difference was found.

Refactor only if the file is touched anyway. Changing `.estimate_prior_sd()` would move EAP prior SDs slightly, so it needs a NEWS line if done.

## 10. Zero-row input (confirmed, low)

`.estimate_thetas(matrix(numeric(0), 0, 3), list(-1, 0, 1))` fails with "no 'dimnames' attribute for array". Top-level callers drop empty respondents and would fail earlier on empty data, so no known user path reaches this. Proposed fix: `FUN.VALUE = c(theta = 0, sem = 0)`.

## Main-session checks (2026-11-10)

Run in the main review session, independent of the `/code-review` pass above. Container, single datasets (eRm `raschdat1`, `pcmdat2`, small simulations).

### Reference comparison (no issue)

- **WLE** (`.theta_wle` via `.fit_cml_thresholds`) against the unpatched `iarm::person_estimates()` on the same psychotools fit, compared after centring: max difference 2e-7 (`raschdat1`, dichotomous) and 1e-5 (`pcmdat2`, polytomous, including score 0).
- **EAP** (`.estimate_thetas(method = "EAP", prior_sd = 1.3)`) against `mirt::fscores(method = "EAP", quadpts = 121, theta_lim = c(-6, 6))` with all mirt parameters fixed at our thresholds and N(0, 1.3²): max difference in EAP and posterior SD below 5e-15 on both datasets.
- `.latent_moments()` on an off-target sample (true mean 4, SD 0.7, 30 items, n = 3000) returned 3.99 / 0.69. `optim(L-BFGS-B)` projects the out-of-box start `c(0, 0)` into the box, which supports finding 5 being speed only.

### 2b. Empty lowest category (confirmed, extends finding 2)

When an item's lowest category is unobserved (pcmdat2 item 1 recoded to use only 1 and 2), `pcmodel()` warns "Minimum score is not zero ... These items are scaled to zero" and returns one threshold for an item still coded 1/2 in the data. `.estimate_thetas(WLE)` then puts 37 respondents at the +6 boundary (0 with item 1 removed), and `RMreliability()` fails with the bare "subscript out of bounds" from `.grid_loglik`. An empty top category is harmless: the item gets one threshold fewer and the data agree with it. So the check proposed for finding 2 should compare each item's threshold count with `max(observed) - min(observed)` and also require `min(observed) == 0`, or the fix has to cover both cases.

### 11. Constant item drops out of the threshold list (confirmed)

For an item with only one observed category, psychotools warns "unidentified items" and **omits** it, so `.fit_cml_thresholds()` returns a list one shorter than `ncol(data)` and every later item is paired with the wrong thresholds. Tested on pcmdat2 (item 2 all 0) and raschdat1 items 1-10 (item 3 all 1):

| Function | Result |
|---|---|
| `RMreliability()` | clear stop: "only one response category" (the only guarded path) |
| `RMlocdepQ3()` | "subscript out of bounds" |
| `RMpersonParameters()`, `RMpersonFit()` | "'names' attribute [10] must be the same length as the vector [9]" |
| `RMitemInfitCutoff()`, `RMitemRestscoreCutoff()` | "All simulation iterations failed. Check your data." |
| `RMlocdepQ3Cutoff()` | **runs**: simulates 9 items (36 pair cutoffs) while `item_names` lists 10, so the cutoff is calibrated for a different scale than the observed Q3 |

`.sparsity_warning()` flags the item but only says estimates "may be unstable". Proposed fix: stop in `.fit_cml_thresholds()` when `length(thr_list) != ncol(data)`, naming the constant item(s), the same message `RMreliability()` already gives. Could share one check with findings 2/2b.

### 12. All-NA respondents crash R in `raschmodel()` (confirmed, latent)

One all-NA row passed to `psychotools::raschmodel()` (psychotools 0.7.7) segfaults in `elementary_symmetric_functions` and kills the R session. `pcmodel()` instead errors with "subscript out of bounds". Every exported caller that reaches `.fit_cml_thresholds()` removes such rows first (tested 11 functions with one all-NA row, dich and poly: no crash). So nothing reaches the crash today, but it would take down a jamovi session if a future caller skipped the filter. Proposed fix: drop all-NA rows inside `.fit_cml_thresholds()`; they carry no information for CML.

Side observation for the person_change.R review: `RMpersonChange()` and `RMretestSD()` with a respondent all-NA at both time points stop with "CML estimation failed (subscript out of bounds). This often indicates sparse data; try estimator = "MML"", which points the user at the wrong cause.

### 13. Stale doc (confirmed, doc only)

The `.theta_wle()` roxygen (line 52) says "Newton-Raphson with adaptive step damping". The code solves with `stats::uniroot()` on `theta_range`.

## Fixes applied (2026-11-10, main session, uncommitted)

All in `R/utils-theta.R` unless noted. Pre-fix copy kept outside the repo only for the session.

| # | Fix |
|---|-----|
| 1 | Clamp side from the sign of the score function at `lo`; raw-score rule kept only as fallback for non-finite values |
| 2, 2b, 11 | New `.check_cml_categories()` called in `.fit_cml_thresholds()`: stops naming items with one category, an empty lowest or an empty interior category. Plus a post-fit check that the threshold list matches `ncol(data)` and each item's observed maximum |
| 12 | `.fit_cml_thresholds()` drops all-NA rows before fitting |
| 3 | `.estimate_thetas(EAP)` sets all-NA rows to NA and leaves them out of the prior-SD estimate; `.theta_eap()` doc corrected |
| 4 | `.estimate_prior_moments()` normalises prior weights over the grid |
| 5 | Second pass of `.latent_moments()` starts at the first-pass answer (new `start` argument) |
| 8 | `.theta_wle()` uses `.test_information()` instead of its own `info_at()` |
| 10 | Named `FUN.VALUE` in the WLE `vapply()` |
| 13 | `.theta_wle()` doc says `uniroot()` |

Not done: 6 (`sem_average` at large sigma, a user decision), 9 (`.estimate_prior_sd` refactor, would move EAP prior SDs, so left for a deliberate change with a NEWS line), and the person_change.R misleading error (belongs to that unit).

**Effect on ordinary data** (raschdat1, pcmdat2, pcmdat2 with 60 NAs, phq9 items 1-9, an off-target and a narrow simulation): thresholds, WLE and EAP identical; latent moments change by at most 1e-4 (off-target SD 0.45847 to 0.45836); largest change in any `RMreliability()` value 5.8e-5.

**Failure cases after the fix:** all four clamp cases go to the correct side; the three reshaped-item cases and a combined one stop with the item named, including through `RMlocdepQ3()`, `RMlocdepQ3Cutoff()` and `RMreliability()`; empty top categories and extreme-only categories still fit; first pass on a narrow sample (true SD 0.1) now 1.29/0.14 instead of snapping to the floor, two-pass result unchanged and equal to an 801-node grid.

**Tests:** new `tests/testthat/test-utils-theta.R` (20 expectations, covers 1, 2, 2b, 3, 10, 11, 12, top-empty alignment, off-target latent moments). NEWS: three bullets under Bug fixes.

## Not yet done for this unit

Per the review plan, still open:

- ~~Reference comparison of WLE/EAP against psychotools/eRm/mirt (dich + poly).~~ Done, see main-session checks.
- `.wle_theta_pool()` and `.rasch_std_residuals()` read only through finding 2.
- Test gaps: no direct tests for `.center_thresholds` or `.wle_theta_pool`. Suggested tests for findings 1, 2, 3 and 10 once fixed.

## /code-review high on the fix diff (2026-11-10)

Verified in the main session. Statuses as above. None of these are fixed yet.

| # | Where | Finding | Status |
|---|-------|---------|--------|
| R1 | item_restscore.R:335 and others | The category guard runs only inside `.fit_cml_thresholds()`. Direct `pcmodel()` + `threshpar()` call sites are unguarded: item_restscore.R:335, conditional_infit.R:368, bootstrap_restscore.R:195, conditional_infit_mi.R:230, item_parameters.R:284 (`.rasch_fit_cml`, used by RMpersonParameters/RMpersonFit/RMpersonChange), plus targeting_plot.R:668 and item_hierarchy.R:162. | confirmed |
| R2 | test-utils-theta.R:6, 77 | Uses eRm (Suggests) without `skip_if_not_installed("eRm")`. Other test files that use eRm have the skip. Fails a noSuggests check. | confirmed |
| R3 | NEWS.md:20 | "and others" overstates the guard's coverage (see R1). "Results for other data do not change" leaves out the changes up to 1e-4 in latent moments and 5.8e-5 in RMreliability values (fixes 4 and 5). | confirmed |
| R4 | utils-theta.R:178 | `.estimate_prior_sd()` still uses unnormalised `dgrid` weights, so the two prior estimators now follow different rules. This is finding 9 from another angle. The snapping risk at the SD floor is unconfirmed. | plausible (= 9) |
| R5 | utils-theta.R:105 | The clamp rule assumes the Warm score function is monotone. A non-monotone one with disordered thresholds could have interior roots that are missed when both ends have the same sign. This predates the fix (`uniroot` finds only one bracketed root). Not reproduced. | plausible, rare |
| R6 | utils-theta.R:498 | `.check_cml_categories()` duplicates the category scan in `.sparsity_warning()` (item_parameters.R:445): one stops, the other warns. Same point as the open item_parameters.R note. | confirmed, maintenance |
| R7 | utils-theta.R:467 | The post-fit alignment check cannot fire once `.check_cml_categories()` has passed. It is a defensive assert costing one `apply()` per bootstrap refit. Keep it with a comment, or drop it. | confirmed, low |
| R8 | utils-theta.R:456 | A third copy of the all-NA-row rule (also at :580 and in `.drop_empty_respondents()`). | confirmed, maintenance |

**R1 reproduction:** pcmdat2 with item 1's category 0 recoded to 1 (empty lowest category). `RMitemRestscore()` runs and returns a table. The only signal is psychotools' "Minimum score is not zero ... scaled to zero" warning, and the person mean behind `Relative_location` uses misaligned thresholds (the review measured mean 1.71, with 37 respondents clamped at +6). `RMpersonParameters()` also runs, with the sparsity warning plus the psychotools warning, and returns WLEs from the reshaped thresholds. These are pre-existing bugs, not regressions from the fix, but the NEWS bullet implies they are covered.

**Suggested resolution for R1, R6 and R8 (needs a decision):** one shared helper (drop all-NA rows quietly, check categories, fit, `threshpar()`, return centred thresholds and optionally the fit and vcov) that every CML call site uses. `.sparsity_warning()` would keep only the warn-level shapes (sparse but non-empty categories, extreme-only categories). This crosses into item_parameters.R and five other files, so it may belong to the item_parameters.R unit.

### Main-session check of R1-R8 (2026-11-10)

- **R1, confirmed, list incomplete.** Besides the seven sites listed, unguarded `pcmodel()`/`raschmodel()` calls are at bootstrap_restscore.R:432, item_restscore.R:411, infitcutoff_plot.R:268, infit_cutoff.R:406, item_restscore_cutoff.R:417 and martin_lof.R:542-629 (dichotomous path and subscale fits). Several refit simulated data, where the observed data already passed the guard through `.wle_theta_pool()`. Not every site necessarily misaligns: sites that only use psychotools' own fitted values may handle reshaped items internally. So each site needs checking before it gets the guard. MML has the same problem (item_parameters.md IP2).
- **R2, confirmed, FIXED.** eRm is in Suggests only. Added `skip_if_not_installed("eRm")` to `pcm_data()` and the latent-moments test. 20/20 pass.
- **R3, confirmed, FIXED.** NEWS now lists the covered functions (each traced from `.fit_cml_thresholds()`, `.wle_theta_pool()`, `.rasch_std_residuals()`, `.marginal_rxx()` or `.wle_psi()`; `RMdimMartinLof()` polytomous and `RMitemCatProb()` were run to confirm), says other CML functions are not yet covered, and has a separate bullet for the at-most-1e-4 shift in the latent moments.
- **R4** = finding 9, left as is.
- **R5, the non-monotonicity is real, but the missed-root case was not found.** On a 0.01 grid over [-10, 10], with heavily disordered random thresholds (1-4 items, 2-5 thresholds, spreads up to ±5), the score function was non-monotone in 2211 of 3000 patterns and had more than one root in 167. Over [-6, 6] with realistic shapes: 0 multi-root cases in 1500 for ordered thresholds (4-9 items), 0 for one reversed pair, and 9 and 1 in 1500 with one or two answered items (disordered, ordered). **In no case did both ends share a sign with roots inside**, which is the only case the clamp affects. When the ends differ in sign, `uniroot()` returns one of several roots, which predates the fix. Plausible, rare. Not fixed.
- **R6, R8** maintenance, go with the R1 decision.
- **R7, plausible.** Two cases where psychotools might still reshape after the pre-check passes, both extreme-only (lowest category used only by all-zero scorers; a dichotomous item whose 1s come only from perfect scorers), fitted aligned and the assert did not fire. Keeping it costs one `apply()` against a full `pcmodel()` fit per refit, so I'd keep it with a comment.
