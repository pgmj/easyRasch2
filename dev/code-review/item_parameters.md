# Code review: R/item_parameters.R

Started 2026-11-10, main session. `/code-review` pass on this file not yet run. Checks in R in the container on eRm `raschdat1` (dichotomous, 30 items), `pcmdat2` (4 items, 0-2) and a mixed set (6 dichotomous + 4 polytomous items, n = 100). Single datasets.

Status key as in `utils-theta.md`. No fixes applied. The user commits the utils-theta.R work first.

Callers: `.rasch_fit_cml()` is used by RMitemParameters, RMpersonParameters, RMpersonFit, RMpersonChange, RMretestSD and RMscoreSE. `.rasch_fit_mml()` is used by the same functions except RMscoreSE. `.sparsity_warning()` is also called from cicc_plot.R, local_dependence.R and residual_pca.R.

## Summary

| # | Line | Finding | Status | Reaches users? |
|---|------|---------|--------|----------------|
| IP1 | 284 | `.rasch_fit_cml()` has no category guard (= R1 in utils-theta.md) | confirmed | yes, silent |
| IP2 | 339 | MML path reshapes items with empty categories the same way | confirmed | yes, silent |
| IP3 | 465 | Sparsity warning advises MML, which does not help with empty categories | confirmed | yes, wrong advice |
| IP4 | 373, 143 | MML SEs are not propagated through the centring | confirmed | yes, SEs 2-19% too large |
| IP5 | 297 | CML SEs and centring | no issue | |
| IP6 | 261, 264 | Stale comments: "(eRm)" header, "delta-method SEs" for CML | confirmed | no, doc only |
| IP7 | 445 | `.sparsity_warning()` duplicates `.check_cml_categories()` (= R6) | confirmed | no, maintenance |

## IP1. No category guard in `.rasch_fit_cml()` (confirmed)

The same psychotools behaviour as findings 2, 2b and 11 of utils-theta.md. `pcmodel()` drops the threshold of an empty interior category, rescales an item with an empty lowest category, and omits a constant item. `.rasch_fit_cml()` only calls `.sparsity_warning()`, so:

- **Empty interior category** (pcmdat2 item 1 recoded to 0/2): `RMitemParameters()` reports I1 with one row, "threshold 1", location 0.01, SE 0.099, for an item whose data run 0-2. `RMpersonParameters()` (CML) puts 56 respondents at the +10 boundary.
- **Empty lowest category** (item 1 recoded to 1/2): `RMpersonParameters()` puts 37 respondents at +10 (also reproduced by the diff review as R1, via `RMitemRestscore()`).
- **Constant item**: `setNames()` fails with "'names' attribute [4] must be the same length as the vector [3]".

## IP2. MML reshapes the same items (confirmed)

`mirt::mirt()` recodes item categories to be consecutive from 0, so `.rasch_fit_mml()` returns one threshold for the 0/2 item and for the 1/2 item, while the data keep the original codes. `RMpersonParameters(estimator = "MML")` puts 56 (interior) and 37 (lowest) respondents at +10, the same counts as CML. A constant item stops with mirt's clear "only one response category" message.

So the problem is how the data are represented downstream, not the estimator. A guard has to run for both estimators.

## IP3. The warning recommends MML for empty categories (confirmed)

`.sparsity_warning()` flags constant items and empty categories and says "CML estimates may be unstable; consider estimator = "MML"". By IP2, MML gives the same wrong person estimates for empty categories. The advice fits sparse but non-empty categories only.

## IP4. MML SEs ignore the centring (confirmed)

With `center = TRUE` (default), MML thresholds are shifted to a grand mean of zero, but `.mml_threshold_se()` propagates only each item's own block of the mirt covariance, on mirt's scale (latent mean fixed at 0). The SE of a centred threshold needs the full covariance, `C J V J' C'` with `C = I - 11'/n`. Computed that way:

| Data | Reported SE | SE of centred thresholds | Reported / centred |
|---|---|---|---|
| raschdat1 (dich) | 0.240-0.339 | 0.216-0.317 | 1.07-1.11 |
| pcmdat2 (poly) | 0.157-0.339 | 0.139-0.291 | 1.04-1.19 |
| mixed | 0.213-0.486 | 0.202-0.455 | 1.02-1.07 |

The centred dichotomous MML SEs (0.216-0.317) are close to the CML SEs (0.216-0.316), as expected for the same contrast. The current Wald CIs for MML are therefore too wide. The per-item mapping itself is correct: `sqrt(diag(J V J'))` equals the reported SEs in all three datasets, so the row alignment of `d` parameters is right, including mixed formats. With `center = FALSE` the current SEs are the correct ones for mirt's scale.

Proposed fix: build `J` over all items, form `J V J'`, and apply `C` when `center = TRUE`. That means `.rasch_fit_mml()` returns the full threshold covariance (or the centring happens inside it), since centring now happens in `RMitemParameters()`.

## IP5. CML SEs and centring (no issue)

`psychotools::threshpar()` already returns thresholds with grand mean 0 (both datasets), and the SEs of centred thresholds equal the reported ones exactly. Dichotomous SEs equal `eRm::RM()` and `psychotools::itempar()` sum-zero SEs (0.2161-0.3159). PCM SEs differ from `eRm::thresholds()` (0.142-0.356 against 0.148-0.337), as the code comment says. eRm's threshold SEs depend on its own constraint, so this is not evidence of an error here. The roxygen claim "CML estimates are already centred" holds.

## IP6. Stale comments (confirmed, doc only)

- Line 261: section header "CML fit + threshold extraction (eRm)". The code uses psychotools.
- Line 264: "(+ delta-method SEs)". The CML SEs come straight from `threshpar(vcov = TRUE)`, with no delta step.

## IP7. Duplicate category scan (= R6, maintenance)

`.sparsity_warning()` and the new `.check_cml_categories()` scan the same shapes, one to warn and one to stop.

## Proposed resolution (needs a decision)

IP1, IP2, IP3 and IP7 together with R1, R6 and R8 from utils-theta.md point to one change. A shared check (the current `.check_cml_categories()`) runs before every CML and MML fit, including the seven direct `pcmodel()` call sites R1 lists. `.sparsity_warning()` then keeps only the warn-level shapes: sparse but non-empty categories and extreme-only categories. Its MML advice applies only to those. Optionally, one shared CML fit helper (drop all-NA rows, check, fit, `threshpar()`, return centred thresholds plus optional vcov) replaces the direct calls, as R1 suggests.

IP4 is independent and local to this file. IP6 is a doc edit.

## Not yet done for this unit

- `/code-review high R/item_parameters.R` in a fresh session.
- Test gaps: no direct test of `.rasch_fit_cml()`; no test of MML SEs. Suggested: SE of centred MML thresholds against `C J V J' C'`; the IP1/IP2 shapes once guarded.

## Shared category check: decisions and frequency simulation (2026-11-10)

User decision: option (b), one shared CML fit helper used by every site. Simulation first.

**Script:** `dev/code-review/category_frequency/category_frequency.R`, results in `category_frequency.rds` (+ dated copy), log `run.log`. Container run: R = 100 observed datasets per cell, B = 50 resamples or splits each, 4 cores, 4.3 min. Design: n 100/200/400, 6 or 10 items, 4 or 5 categories, theta ~ N(0, 1) (targeted) or N(-1, 1) (one logit off target). Ordered thresholds: item location U(-1, 1), spread U(0.8, 1.4) over [-1.5, 1.5], jitter SD 0.2. One threshold generator only, so treat rates as indicative.

Columns, in %: **obs_fail** = observed data the new check would stop. **obs_any_empty** = any empty category, including an empty top. **par_any_empty** = parametric `*Cutoff()` simulations (CML thresholds + resampled WLE pool + `sim_partial_score()`) the existing validation discards ("not all categories represented"). **par_fail** = those the new check would stop. **boot_fail** = row-bootstrap resamples the new check would stop. **split50/70** = random 50/50 or 70/30 DIF split with any empty category in a group (eRm::LRtest then drops the item), and the mean number of items dropped. par/boot/split are among datasets that passed the check.

| n | items | cats | mu | obs_fail | obs_any_empty | par_any_empty | par_fail | boot_fail | split50 % (items) | split70 % (items) |
|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 6 | 4 | 0 | 1 | 2 | 7 | 3 | 6 | 37 (0.45) | 56 (0.79) |
| 200 | 6 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 5 (0.05) | 16 (0.17) |
| 100 | 10 | 5 | 0 | 7 | 11 | 29 | 14 | 19 | 84 (1.67) | 92 (2.54) |
| 200 | 10 | 5 | 0 | 0 | 0 | 4 | 2 | 4 | 24 (0.29) | 45 (0.63) |
| 400 | 10 | 5 | 0 | 0 | 0 | 0 | 0 | 0 | 2 (0.02) | 10 (0.10) |
| 100 | 6 | 4 | -1 | 0 | 42 | 56 | 0 | 0 | 90 (1.80) | 94 (2.20) |
| 200 | 6 | 4 | -1 | 0 | 13 | 25 | 0 | 0 | 69 (1.04) | 77 (1.30) |
| 200 | 10 | 5 | -1 | 0 | 47 | 67 | 0 | 0 | 94 (2.61) | 96 (3.13) |
| 400 | 10 | 5 | -1 | 0 | 18 | 31 | 0 | 0 | 67 (1.15) | 77 (1.53) |

(Selected rows. Full 24-cell table in the rds `summary`.)

**What it shows**

1. **The new check rarely stops observed data:** at most 7% at n = 100 with 10 five-category items, and 0-1% from n = 200. Off target, the empty categories are top categories, which the check allows, so it stopped 0-3%.
2. **The existing simulation validation is the bigger issue.** It discards any simulated dataset with an empty category, including an empty top category. Off target that is 56-93% of iterations at n = 100, 25-67% at n = 200 and 3-31% at n = 400. The returned null is conditioned on every category being observed, which off-target data often do not satisfy themselves. Under the new rule (empty top allowed) discards would be 0-14%. This touches all six `*Cutoff()` functions (infit, restscore, Q3, gamma LD, residual PCA, CFA). Whether the conditioning moves the cutoffs, and whether allowing an empty top is valid for each statistic, needs its own null check. Not decided. Record for the Cutoff units.
3. **The shared helper never fires inside the parametric refits,** because the existing validation runs first and is stricter. So the point-2 worry (new failures biasing the null) does not arise there. Row-bootstrap refits (RMreliability, RMreliabilityCurve, RMitemRestscoreBoot) would lose up to 19% of resamples at n = 100 and 0-4% from n = 200.
4. **DIF splits empty a category very often.** A 50/50 split leaves an empty category in some group in 37-84% of targeted datasets at n = 100, and in 22-100% off target at every n tested. On average 0.2-4.6 items are then dropped by `eRm::LRtest()`. `RMdifLR()` reports such items as `Flagged = no` (see below), so this is common, not an edge case.

## RMdifLR and empty categories (checked 2026-11-10, for the dif_lr unit)

- A category empty in one group (any position, or a constant item): `eRm::LRtest()` drops the item from the full and group models ("Full and subgroup models are estimated without these items!"). `RMdifLR()` shows NA group columns, an `All` location from the separate `eRm::PCM(data)` fit that kept the item, and `Flagged = no`. Only eRm's raw warning signals it. The other items' `All` values also come from that fit, so they are on a slightly different centring than the group columns.
- An empty interior category in the whole sample: eRm warns "Estimation may not be feasible", fits anyway, and returns SEs of 13-210 logits. No stop.
- Needed: the shared check on the full data, and an explicit "not tested (empty category in group X)" status. Collapsing categories consistently across groups is the common alternative and a design decision.

## Shared helper, step 1 applied (2026-11-10, uncommitted)

- `R/utils-theta.R`: new `.cml_fit(data, model = c("pcm", "auto"), hessian = TRUE)`, the single psychotools entry point. It drops all-NA rows, runs the check, then fits (`"auto"` = `raschmodel()` for all-dichotomous data). `.check_cml_categories()` renamed `.check_response_categories()`, since it now also guards MML. Its error has class `easyRasch2_category_error`. `.fit_cml_thresholds()` now goes through `.cml_fit()`, and its post-fit assert has a comment (R7).
- `R/item_parameters.R`: `.rasch_fit_cml()` uses `.cml_fit()`, and its `tryCatch()` passes the category error through without the "try MML" advice. Gotcha: a separate `easyRasch2_category_error =` handler does not work, because `tryCatch()` nests handlers, so the generic `error` handler re-caught the re-raised error. The fix is one handler that tests the class. `.rasch_fit_mml()` runs the check before `mirt()` (IP1, IP2). IP4: `.mml_threshold_vcov()` (replaces `.mml_threshold_se()`) keeps the full threshold covariance, both fits return `vcov_thr`, and `RMitemParameters()` computes SEs as `sqrt(diag(C V C'))` when `center = TRUE`. IP6 comments fixed. Roxygen: SE paragraph and a "Response categories" paragraph; `man/RMitemParameters.Rd` regenerated with roxygen2 8.0.0 (matches `Config/roxygen2/version`).
- **Baseline comparison** (raschdat1 items 1-15, pcmdat2, pcmdat2 with 7% NA, mixed 6 dich + 4 poly, phq9 items 1-9; RMitemParameters CML/MML x long/wide/no-centre, RMpersonParameters WLE/EAP/MML, RMpersonFit, RMscoreSE, RMpersonChange, RMretestSD): everything identical (< 1e-10) except `RMitemParameters(estimator = "MML")` with centring, where only `se`/`ci_*` moved. MML SE before/after/CML: raschdat1 0.241-0.341 / 0.214-0.310 / 0.215-0.306; pcmdat2 0.157-0.339 / 0.139-0.291 / 0.142-0.356; mixed 0.213-0.486 / 0.202-0.455 / 0.202-0.454.
- **Failure cases:** empty middle, empty lowest and constant item now stop with the item named through RMitemParameters, RMpersonParameters, RMpersonFit, RMpersonChange, RMretestSD (CML and MML) and RMscoreSE. An empty top category still fits.
- **Tests:** 3 new blocks in `test-item_person_parameters.R`. Full suite: 0 fail, 1662 pass, 7 skip, 1 warning (the same mirt EM warning in test-reliability.R:100 as before). NEWS: coverage bullet extended, new MML SE bullet.
- **Still to do:** step 2 sites (restscore, infit, infit MI, restscore boot, infit plot, hierarchy, targeting, Martin-Löf dichotomous + subscales), `RMdifLR()` full-data check, then shrinking `.sparsity_warning()` (its MML advice is still wrong for empty categories, but those now stop before it matters).

## /code-review high on the step-1 diff (2026-11-10), verified

Run as a forked `/code-review high` on the uncommitted diff (plus the committed utils-theta work). Ten findings, labelled D1-D10. Checked in the main session with scratch scripts (not kept).

| # | Finding | Status | Severity |
|---|---|---|---|
| D1 | NaN in the threshold vcov spreads to every SE through `J V J'` and `C V C'` | confirmed mechanically, low in practice | low |
| D2 | `RMscoreSE(method = "EAP")` calls `mirt()` with no category check, NEWS claims coverage | confirmed | user-facing |
| D3 | Remaining direct psychotools sites (targeting, hierarchy, infit, restscore, restscore boot, infit plot, infit MI, Martin-Löf dich) | confirmed, = planned step 2 | user-facing |
| D4 | Bootstrap resamples that lose a sparse category now error and are dropped | confirmed | moderate |
| D5 | `.estimate_prior_sd()` prior weights not normalised | confirmed, edge case | low |
| D6 | mirai export lists lack `.cml_fit`, `.check_response_categories`, `.test_information` | not a user bug | none |
| D7 | All-NA input gives an arbitrary `prior_sd` attribute | confirmed, unreachable from exported functions | low |
| D8 | `.restscore_gamma()` has its own category assert | reuse point confirmed, coverage claim wrong | maintenance |
| D9 | Redundant all-NA row drop before the check, row copy in `.cml_fit()` | confirmed | cosmetic |
| D10 | Dense O(n^3) centring products | not worth changing | none |

**D1.** `J %*% V %*% t(J)` with one NaN in `V` gives 36 of 36 NaN entries (6 x 6 test), and `diag(C V C')` gives 6 of 6. The old per-item blocks kept a NaN inside one item. For MML with `center = TRUE` the spread is correct in substance: a centred threshold's variance includes `mean(V)`, so one undefined covariance leaves every centred SE undefined. The review's suggested `diag(V) - 2 * rowMeans(V) + mean(V)` spreads NaN the same way. mirt usually fails the whole information matrix rather than one block, so partial NaN is rare. Possible fixes: build `J V J'` block by block (only matters for `center = FALSE`), and skip the `C` product for CML, whose thresholds are already sum-zero (IP5), so it changes nothing there except NaN spread.

**D2.** Simulated PCM data, 400 x 5 items, 0-3. Clean data: `raw_score` 0-15. Same data with category 1 of V1 merged into 2: EAP runs without error, `raw_score` 0-14. mirt renumbered V1 to 0-2, so every row of the score table is labelled with a score the user's data do not use. WLE on the same data stops with the named-item message. The NEWS bullet lists `RMscoreSE()` "with both CML and MML", which is wrong for EAP. Fix: call `.check_response_categories()` in `.scoreSE_eap()` (or once in `RMscoreSE()` before the method branch), or narrow the NEWS bullet.

**D3.** The planned step 2. The review's list matches it: `targeting_plot.R` `.estimate_thresholds_cml`, `item_hierarchy.R`, `conditional_infit.R`, `item_restscore.R`, `bootstrap_restscore.R`, `infitcutoff_plot.R`, `conditional_infit_mi.R`, dichotomous `martin_lof.R`. NEWS already says "Other functions that fit the model by CML are not yet covered".

**D4.** Comes from the committed utils-theta fix (the stop in `.fit_cml_thresholds()`), not this diff. Before that fix such resamples ran with misaligned thresholds, which was worse. One item with 2 of 400 responses in category 1: 12.9% of row-bootstrap resamples lose the category. `RMreliability()` drops them, and the caption shows the reduced count ("N bootstrap resamples"), so it is not fully silent, but the reason is not given and the doc says only "fail to converge". The CI is then conditioned on resamples containing every category. Same in `RMreliabilityCurve()` and `RMitemRestscoreBoot()`. Matches point 3 of the frequency simulation (up to 19% at n = 100). Options: (a) message the count and reason, (b) collapse the empty category in the resample for that iteration, (c) leave as is and document. Needs a decision.

**D5.** `RMpersonParameters()` passes its `theta_range` (default `c(-10, 10)`) to the EAP grid, so node spacing is 0.167 with 121 nodes. At sigma = 0.05 the Riemann sum of `dnorm * dgrid` is 1.34 (1.01 on `c(-6, 6)`). Only matters when the marginal likelihood favours a SD near the floor, i.e. near-zero latent variance. At realistic SDs the sum is 1 to many digits, so normalising like `.estimate_prior_moments()` would not move normal results. One-line fix.

**D6.** Functions exported to mirai keep the package namespace as their environment, so their nested calls resolve in the namespace loaded on the daemon, not in the export list. With an installed package everything is found. The export list only matters for top-level names in the daemon expression. Under `load_all()` with a stale installed version the daemons run old code anyway, export list or not. No user-facing issue.

**D7.** `.estimate_thetas(all-NA, method = "EAP")` returns `prior_sd` 4.99999 (the top of the `optimize()` interval). `RMpersonParameters()` drops empty respondents first, so exported functions do not reach this. Return `NA` when no row is answered, if touched.

**D8.** The local assert compares each item's max response with its threshold count after `pcmodel()`. The review says it misses an empty lowest category, but psychotools then returns one threshold fewer while the max stays at the original code, so the assert fires. The reuse point holds: replace with `.cml_fit()` in step 2, which gives the named-item message.

**D9.** `.check_response_categories()` drops NA per column, so the `rowSums(!is.na())` subset before it in `.rasch_fit_mml()` does nothing. `.cml_fit()` copies the matrix even when `keep` is all `TRUE`. Cosmetic, `if (!all(keep))` is enough.

**D10.** 400 thresholds is about 6e7 flops per product, well under a second, and the call is once per fit, not per bootstrap iteration. Not worth changing.

**Proposed for this unit:** D2 (check in RMscoreSE EAP) and D5 (normalise) now, D1 CML-skip optional. D3 and D8 in step 2. D4 needs a decision. D6, D7, D9, D10 no action.

### Second opinion on D1-D10 (main session, 2026-11-10)

- **D4 is not new.** Rerun with the pre-commit `utils-theta.R` (b013f7e) swapped into the namespace: on a 400-row item with 2 responses in category 1, 14.3% of 300 row-bootstrap resamples lose the category, and **the old code failed exactly the same 14.3%** (in `.latent_moments()`, "subscript out of bounds"). No resample that lost the category survived under the old code. The fix changed the error text, not which resamples are dropped. `RMitemRestscoreBoot()` does not use the committed check (it calls `pcmodel()` directly, a step-2 site), so the claim that it is affected by the commit is also wrong. The underlying issue (a CI conditioned on every category being present, reason not reported) is real and predates the review.
- **D2 is right but too narrow.** Four mirt fits run without the check, not one: `RMscoreSE(method = "EAP")` (score_se.R:295), `RMlocdepQ3()` (local_dependence.R:522, 532), the RMU path of `RMreliability()` (reliability.R:380) and `RMtargeting()` (targeting_plot.R:697). Add them to the step-2 inventory with the psychotools sites.
- **D5 confirmed:** `sum(dnorm(grid, 0, 0.05)) * dgrid` = 1.340 on the 121-node [-10, 10] grid, 1.000 at SD 1.

## Point 1, step 2 and D4a applied (2026-11-10, uncommitted)

**Point 1.** D5: `.estimate_prior_sd()` normalises its prior weights (log-sum-exp of `lprior`, `dgrid` dropped). EAP estimates on the baseline datasets move by at most 1.2e-6. D1: `RMitemParameters()` applies the centring product to MML only, so a partial NaN in a CML covariance stays local. CML SEs identical. NEWS bullet narrowed, then superseded by the step-2 wording.

**Step 2.** Observed-data fits switched to `.cml_fit()`: RMitemInfit (conditional_infit.R:368), RMitemInfitMI (conditional_infit_mi.R:230), RMitemRestscoreBoot full fit (bootstrap_restscore.R:195), RMitemInfitPlot (infitcutoff_plot.R:268), RMitemHierarchy (item_hierarchy.R:162), RMitemRestscore (item_restscore.R:335), RMtargeting CML (targeting_plot.R:668), RMdimMartinLof dichotomous sampling parameters (martin_lof.R:542, `model = "auto"`). `.check_response_categories()` added before every mirt fit (`.q3_residual_matrix()` MML, RMreliability RMU, RMscoreSE EAP, RMtargeting MML) and before `eRm::PCM()/RM()` in RMdifLR. RMitemInfitMI's per-imputation handler passes the category error through (it was hidden behind "Model fitting failed for all 2 imputed datasets").

Deliberately left as direct psychotools calls, each a refit inside a simulation or resample that uses only the fit object itself (logLik or iarm), which psychotools keeps self-consistent: `compute_ml_statistic()` and the subscale WLE fits (martin_lof.R), `run_single_boot_restscore()` (bootstrap_restscore.R:432), the infit and restscore Cutoff refits (already behind the stricter "not all categories represented" validation), and the complete-case refit in RMitemRestscore (item_restscore.R:411), whose `.restscore_gamma()` refuses reshaped data inside a `tryCatch()`. Not covered: RMdifTree (psychotree, reports node rescaling itself).

Baseline (13 functions x 4 datasets: raschdat1 items 1-12, pcmdat2, pcmdat2 with 7% NA, phq9 items 1-9): all 52 comparisons identical, including the three pre-existing input-limit errors. Failure cases: empty interior category (pcmdat2 item 1) and a constant dichotomous item stop with the item named in all 14 paths tested.

**D4a.** New `.notify_dropped_resamples()` (utils-validation.R), called in RMreliability and RMreliabilityCurve after the bootstrap. It informs "N of B bootstrap resamples were dropped", how many lacked a category and for which items, and the most common other reason. Silent when nothing failed. Example (item with 2 of 400 responses in a category, 60 resamples): RMreliability 8 dropped, RMreliabilityCurve 7, both naming I1. Roxygen updated ("silently dropped" line replaced), Rd regenerated with roxygen2 8.0.0.

**Tests:** new `tests/testthat/test-category-check.R` (16 expectations). Full suite: 0 fail, 1678 pass, 7 skip, 1 warning (the usual mirt EM warning, test-reliability.R:100). NEWS: coverage bullet rewritten ("every function that fits a Rasch or partial credit model", RMdifTree excepted), D4a bullet.

**Open:** shrink `.sparsity_warning()` (it still warns "consider MML" just before the category error in four functions); RMdifLR "not tested" status for items eRm drops within a group (dif_lr unit); the Cutoff simulations' any-empty validation (Cutoff units).

## /code-review high on the step-2 diff (2026-11-10), not yet verified

Ten findings, most severe first. Only E1 was reproduced by the reviewer. Nothing changed yet.

**Wrong results without any warning**

**E1.** `run_single_boot_restscore()` (bootstrap_restscore.R:458) reads `iarm::item_restscore()` rows by position. When a resample has a single-category item, psychotools drops it, so every later item gets its neighbour's results and the last becomes NA. Reproduced with 6 items and V3 constant: V3 got V4's values, V4 got V5's, V6 became NA. An item with 2 endorsements in 400 rows and `samplesize = 300` is constant in about 22% of resamples. This contradicts the step-2 note that this refit is self-consistent. Fix: match on `rownames(res_mat)`, or treat the resample as failed when `nrow(res_mat) != n_items`.

**E2.** The complete-case refit in RMitemRestscore (item_restscore.R:411) still calls `pcmodel(cc)` directly, and `tryCatch(error = NULL)` hides failures. If the complete-case subset loses a category the full data have, the function quietly falls back to iarm's rounded Observed/Expected, which are then compared with an unrounded null. Fix: run `.check_response_categories(cc)` and stop or tell the user.

**E3.** The mirai export lists (reliability.R:728, reliability_curve.R:745) lack `.cml_fit` and `.check_response_categories`. Under `load_all()` with an older installed version, every resample fails with "could not find function" and the notifier reports them all as dropped for "another reason". Fix: add both helpers, or drop the hand-maintained lists.

**Messages that mislead or contradict each other**

**E4.** `.sparsity_warning()` (item_parameters.R:304) fires just before the category error and suggests `estimator = "MML"`, which now stops on the same data. Already listed as open above.

**E5.** Roxygen (item_parameters.R:232), NEWS ("MML SEs and CIs get smaller") and the IP4 test `all(raw$se >= res$se)` say centring always shrinks MML SEs. Not guaranteed: Var(t_i - tbar) < Var(t_i) only when 2 Cov(t_i, tbar) > Var(tbar). Fix: "usually smaller", drop the `all()` check.

**E6.** The once-per-session asymptotic warning in RMitemRestscoreBoot (bootstrap_restscore.R:182) fires before the category check. If that call stops on the check, the warning is spent and the corrected call that returns flags shows none. Fix: move it after the full-sample fit.

**E7.** `.notify_dropped_resamples()` (utils-validation.R:775) says "The interval uses the remaining 1" and is followed by "band not drawn". It also says "interval" in RMreliabilityCurve, where the output is a band.

**Design and fragility**

**E8.** The notifier (utils-validation.R:771) finds dropped items by parsing the error message text. It breaks if the wording changes and truncates names like `Q1 (rev)`. The test builds the message by hand, so it would not catch a rewording. Fix: carry item names as a field on the classed condition.

**E9.** Resamples that lose a sparse category are still dropped (reliability.R:496), so the CI or band comes only from resamples that kept it (about 13-14% dropped with one 2-of-400 category). Same issue as D4 above, which predates the review. A real fix belongs in the resampling step (collapse or stratify). Needs a decision.

**E10.** The pre-fit check exists in three forms (inside `.cml_fit()`, inline with a row filter in `.rasch_fit_mml()` at item_parameters.R:376, and bare elsewhere). Maintenance cost, not a bug. Overlaps D9.

**Proposed order:** E1-E3 first (wrong output without any warning), then E4-E7. E9 needs a decision.

### Main-session check of E1-E10 (2026-11-10)

- **E1, not reproduced as stated, fix still right.** With iarm 0.4.3, `iarm::item_restscore()` on a fit where psychotools dropped a constant item returns a 5-row `noquote` character matrix, with item names as row names. The worker's `res_mat[seq_len(6), 1L]` then errors ("subscript out of bounds"), so the resample is dropped, not shifted. Seen directly (6 dichotomous items, V3 constant) and in 40 resamples of a 2-in-400 item: no shifted rows. A data.frame result would shift silently, so the reviewer may have had another iarm version. Either way: match rows on `rownames(res_mat)` and give a dropped item NA for that resample. That keeps the other items' results, which the current code throws away with the whole resample. The step-2 note's "self-consistent" claim was wrong for this site, because the worker reads by position.
- **E2, confirmed, mechanism refined.** pcmdat2 with I2 missing for every respondent with I1 = 1: complete cases have I1 counts 76/0/80. `pcmodel(cc)` reshapes I1, `.restscore_gamma()` refuses (its max/threshold check), and the function falls back to iarm's rounded values for **all** items. I1 gets NA Expected, Difference and p. The only signals are psychotools' "null categories" and R's "NAs introduced by coercion" warnings. Proposal: run `.check_response_categories(cc)`. On failure, inform naming the item ("I1 cannot be tested on complete cases: category 1 occurs only in incomplete responses") and keep the other items. Stopping would be stricter than the full-data rule. Decision.
- **E3, confirmed, dev-only.** The export lists are complete by intent, so add `.cml_fit` and `.check_response_categories` to all four (reliability.R:728, reliability_curve.R:744, local_dependence.R:1422, residual_pca.R:803), not just the two named. Installed users are unaffected (D6).
- **E4** known (shrink `.sparsity_warning()`).
- **E5, confirmed in principle.** Var(t_i - tbar) < Var(t_i) needs 2 Cov(t_i, tbar) > Var(tbar). It held in all three test datasets (ratios 1.02-1.19), but it is not guaranteed. Reword to "usually smaller" in roxygen and NEWS, and drop the `all()` expectation.
- **E6, agree.** Move the once-per-session warning after the full-sample `.cml_fit()`.
- **E7, agree.** With fewer than 2 survivors the notifier should not say "the interval uses the remaining 1". Say "results" or "band" per caller.
- **E8, agree.** Carry the item names on the classed condition (`items` field), and have the bootstrap workers return them as an attribute of the failure string (attributes survive mirai serialization). Item names containing " (" are truncated today.
- **E9** = D4, decided as (a).
- **E10, partly.** The row filter in `.rasch_fit_mml()` is redundant, since the check ignores NA cells and all-NA rows cannot change category counts. Remove it. The bare calls before mirt/eRm are the intended form.

### E1-E10 fixes applied (2026-11-10, uncommitted)

User decision: E2 stops (the user decides how to handle the category).

- **E1:** `run_single_boot_restscore()` matches iarm rows by `rownames(res_mat)`. An item that is constant in a resample gets NA there, and `table()` with fixed levels leaves it out of that item's percentages, while the other items keep the resample. If iarm ever returns unnamed rows for a reduced item set, the resample errors instead of being read by position. Test: V3 constant, V6's diff equals iarm on the same resample.
- **E2:** new `.check_complete_case_categories()` (utils-theta.R) stops with "This test uses complete cases only (n of N respondents), and in them: ..." in RMitemInfit, RMitemRestscore, RMitemRestscoreBoot and RMitemInfitPlot. **RMitemInfit had the same problem in a worse form:** with I1 category 1 only in incomplete rows (pcmdat2), it returned infit MSQ 0 for all four items with only psychotools' "null categories" warning. The comment at the RMitemRestscore complete-case refit was updated.
- **E3:** `.cml_fit` and `.check_response_categories` added to all four mirai export lists (reliability, reliability_curve, local_dependence, residual_pca).
- **E4:** `.sparsity_warning()` keeps only the extreme-only warning. The constant/empty part (with its MML advice) only ever preceded the category error. `@param estimator` in RMitemParameters reworded.
- **E5:** "usually smaller" in roxygen and NEWS. The `all(raw$se >= res$se)` expectation was dropped.
- **E6:** the once-per-session asymptotic warning in RMitemRestscoreBoot now follows the complete-case check and the full-sample fit.
- **E7:** the notifier takes `what = "interval"/"band"` and says "Too few remain for a ..." below 2 survivors.
- **E8:** the category error carries `items`. The bootstrap workers return `structure(conditionMessage(e), items = e$items)`, and the notifier reads the attribute (names with " (" survive).
- **E10:** the redundant row filter in `.rasch_fit_mml()` was removed.
- **E9** unchanged (= D4, decided (a)).

Baselines after the fixes: step-2 set 52/52 and item/person set 60/60 identical to before. Tests: `test-category-check.R` extended (E1, E2, E7, E8, items field). The three affected test files pass 116/116. NEWS: E2 bullet, E1 bullet, MML-advice removal, "usually". Rd regenerated (RMitemParameters, RMitemRestscoreBoot) with roxygen2 8.0.0.

Full suite after the E fixes: 2 failures, both tests that expected the removed "sparse/zero-variance" warning (test-local_dependence.R:409, test-residual_pca.R:166, constant item). Rewritten to expect the category error naming the item. Those two files then pass 119/119. Rest of suite: 1687 pass, 7 skip, 1 warning (the usual mirt EM warning).

## /code-review high on the E-fix diff (2026-11-10), not yet verified

Ten findings, most severe first. The reviewer confirmed two points before writing them up: `psychotools::threshpar()` returns sum-zero thresholds whose covariance rows sum to 0, so skipping the centring for CML is safe, and an item with an empty lowest category in a resample stays aligned in the iarm refit in `run_single_boot_restscore()`. Nothing else was reproduced. Nothing changed yet.

**Wrong results or wrong advice**

**F1.** The complete-case check runs before the full-data check (conditional_infit.R:361, item_restscore.R:322, bootstrap_restscore.R:198, infitcutoff_plot.R:255). A category that is empty in the full data is then blamed on complete cases, and the message suggests imputing, which cannot fix it. Example: pcmdat2 with I1 category 1 recoded to 2 and one NA anywhere. Fix: run `.check_response_categories()` on the full data first (see F9).

**F2.** RMitemRestscoreBoot (bootstrap_restscore.R:345) now gives each item its own denominator (resamples where iarm dropped the item or returned NA are left out), but the caption still names `actual_iterations` as the base. An item that is NA in every resample gets percent NaN, `percent > cutoff` (line 376) is NA, and all-NA rows enter the kable. Example: 1 endorsement in 400 rows, `samplesize = 50`, `iterations = 5` (constant in about 88% of resamples).

**F3.** The MML centring product C V C' (item_parameters.R:162) spreads one non-finite entry in `vcov_thr` to every threshold SE and CI. Before the diff only the affected item was NaN. Fix: centre only when `all(is.finite(vcov_thr))`, else keep the uncentred per-item SEs.

**F4.** `test-category-check.R:28` calls functions that need Suggests packages without skipping: iarm (RMitemInfit, RMitemRestscore, RMitemRestscoreBoot, RMitemInfitMI, and the E1 test at line 127) and ggdist (RMreliability). `mid_empty()` only skips on eRm. Fails under a noSuggests check. Fix: `skip_if_not_installed("iarm")` / `("ggdist")`.

**Messages and documentation**

**F5.** The E7 notifier (utils-validation.R:188) prints "Too few remain for a interval.", and the regex at test-category-check.R:93 accepts the typo. Fix: "an interval" / "a band", or pass the article with `what`.

**F6.** The E6 move was not far enough. The once-per-session asymptotic warning (bootstrap_restscore.R:205) still fires before the bootstrap, so a call that then stops with "All bootstrap iterations failed" (line 298) spends it. Fix: move it after `successful` is non-empty.

**F7.** The new roxygen (bootstrap_restscore.R:80) says only an item that is constant in a resample is left out. An item whose interior category is empty in a resample is also left out (iarm returns NA Expected), and `@return` does not mention the NA item_restscore and diff rows in `output = "raw"`.

**Efficiency and maintenance**

**F8.** `.fit_cml_thresholds()` (utils-theta.R:463) calls `.cml_fit(data, model = "auto")` with the default `hessian = TRUE` but uses only point estimates. This is the busiest path (RMreliability and RMreliabilityCurve resamples, Q3 and residual-PCA simulations, polytomous Martin-Löf). Fix: `hessian = FALSE`.

**F9.** The complete-case check plus a separate `.cml_fit()` call is written out at four call sites, all in the F1 order. Fix: a `complete_cases = TRUE` argument on `.cml_fit()` that runs the full-data check first. Overlaps E10 and D9.

**F10.** `.sparsity_warning(data, is_poly)` (item_parameters.R:489) no longer uses `is_poly`, but its callers (item_parameters.R:304, cicc_plot.R:287, residual_pca.R:209, local_dependence.R:318) still pass it. Fix: drop the argument or call `.extreme_only_categories()` directly.

**Proposed order:** F1 (with F9), F3, F4, F5 first. They are small and contained, and F1 gives wrong advice.

### Main-session check and fixes of F1-F10 (2026-11-10, uncommitted)

All ten confirmed by reading the code. F8 was measured: `hessian = FALSE` gives identical thresholds (max diff 0) and is 1.4x (pcmdat2) to 5x (raschdat1) faster.

- **F1 + F9:** `.check_complete_case_categories()` now runs `.check_response_categories()` on the full data first, so a category empty in the full data gets the plain message, not "complete cases only ... impute". The four call sites are unchanged (call order no longer matters), so no new `.cml_fit()` argument was needed. Test: mid_empty plus one NA gives the plain message.
- **F2:** RMitemRestscoreBoot adds `n_tested` per item, `percent` is NA (not NaN) when 0, the kable filter skips NA, and the caption lists items tested in fewer resamples. Test: one endorsement in 400, samplesize 50, so V3 has `n_tested` < 5, V1 = 5, and no NaN. `@return` documents `n_tested` and the NA raw rows (F7).
- **F3, partly as proposed:** when `vcov_thr` has non-finite entries, the per-item SEs on mirt's scale are kept **with a warning** that they are not adjusted for centring. Silently mixing scales would mislabel the SEs.
- **F4:** skips for iarm, ggdist and ggplot2 in test-category-check.R.
- **F5:** "an interval" / "a band", with the test regex tightened.
- **F6:** the once-per-session warning moved after `actual_iterations` (at least one success).
- **F7:** roxygen says "cannot be tested in a resample (constant, or with an empty category there)".
- **F8:** `.fit_cml_thresholds()` uses `hessian = FALSE`.
- **F10:** `.sparsity_warning(data)`, with the four callers updated.

Baselines: step-2 set 52/52 and item/person set 60/60 identical (the new `n_tested` column excluded from the comparison). Affected test files 289/289. NEWS: RMitemRestscoreBoot bullet extended. Rd regenerated (RMitemRestscoreBoot).

Full suite after the F fixes: 1 failure, test-rounding.R:88, which requires every numeric column passed to `.round_display()` to have a rule or be allow-listed. The new integer `n_tested` was added to `allowed_unrounded`, next to `n`. The rounding tests then pass. Rest: 1694 pass, 7 skip, the usual mirt warning.

## /code-review high on the F-fix diff (2026-11-10), not yet verified

Nine findings, labelled G1-G9. The reviewer checked two points in R before writing them up. First, psychotools shifts an item with an empty lowest category in a resample down to start at 0, and the item-restscore result is identical to recoding it by hand, so that item is still tested and its numbers are valid (G1 is a documentation error only). Second, CML thresholds are already sum-zero, so centring SEs for MML only is correct. No tests were run. Nothing changed yet.

**Behaviour and documentation that disagree**

**G1.** RMitemRestscoreBoot (bootstrap_restscore.R:364). The roxygen (line 85), NEWS and the caption note say an item that is constant or has an empty category in a resample is left out of its percentages. Only constant items and empty interior categories give NA. An item with an empty lowest category is shifted by psychotools, tested and counted in `n_tested`, so the documented base is wrong for that case.

**G2.** RMitemRestscoreBoot (bootstrap_restscore.R:512). When an item is constant in a resample, iarm drops it and applies BH over K-1 items. Those resamples are now pooled with K-item resamples, so every item's flag rate moves slightly (example: 8 items, one constant in 30% of resamples, the other 7 get BH over 7 tests there). Before the E1 fix such resamples were dropped whole, so the families were not mixed.

**G3.** `.notify_dropped_resamples()` (utils-validation.R:199) says the interval or band "is therefore based on resamples that contain every category". A resample without an item's top category passes `.check_response_categories()` and is kept with one threshold fewer, so the message is wrong and the selection is asymmetric.

**G4.** RMitemInfitMI (conditional_infit_mi.R:255). The comment assumes an empty or single category is the same in every imputation, and one imputation's category error now aborts the function. That holds for pmm but not for polr, polyreg or norm, where an imputation can fill a rare interior category in some imputations and not others. Before, that imputation was skipped with a warning and the rest pooled.

**G5.** `.sparsity_warning()` (item_parameters.R:313) still runs before the category check in `.rasch_fit_cml()`, RMlocdepQ3, RMdimResidualPCA and RMitemICCPlot. Data about to stop with a category error can first give the extreme-only warning (empty interior category with the top category used only by perfect scorers). The new tests hide this by wrapping `expect_error()` in `suppressWarnings()`. Fix: warn after `.cml_fit()` / `.fit_cml_thresholds()`.

**Comments and cleanups**

**G6.** Stale comments at local_dependence.R:314 and residual_pca.R:206 still say the CML path warns about sparse or zero-variance categories and points to MML (removed in E4).

**G7.** The full-data category check runs twice: `.check_complete_case_categories()` calls `.check_response_categories()` (F1), then `.cml_fit()` runs it again. Sites: bootstrap_restscore.R:203, conditional_infit.R:361/370, item_restscore.R:322/337, infitcutoff_plot.R:255/270. Fix: let `.check_complete_case_categories()` check only complete cases and call `.cml_fit()` first, or the reverse.

**G8.** `.cml_fit()` (utils-theta.R:498) is called "the single entry point for psychotools CML fits", but the bootstrap worker (bootstrap_restscore.R:493), Martin-Löf `compute_ml_statistic()` and subscale fits (martin_lof.R:561-629), the Cutoff simulation workers and the complete-case refit (item_restscore.R:414) still call `pcmodel()`/`raschmodel()` directly, and each mirt path adds its own `.check_response_categories()`. A new call site without the manual check would bring back the silent threshold misalignment. Proposal: a matching `.mml_fit()` wrapper, and route the remaining direct psychotools calls through `.cml_fit()` (with a no-check option for simulated-data workers). **Deferred by the user to a separate step, after G1-G7 and G9 are fixed and committed.**

**G9.** RMreliability (reliability.R:496) and RMreliabilityCurve (reliability_curve.R:659). With fewer than 2 resamples left, the notifier says "Too few remain for an interval" and the existing warning "Fewer than 2 bootstrap iterations succeeded; CI not reported." fires too. Fix: drop one of them.

**Proposed order:** G1-G5 (what users see or what results mean), then G6, G7, G9. G2 and G4 need a decision (G2: drop resamples with a constant item again, or accept mixed BH families and document. G4: skip the failing imputation again, or keep the stop). G8 after the commit.

### Main-session check and fixes of G1-G9 (2026-11-10, uncommitted)

Checked in R with iarm 0.4.3 and psychotools 0.7.7 (pcmdat2, raschdat1):

| Case in a resample | psychotools | iarm `item_restscore()` | Effect in RMitemRestscoreBoot |
|---|---|---|---|
| empty lowest category | shifts item to 0 | valid, identical to hand recode | tested, counted |
| empty interior category | drops threshold | row kept, expected and p NA; others BH over K | item NA, others keep the resample |
| empty top category | fine | valid | tested |
| constant, polytomous | drops item | **errors** ('names' attribute) | whole resample fails |
| constant, dichotomous | drops item | K-1 rows, BH over K-1 | item NA, others BH over K-1 (G2) |

- **G1 confirmed, and broader:** roxygen, the worker comments, the caption note and NEWS now say what each case does. The caption parenthetical was dropped rather than listing the cases.
- **G2 confirmed, dichotomous only** (polytomous constant still drops the whole resample, as before E1). Needs a decision.
- **G3 confirmed, fixed:** the message now says the interval "leaves out resamples that lost an item's lowest or a middle category". The roxygen for RMreliability, RMreliabilityCurve `boot_iter` and the helper were changed to match.
- **G4 confirmed by reading:** the comment holds for pmm, not for polr/polyreg. Needs a decision.
- **G5 confirmed, fixed:** `.sparsity_warning()` moved after the fit in `.rasch_fit_cml()`, RMlocdepQ3 (CML), RMdimResidualPCA and RMitemICCPlot. New test: an empty interior category whose top category only a perfect scorer uses gives the category error and no warning in all four functions.
- **G6 fixed** with G5. The stale comments are gone, and no others remain.
- **G7 confirmed, not changed:** the check is a per-item table, which costs nothing next to the fit. Removing the duplicate would make the call order matter again, the problem F1/F9 removed. Recommend leaving as is.
- **G8 deferred** (user).
- **G9 confirmed, fixed:** with fewer than 2 left the notifier gives only the count. The caller's warning says the interval or band is not reported. Test updated.

Rd regenerated (RMitemRestscoreBoot, RMreliability, RMreliabilityCurve). test-category-check.R 43/43. Full suite: 0 failures, 7 skips, the usual mirt warning.

### G2 and G4 decisions applied (2026-11-10, uncommitted)

User's decisions: G2 (a), G4 (b).

- **G2 (a), accept and document.** In a resample where a dichotomous item is constant, iarm drops it, and the BH adjustment covers the K-1 items tested there. That is correct for the tests actually run, and dropping those resamples would condition on the sparse item being varied. Documented in the RMitemRestscoreBoot roxygen ("@details") and the NEWS bullet. No code change.
- **G4 (b), skip the imputation.** In `RMitemInfitMI()` the handler now collects category errors instead of re-raising them. After the loop:
  - every imputation lacked the category: stop with the first category error (full message with remedy; the pmm case, where the observed data lack it);
  - otherwise one warning: "k of m imputed datasets lacked a response category ... (I1 in k). They were skipped, so the pooled results come only from imputations that filled that category. Consider merging the sparse category with an adjacent one.";
  - fewer than 2 left: the existing "Only n imputation(s) succeeded" stop, with that note appended;
  - the generic "failed and were excluded" warning counts only the other failures.
  The roxygen (failure paragraph) and NEWS (new bug-fix bullet) describe it. `RMitemInfitCutoffMI()` already skipped a failing imputation with a warning, so the two now agree.
- **Tests:** the existing MI test (renamed "stops when every imputation lacks the category") still gets the category error. New test: m = 4 pmm imputations with category 1 put into imputations 1-3 by hand. It returns a data frame, with the "1 of 4 ... lacked" warning naming "I1 in 1" and no generic failure warning.

Rd regenerated (RMitemRestscoreBoot, RMitemInfitMI). With NOT_CRAN=true: test-category-check 47/47, test-infit_mi 40/40, test-bootstrap_restscore 11/11, test-rounding 6/6. The full suite was not rerun after these changes. It had passed after the G1-G9 fixes, and these changes touch only RMitemInfitMI and documentation.

Open after this round: G8 (route the remaining direct pcmodel/raschmodel and mirt calls through `.cml_fit()` and a new `.mml_fit()`), deferred by the user until after the commit. Also RMdifLR "not tested" status (dif_lr unit) and the Cutoff any-empty validation (Cutoff units).

## G8 applied (2026-11-10, after the user's commit, uncommitted)

Every psychotools and mirt fit in the package now goes through `.cml_fit()` or the new `.mml_fit()` (both in utils-theta.R). A test fails if a direct call appears anywhere else.

- **`.cml_fit(..., check = TRUE)`:** new `check` argument. `check = FALSE` skips the category check (all-NA rows are still dropped). Two kinds of refit use it, because they rely on psychotools reshaping an item:
  - the RMitemRestscoreBoot resample worker, where an item with an empty middle category becomes NA, an empty lowest category is shifted, and a constant dichotomous item is dropped (G1/G2);
  - `compute_ml_statistic()` in Martin-Löf, for the observed data and the simulated datasets. The observed data were already checked by `extract_ml_sampling_params()`.
- **`.mml_fit(data, ..., check = TRUE)`:** runs `.check_response_categories()` and then `mirt::mirt(data, model = 1, itemtype = "Rasch", verbose = FALSE, ...)`. No rows are dropped, since callers match `fscores()` and residuals to the input rows. It replaces the check-then-mirt pairs in `.rasch_fit_mml()`, `.q3_residual_matrix()` (both branches), RMreliability, RMscoreSE (EAP) and RMtargeting `.estimate_thresholds_mml()`.
- **Sites moved to `.cml_fit()`:**
  - bootstrap_restscore worker (`hessian = FALSE, check = FALSE`);
  - RMitemInfitCutoff sim worker (default hessian, as before);
  - RMitemRestscoreCutoff sim worker (`hessian = FALSE`);
  - RMitemRestscore complete-case refit with a cutoff (`hessian = FALSE`);
  - Martin-Löf `compute_ml_statistic()` (full and subscale fits);
  - `compute_subscale_wle_and_correlations()`.
  In the Cutoff workers the polytomous validation (every category present) is stricter than the check. **Corrected after H3:** the dichotomous validation (8 positive responses) lets an item answered 1 by everyone through, and the check then stops that iteration. It was discarded before G8 too, so results do not change.
- **Mirai export lists:** `.cml_fit` and `.check_response_categories` were added for bootstrap_restscore, infit_cutoff, item_restscore_cutoff and martin_lof. `.mml_fit` was added for the Q3 cutoff.
- **Not in scope:** `eRm::PCM()/RM()` in RMdifLR (already checked, dif_lr unit) and psychotree in RMdifTree.

**Gotcha found on the way (dev only):** mirai puts the exported objects in the daemon's global environment. A package function serialized to a daemon still looks names up in its namespace first, and in the daemon that is the **installed** easyRasch2 (1.4.0 in the container). New helpers that are absent from the installed version fall through to the exported copies and work. A helper whose signature changed resolves to the installed old version. A first version gave `compute_ml_statistic()` a `check` argument, and every parallel Martin-Löf simulation failed with "unused argument (check = FALSE)", because the daemon called the installed 1.4.0 copy. Fixed by keeping its signature and setting `check = FALSE` inside. With a current installed package this does not arise. It does mean parallel tests under `load_all()` test a mix of the dev code and the installed package. **Corrected after H2:** this understated it. Under `load_all()`, a parallel worker runs the installed version of *every* helper it calls that the installed build also has, so the parallel rows of the G8 comparison below tested installed code for those helpers. Those rows were rerun with each version installed (see the H section).

**Verification:** the committed version (HEAD, exported with `git archive` to the scratchpad) and the working tree, on the same data:
- the step-2 baseline: 52/52 identical;
- the item/person baseline: 60/60 identical;
- a new G8 baseline: 60/60 identical. It covers RMitemInfitCutoff and RMitemRestscoreCutoff (sequential and parallel), RMitemRestscore with a cutoff on NA data, RMitemRestscoreBoot (sequential and parallel), Martin-Löf (sequential and parallel), the MML Q3 cutoff (parallel), RMtargeting's sparse MML path, RMscoreSE EAP and RMitemParameters MML. It includes a set where resamples lose an item's lowest category and a dichotomous set where I12 was tested in only 34 of 40 resamples.

New tests in test-category-check.R: `.mml_fit()` and `.cml_fit()` stop on an empty interior category, `check = FALSE` returns a pcmodel, and a source scan finds no direct `pcmodel()`, `raschmodel()` or `mirt::mirt()` call outside utils-theta.R (skipped when the R sources are absent). Full suite: 0 failures, 1713 pass, 7 skips, 1 known mirt warning. No NEWS entry, since no results change.

## /code-review high on the G8 diff (2026-11-10), verified in the main session

| # | Where | Finding | Status |
|---|---|---|---|
| H1 | bootstrap_restscore.R:503 | The daemon worker calls `.cml_fit(..., check = FALSE)`. It resolves to the installed namespace, and HEAD's `.cml_fit` has no `check` argument. | **Confirmed, dev only** |
| H2 | martin_lof.R:566 | Same mechanism in `compute_ml_statistic()`. | **Not confirmed as stated** |
| H3 | infit_cutoff.R:406, item_restscore_cutoff.R:417, this log | "The validation is stricter than the check" is false for dichotomous data. | **Confirmed, comments and log only** |
| H4 | test-category-check.R:249 | The source-scan test misses eRm, psychotree, `do.call(psychotools::pcmodel, ...)`, unqualified calls and code after a `"#"` string. | **Confirmed, low** |
| H5 | martin_lof.R:901, bootstrap_restscore.R:582 | `.check_response_categories` is exported to daemons that never call it. | **Confirmed, low** |
| H6 | utils-theta.R:543 | `.mml_fit()` has an unused `check` argument. Hard-coded `verbose = FALSE` collides with `...`. | **Confirmed, design** |
| H7 | item_restscore.R:414 | The complete-case refit runs the category check again after `.check_complete_case_categories()`. | **Confirmed, but not in a loop** |
| H8 | martin_lof.R:626 | Repeats the model-choice expression stored as `ml_model` in `compute_ml_statistic()`. | **Confirmed, low** |

**H1.** HEAD (dab3a66, 1.4.0.9001) was installed to a scratch library with `git archive` + `R CMD INSTALL`, and the working tree was loaded with `load_all()`. A daemon reports `formals(.cml_fit)` as `data, model, hessian`. A namespace-environment closure sent with `.cml_fit = .cml_fit` exported also sees the installed version. `RMitemRestscoreBoot(d, iterations = 8, parallel = TRUE, n_cores = 2)` stops with "All bootstrap iterations failed". This only happens while the installed build is older than the loaded tree. It disappears once the working tree is installed, and CRAN users always have a single version. The container never showed it because its installed 1.4.0 has no `.cml_fit`, so the exported copy was used.

**H2.** The same setup gives identical Martin-Löf results whether HEAD is installed or not (T_obs 7.89, p .667, 20/20 draws). `run_single_ml_iteration` resolves `compute_ml_statistic` from the installed namespace too. HEAD's copy is internally consistent, so nothing fails. The daemon simply runs HEAD's code. The real consequence is wider than the finding: **under `load_all()`, a parallel worker runs the installed versions of every helper it calls that also exists in the installed build.** The exported copies only matter for helpers the installed build lacks. So the G8 parallel baselines (and any parallel `load_all()` test) compared installed code against installed code for those callees, not the working tree. The G8 gotcha paragraph above understates this.

**H3.** For dichotomous data the Cutoff validation only requires at least 8 positive responses per item. An item answered 1 by everyone passes it and then stops in `.check_response_categories()`. The iteration is still discarded (before G8, pcmodel dropped the item and `data.frame()` failed), so results do not change. Only the error string differs. The comment at both sites and the G8 bullet "the check never fires there" are wrong for dichotomous data.

**H4.** Of the misses, RMdifLR's `eRm` fit is guarded explicitly (dif_lr.R:226). RMdifTree's psychotree fits have no category check (dif_tree unit). The rest are hypothetical future code.

**H7.** Line 414 is in `RMitemRestscore()`, which no Cutoff loop calls (RMitemRestscoreCutoff uses `.restscore_gamma()` directly). It costs one extra pass per call, not per iteration.

**Proposed** (nothing applied):
- H1/H2: no code change. Reinstall the working tree before any parallel check, and correct the G8 gotcha paragraph. Optional: a dev-note line, or a test helper that skips parallel tests when the installed version differs from the loaded one.
- H3: fix the two comments to say the polytomous validation is stricter, and that a dichotomous item with no 0 responses passes validation and then stops in the check (iteration discarded either way).
- H4: widen the regex to `(psychotools::)?(pcmodel|raschmodel)\b` and `eRm::(PCM|RM)\(`, with an allow-list for dif_lr.R. Leave psychotree for the dif_tree unit.
- H5, H8: drop the dead exports and pass `ml_model` in. Both are cheap.
- H6: drop `check` from `.mml_fit()` until a caller needs it. Leave `verbose` alone.
- H7: `check = FALSE` at line 414, optional.

### H1-H8 applied in the main session (2026-11-10, uncommitted)

Checked against the code. H1 and H2 match the mechanism I saw in the G8 Martin-Löf failure, and H3 matches the dichotomous validation (`colSums < 8` only).

- **H1/H2:** no package code change. The G8 gotcha paragraph and the Cutoff bullet above are corrected. The parallel comparison was then redone properly. HEAD (1.4.0.9001) and the working tree were each installed to their own scratch library (`R CMD INSTALL -l`), and the G8 baseline was run with `R_LIBS` pointing at the matching library, so the daemons run that version's own code. Result: 60/60 identical, all parallel paths included (RMitemRestscoreBoot, RMitemInfitCutoff, RMitemRestscoreCutoff, Martin-Löf, MML Q3 cutoff). **Recipe for parallel checks in the container:** install the tree to a scratch library and set `R_LIBS` to it, since daemons inherit it.
- **H3:** the comment at both Cutoff sites now says that the polytomous validation is stricter and that a constant dichotomous item stops in the check (discarded either way).
- **H4:** the scan now strips strings before comments and matches unqualified `pcmodel`/`raschmodel`, `mirt::mirt` and `eRm::PCM`/`RM` (dif_lr.R allow-listed). Checked against the pre-G8 sources: it finds all 16 old direct calls. psychotree is left for the dif_tree unit.
- **H5:** `.check_response_categories` was removed from the Martin-Löf and RMitemRestscoreBoot daemon exports (both call `.cml_fit(check = FALSE)`).
- **H6:** `check` removed from `.mml_fit()`, and `verbose` left alone.
- **H7:** `check = FALSE` on the complete-case refit in RMitemRestscore, with a comment that `.check_complete_case_categories()` has checked `cc`.
- **H8:** not changed. Two one-line expressions in separate functions, and passing `ml_model` in would add an argument to a function for no behavioural gain.

Tests: test-category-check 51/51. Full suite (with `R_LIBS` at the installed working tree): 0 failures, 1713 pass, 7 skips, 1 known mirt warning.

## Cutoff empty-category study: script written (2026-11-11, not yet run)

`dev/code-review/cutoff_empty_categories/cutoff_empty_categories.R`, for the open point 2 above. Theory first: the observed data always use every category their own CML fit has (that defines the threshold count), so the current discard rule conditions the null on an event the observed data satisfy. It is then right in principle and only wasteful. Three arms on the same simulated datasets, for RMitemInfitCutoff (default DGP "resample") and RMitemRestscoreCutoff (default "conditional"):
- **A** current rule (fewer iterations than requested);
- **B** allow empty top categories (only `.cml_fit()`'s check);
- **C** current rule topped up to the requested count (at most 10 x the draws).

Outcomes: iterations used, runs with none, familywise false-flag rate under the default p-value flagging and under interval flagging. Design: n 100/200/400, 6 items x 4 categories or 10 x 5, theta mean 0 or -1 (12 cells). It uses the package's own workers (patched only at the discard line) and RMitemInfit/RMitemRestscore. A fidelity check stops the run unless arm A reproduces both *Cutoff() functions exactly.

Container smoke run (R = 4, 40 iterations, 0.9 min): fidelity ok for both. A kept 13-40 of 40 iterations (as low as 7 for restscore), B and C 40. Interval flagging gives high familywise rates in every arm (it has no multiplicity correction, 6-10 items at 95% each), so the p-value rate is the main outcome. Host run: `nohup Rscript dev/code-review/cutoff_empty_categories/cutoff_empty_categories.R 300 10 400 > dev/code-review/cutoff_empty_categories/run.log 2>&1 &`, estimated 3-5 h on 10 cores.
