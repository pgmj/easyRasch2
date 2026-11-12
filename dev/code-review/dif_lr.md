# Code review: R/dif_lr.R (RMdifLR)

Started 2026-11-10, after the item_parameters.R unit (see item_parameters.md, "RMdifLR and empty categories" and the frequency simulation, for the background).

## Decision (user, 2026-11-10)

Replace eRm with our own Andersen LR test on psychotools. Stop with an informative message when an item's category is unused in a group, rather than reporting the item as "not tested". Add a proper item-level SE.

## Feasibility check

Own LR, `2 * (sum of group logLik(pcmodel) - logLik(full))`, with df the extra free parameters, against `eRm::LRtest()`. pcmdat2 with 2 and 3 groups, raschdat1 (dichotomous), and pcmdat2 with 5% NA: LR, df and p equal to 4 decimals in all four.

## Applied (uncommitted)

- **Fits:** `.check_response_categories(data)` (full data first, plain message), then the new `.check_group_categories(data, dif_factor)`, then `.cml_fit()` on the full data and on each group. LR, df and p are computed directly. eRm's `LRtest()`, `PCM()/RM()`, `thresholds()`, `extract_lr_locations()`, `extract_overall_locations()`, `one_fit_long()` and `parse_threshold_names()` are gone. The `requireNamespace("eRm")` stop is removed.
- **`.check_group_categories()`:** every category an item uses in the full data must be used in every group, including the top category, which the full-data check allows. psychotools would otherwise fit that group with fewer thresholds. It stops with a classed `easyRasch2_category_error` (`items` field). Message: "Item locations cannot be compared across DIF groups for: I2 (category 2 unused in group 'B'), I3 (no responses in group 'A'). Every response category of an item must be used in every group. Merge a sparse category with an adjacent one, combine groups, or leave the item out. RMplotTile(data, group = dif_var) shows the counts." This replaces eRm's silent dropping, which left such items "Flagged = no".
- **`lr_fit_summary()`:** thresholds from `threshpar(vcov = TRUE)`, sum-zero within each fit. The item location is the mean of its thresholds, and the item SE is `sqrt(w' V w)` over the item's block (`w = 1/m`). Threshold labels keep eRm's style ("c1", "c2"; "1" for RM), so the output columns are unchanged.
- **Caption:** "complete cases" changed to "incomplete responses retained" (CML keeps respondents with missing responses, eRm did too).
- **utils-iarm.R:** the eRm `Rm` branches of `iarm_person_estimates()` were removed (no caller passes eRm fits), and the header now names its only user (RMdimMartinLof). Martin-Löf subscale WLEs and correlations are identical to HEAD.
- **Docs:** RMdifLR roxygen (description, `model`, details with a "Response categories" paragraph, Andersen 1973 reference). NEWS: one Changes bullet (results move for polytomous data) and one bug-fix bullet (silent item dropping, caption).
- **Tests:** test-dif_lr.R drops the 10 `skip_if_not_installed("eRm")` lines its own data never needed, and gains four tests:
  - LR and df equal `eRm::LRtest()` (2 and 3 groups, NA, dichotomous);
  - dichotomous locations and SEs equal eRm, and polytomous thresholds equal eRm up to a constant per group;
  - item location and SE equal `psychotools::itempar(vcov = TRUE)`;
  - the group category stop (top category, lowest category, no responses; class and `items`).
  The source-scan test in test-category-check.R now also forbids `eRm::PCM/RM/LRtest` everywhere (no allow-list).

## Before/after (HEAD vs working tree; pcmdat2 2 and 3 groups, pcmdat2 7% NA, raschdat1, raschdat1 5% NA, phq9)

- LR, df, p: identical in all six.
- Dichotomous: locations, SEs and MaxDiff within 4e-5. Flags identical.
- Polytomous: each group column shifted by a constant (0.29-0.41 logits), equal within 8e-5 after removing it. MaxDiff changed by up to 0.094 logits (3 groups). Flags were the same in these examples but can change.
- Threshold SEs differ by up to 0.047. A threshold's SE depends on the identification (eRm makes the item-category parameters sum to zero, psychotools the thresholds).
- Item SEs: now 0.38-0.67 of the old mean-of-threshold-SEs value. Adjacent thresholds of an item are negatively correlated, so the SE of their mean is smaller. Checked equal to psychotools `itempar()` SEs (difference ~1e-16).

test-dif_lr 58/58, test-category-check 51/51, test-rounding 6/6. Full suite (working tree installed to a scratch library and `R_LIBS` set): 0 failures, 1731 pass, 7 skips, 1 known mirt warning.

## Not yet done for this unit

- The rest of the per-unit review (edge cases such as one tiny group or 2 items, docs vs code, conventions). The `cutoff = 0.5` default is a fixed rule-of-thumb flag, which conflicts with the no-fixed-cutoffs principle. Raise with the user, do not change unasked.
- The ggplot error bars at item level now use the new item SE, so they are narrower.

## /code-review high on the uncommitted changes (2026-11-11)

User-run `/code-review high` over the uncommitted RMdifLR move (dif_lr.R, utils-iarm.R, tests, man, NEWS). Findings DL1-DL3 were confirmed by the review agent with small runs (package loaded from source). DL4-DL9 come from reading the code and have not been re-checked by me in R. Status key as in `utils-theta.md`. No fixes applied.

| # | Line | Finding | Status | Reaches users? |
|---|------|---------|--------|----------------|
| DL1 | 241 | No extreme-score-only category check within DIF groups | confirmed | yes, silent |
| DL2 | 337 | Caption says "incomplete responses retained" when none are | confirmed | yes, caption only |
| DL3 | 452 | Threshold labels c1..cK sort as text for 10+ thresholds | confirmed | yes, output order |
| DL4 | 246 | Fit failures give bare psychotools errors, NaN SEs possible | plausible | yes, unhelpful error |
| DL5 | 441 | No threshold-count vs observed-category guard in `lr_fit_summary()` | plausible | not with current checks |
| DL6 | 410 | Threshold/vcov extraction duplicates `.rasch_fit_cml()` | plausible | no, maintenance |
| DL7 | 244 | `.cml_fit()` uses `model = "pcm"` for dichotomous data | plausible | no, speed only |
| DL8 | utils-theta.R 552, NEWS.md 43 | Docs still describe eRm fits | confirmed (reading) | no, doc only |
| DL9 | test-dif_lr.R 270 | `itempar()` equivalence test only covers equal category counts | plausible | no, test coverage |

### DL1. Extreme-score-only categories within a group (confirmed)

`.check_response_categories()` checks the full data for categories used only by respondents with a zero or perfect score. Nothing does this per DIF group. Example: dichotomous data, 120 rows, groups A and B, where I1 = 1 in group B only for two perfect scorers. Group B gives I1 = 11.77 (SE 165). Through the sum-zero centring the other items move to about -3 (SE 41), and every item is Flagged = TRUE with MaxDiff 2.3 to 12. No warning or error. `RMitemParameters()` warns via `.extreme_only_categories()`. Proposed fix: run the same check per group inside `.check_group_categories()`.

### DL2. Caption wording (confirmed)

`has_na_lr` (line 171) is TRUE when only `dif_var` has NAs, or when the only missing responses are in all-NA rows that are dropped. Complete 200x6 data with 3 NAs in `dif_var` gives "n = 197 of 200 respondents (incomplete responses retained)". Complete data plus one all-NA row gives "n = 199 of 200 (incomplete responses retained)". Proposed fix: compute `anyNA(data)` after the NA `dif_var` rows and empty rows are dropped.

### DL3. Threshold label order (confirmed)

Labels are unpadded "c1".."cK", and `build_wide_table()` and the ggplot sort them as text. A 0-10 NRS item with `level = "threshold"` comes out c1, c10, c2, ..., c9 in the dataframe, kable and colour legend. Proposed fix: sort by the integer index (keeps the eRm-style labels) or zero-pad.

### DL4. Error handling around the fits (plausible)

The old `tryCatch` ("LRtest failed ... inspect with RMplotTile(data, group = dif_var)") was removed and the group fits are not wrapped. A psychotools failure in one group now gives a raw message naming neither the group nor a remedy. `.rasch_fit_cml()` still wraps `.cml_fit()` with advice. A singular Hessian gives NaN SEs with a "NaNs produced" warning from `sqrt(diag(V))` at line 455.

### DL5. Threshold-count guard (plausible)

`lr_fit_summary()` checks that item names match across fits, not that each item's threshold count matches its observed categories, as `.fit_cml_thresholds()` does. The category checks should prevent a reshaped item, so this is a safety net. If psychotools ever returned fewer thresholds for an item in one group, labels and wide-table rows would be mispaired without an error.

### DL6. Duplicated extraction (plausible)

`lr_fit_summary()` repeats the `threshpar(vcov = TRUE)`, `unname(as.matrix(attr(tp, "vcov")))` and `rep.int(seq_along, lengths)` block split from `.rasch_fit_cml()` in item_parameters.R. A shared helper returning `thr_list`, `se_list` and `vcov_thr` would serve both and carry the DL5 guard once.

### DL7. Model choice for dichotomous data (plausible)

All G+1 fits use `.cml_fit()`'s default `model = "pcm"`. `model = "auto"` (raschmodel) gives the same estimates about twice as fast for dichotomous data, and `.fit_cml_thresholds()` already uses it.

### DL8. Stale eRm docs (confirmed by reading)

The `.check_response_categories()` roxygen (utils-theta.R line 552) still says `eRm::LRtest()` drops such items and that the check runs before those fits. The unreleased NEWS bullet at line 43 still lists "(CML, MML or eRm)" fits. No package code fits with eRm any more.

### DL9. itempar test scope (plausible)

The test equating item locations with `psychotools::itempar()` passes because `make_polytomous()` gives all items the same categories. itempar centres item means to sum zero, RMdifLR centres thresholds to sum zero, and these differ by a constant when items have different category counts (e.g. mixed 0-1 and 0-4 items). Either add a mixed-category case that compares up to a constant, or narrow the test title.

### Main-session check and fixes of DL1-DL9 (2026-11-11, uncommitted)

All nine checked against the code. DL1 was also confirmed on real data. The new check stopped my own test case, a random 3-way split of pcmdat2: in one group I3's category 0 was used only by zero scorers. There eRm returned I3.c1 = -13.5 (SE 276), and every other SE in that group was ~55. pcmdat2's I3 has 14 zeros, only 7 from respondents without a zero score.

- **DL1 fixed:** `.check_group_categories()` now also runs `.extreme_only_categories()` per group, skipping items already reported as unused there. It stops with "I1 (category 1 used only by respondents with a zero or perfect score in group 'B')", and the message explains why. RMdifLR also calls `.sparsity_warning(data)` on the full data after the fits, like the other CML functions. Roxygen and NEWS updated (and the eRm comparison corrected: eRm dropped unused-category items, but gave diverging estimates for extreme-only ones).
- **DL2 fixed:** `has_na_lr <- anyNA(data)` after the `dif_var`-NA and empty-row drops. Test covers all three cases.
- **DL3 fixed:** threshold rows are ordered by the integer after "c", and the ggplot `Threshold` is a factor in numeric order. Labels unchanged. Test with 0-10 items.
- **DL4 fixed:** full and group fits are wrapped. A non-category error becomes "CML estimation failed in group 'B' (...). Small or sparse groups are the usual cause; inspect the counts with RMplotTile(data, group = dif_var)." Category errors pass through. Non-finite SEs become NA with one warning naming the group, instead of "NaNs produced".
- **DL5 fixed:** `lr_fit_summary()` checks the per-item threshold counts against the full-data top categories (`n_thr`) and stops with an internal error on mismatch.
- **DL6 not changed:** a shared extraction helper would mean editing `.rasch_fit_cml()` in the already-reviewed item_parameters unit for a 6-line duplication. DL5's guard sits in `lr_fit_summary()`. Can be revisited if a third user appears.
- **DL7 fixed:** `.cml_fit(model = "auto")`. Dichotomous results move by at most 7e-9, and polytomous results are identical.
- **DL8 fixed:** `.check_response_categories()` roxygen (no eRm, points to `.mml_fit()` and `.check_group_categories()`). NEWS "(CML or MML)".
- **DL9 fixed:** the itempar test is now titled for equal category counts. A new test with mixed 0-1/0-2/0-3 items compares locations up to a constant (SEs depend on the identification, so not compared).

Tests: the LR-equality test's 3-group case now leaves I3 out (see DL1). test-dif_lr 66/66, test-category-check 51/51, test-rounding 6/6. DIF baseline vs the previous working tree: identical except dichotomous ≤ 7e-9 (DL7). No baseline split hit DL1. Full suite (working tree installed, `R_LIBS` set): 0 failures, 1739 pass, 7 skips, 1 known mirt warning. Rd regenerated (RMdifLR).

Still open for the unit: other edge cases (e.g. 2 items), and the fixed `cutoff = 0.5` default (raise with the user).

### Docs: cutoff and artificial DIF (2026-11-11, uncommitted)

User decision: keep `cutoff = 0.5` but document its weakness, with references. Drop the references to the old easyRasch package.

- Removed the "replaces the four legacy helpers (RIdif...)" sentence and the "legacy easyRasch package" cell_spec paragraph.
- `@param cutoff`: the default is a convention (Jin et al., 2013; Kim et al., 2017; Pickard et al., 2006), not a justified value.
- Details, "DIF size and the cutoff": the cutoff depends on items, categories and DIF type (Van de Winckel et al., 2022). Judge size against practical consequences (Henninger et al., 2025), ideally by simulating person locations with and without the DIF item split (Johansson et al., 2026, OSF preprint 10.31219/osf.io/3htzc_v2). Threshold-level DIF, MaxDiff uncertainty, and LR test vs size.
- Details, "Real and artificial DIF": real DIF in one item induces opposite artificial DIF in the others (Hagquist & Andrich, 2017, checked against the PDF in articles/DIF). In RMdifLR this happens through the per-group centring. Resolve one item at a time, largest first, and do not overinterpret one run.
- `@references` (all checked against the Johansson et al. reference list or the PDFs) and `@seealso` RMdifGamma/RMdifGammaCutoff and RMdifTree.
- **Found on the way:** the Henninger et al. (2025) reference in dif_tree.R (roxygen and header comment) and dif_partgam.R named Debelak as third author. The Behaviormetrika PDF says Sengewald, M.-A. Corrected (the 2023 EPM paper is correctly Henninger, Debelak & Strobl). Rd regenerated: RMdifLR, RMdifTree, RMdifGammaCutoff. NEWS Docs bullet added.

### Artificial DIF in the other DIF docs (2026-11-11, uncommitted)

User decision: add the real/artificial DIF note to the other DIF functions, and leave `RMdifTree(purification = "none")` as the default.

- **RMdifGamma:** Details paragraph "Real and artificial DIF". Partial gamma conditions on the total score (item included), so real DIF induces opposite artificial DIF (Hagquist & Andrich, 2017). The simulated cutoffs do not guard against this, since the null has no DIF. Resolve one item at a time.
- **RMdifGammaCutoff:** one paragraph saying the null has no DIF, so the cutoffs do not catch artificial DIF, pointing to RMdifGamma.
- **RMdifTree:** Details paragraph. The splits come from a joint parameter test. The per-node MH and partial gamma condition on the total score, so B/C classifications can be artificial. `purification = "iterative"` limits this (Bjorner et al., 1998).
- Hagquist & Andrich (2017) added to the references of all three. NEWS Docs bullet extended. Rd regenerated: RMdifGamma, RMdifGammaCutoff, RMdifTree.

## Edge cases, docs vs code, conventions (2026-11-11, uncommitted)

Script: scratchpad `dif_edge.R`. Fine as they were:
- 2 polytomous items;
- 6 groups;
- numeric, logical and factor-with-unused-level `dif_var`;
- matrix and tibble input;
- item names with spaces, dots and parentheses;
- mixed 0-1 and 0-2 items (item and threshold level);
- `model = "PCM"` on dichotomous data;
- all-NA `dif_var`, an all-NA item, non-integer data;
- kable with `cutoff = NULL` and `sort`, and the threshold ggplot.

A tiny group (n = 8) or a one-person group gives the per-group category message.

Found and fixed:
- **2 dichotomous items failed** ("n > 1L is not TRUE"). This was a DL7 regression: `raschmodel()` errors in its Hessian with two items, while `pcmodel()` gives the same estimates. `.cml_fit(model = "auto")` now uses raschmodel only with 3+ items. With `hessian = FALSE`, as in `.fit_cml_thresholds()`, it already worked. Unreleased, so no NEWS entry.
- **A `dif_var` level named like an output column:** "All" silently replaced the whole-sample column (reshape "first taken" warning only), "Item" errored obscurely, and "MaxDiff" duplicated a column. Now stops naming the clashing levels (Item, Threshold, All, MaxDiff, Flagged, SE_*). In the released version too, so it has a NEWS entry.
- **Duplicated item names** mixed up rows, and 1 item gave a confusing extreme-score message. Both now stop up front. `validate_response_data()` checks neither, so they are local to RMdifLR. A package-wide check belongs to the utils-validation.R unit.
- **Docs vs code:** `@param output` said the data.frame always has `Flagged` (absent with `cutoff = NULL`), and `@param sort` said kable only (it also sorts the data.frame). Both corrected.
- **Conventions:** kable rounding now goes through `.round_display()` (3 dp, unchanged), so test-rounding.R actually covers RMdifLR. Before, it silently saw nothing. Output names and default (`"kable"`), caption `.n_caption()` and `.format_p()` are as elsewhere. `sort` is logical, as in the package-wide open `sort` question.

New tests: 2 dichotomous items with finite SEs, the 1-item error, duplicate names, and reserved levels. DIF baseline identical to before these fixes. test-dif_lr 72/72. Full suite (installed working tree): 0 failures, 1745 pass, 7 skips, 1 known mirt warning.

Remaining test gaps (proposed, not written): bold rendering of flagged rows in the kable, item-level ggplot content (error bars use the item SE), `conf` changing the bar width, and the DL4 fit-failure message (hard to trigger now that the checks run first).

Noted for other units: RMreliability fails on 2 items inside mirt ("Too few degrees of freedom"), for the reliability unit.

**Unit status:** dif_lr done apart from the proposed test gaps.

## Second /code-review high (2026-11-11)

User-run `/code-review high` over the two unpushed commits plus the uncommitted working tree (`git diff @{upstream}`), after the fixes above. The review agent confirmed LR, df and p still equal `eRm::LRtest()` on pcmdat2. DL10 and DL12 were confirmed by the agent with small runs (package loaded from source). DL11 and DL13-DL16 come from reading the code and have not been re-checked in R. No fixes applied.

- **DL10 (confirmed), line 298, reserved group names.** The check covers only a fixed list (Item, Threshold, All, MaxDiff, Flagged, SE_*). A `dif_var` with an empty-string level passes it, and RMdifLR then fails in `build_wide_table()` with "undefined columns selected" (reproduced on pcmdat2). Proposed fix: also stop on `!nzchar(groups)`, or build the wide table on internal keys (g1..gK) and add the labels at the end.
- **DL11 (plausible), line 365, `.sparsity_warning(data)` cannot fire.** A category used only by extreme scorers in the full data is unused or extreme-only in every group, so `.check_group_categories()` stops first. The user gets the per-group message advising "combine groups", which cannot help when the problem is in the full data. Proposed fix: drop the call, or check the full data for extreme-only categories first with the plain full-data message, as `.check_response_categories()` does at line 341.
- **DL12 (confirmed), utils-theta.R line 524, 2-item guard too broad.** The DL7 regression fix sends 2-item data to `pcmodel()` also when `hessian = FALSE`, where `raschmodel()` works (checked). 2-item dichotomous subscales in `RMdimMartinLof()` (`compute_ml_statistic()` per simulation iteration) and `.fit_cml_thresholds()` in the bootstrap loops now use the about 2x slower fit. Proposed condition: `ncol(mat) > 2L || !hessian`.
- **DL13 (plausible), line 536, per-group item maxima in the extreme-only check.** The maxima come from the group's own data. If a group never uses an item's top category, a respondent there can count as perfect without being perfect on the full scale, so other items may be wrongly named as extreme-only in that group. Merging the first item's category would then clear them, so the message can send the user to recode items needlessly. Proposed fix: pass the full-data item maxima to `.extreme_only_categories()`.
- **DL14 (plausible), line 585, duplicated threshold extraction.** Same as DL6. The new threshold-count guard and the NA-SE handling exist only in the RMdifLR copy, so `RMitemParameters()` can still give "NaNs produced" from `sqrt(diag())` and has no alignment guard. Strengthens the case for the shared helper deferred under DL6.
- **DL15 (plausible), line 347, category check repeated.** `fit_lr()` calls `.cml_fit()` with the default `check = TRUE`, so `.check_response_categories()` runs again on the full data and on every group (G+1 extra times) after lines 341-342. Proposed fix: `check = FALSE`, as in the other already-checked refits (item_restscore.R, martin_lof.R).
- **DL16 (plausible), line 704, mixed precision in the kable.** Flagged rows are reformatted with `format(nsmall = 3)`, unflagged rows use `as.character()`, so one column can show **0.200** and 0.2. Proposed fix: format the numeric columns once with `formatC(x, format = "f", digits = 3)` before adding the bold.

### Main-session check and fixes of DL10-DL16 (2026-11-11, uncommitted)

- **DL10 confirmed, fixed:** an empty-string level now stops with the clash message ("... or are empty"). A level with spaces or hyphens works (checked). Test added.
- **DL11 confirmed, fixed:** a full-data extreme-only category gave the per-group message twice ("combine groups" cannot help). RMdifLR now checks the full data first with `.extreme_only_categories()` and stops with the plain full-data message plus "RMdifLR() stops here, since these thresholds cannot be compared across groups" (classed, `items`). The dead `.sparsity_warning()` call is removed. Test added (message has no "in group").
- **DL12 confirmed, fixed:** `.cml_fit(model = "auto")` uses raschmodel when `ncol > 2 || !hessian`. 200 fits of 2 dichotomous items took 0.08 s with raschmodel and 0.13 s with pcmodel. This restores the pre-DL7 route for the Hessian-free 2-item fits (Martin-Löf subscales, `.fit_cml_thresholds()`). Test: class raschmodel without a Hessian, pcmodel with one.
- **DL13 plausible, then confirmed by a constructed case, fixed:** `.extreme_only_categories()` gains an optional `item_max` (default unchanged), and `.check_group_categories()` passes the full-data maxima. Test: group B with I1's top category unused and (1, 2, 2) rows. Judged against the group's own maxima, I2's category 2 is wrongly extreme-only. Against the full-data maxima, only I1 is named.
- **DL14 not changed, user decision:** shared threshold-extraction helper for `.rasch_fit_cml()` and `lr_fit_summary()`, so RMitemParameters also gets the alignment guard and NA-instead-of-NaN SEs. See DL6.
- **DL15 fixed:** `fit_lr()` uses `check = FALSE` (full data and every group checked just before).
- **DL16 not confirmed:** in `render_lr_kable()` the numeric column is formatted once (`format(col, nsmall = 3)`) before the `ifelse()`, so both branches use the formatted string. A flagged table shows 0.020 and 0.325 alike (checked). With nothing flagged the columns stay numeric and kable formats them per column. No change.

DIF baseline identical. test-dif_lr 81/81, test-category-check 51/51, test-utils-theta 20/20. Full suite (installed tree): 0 failures, 1754 pass, 7 skips, 1 known mirt warning. NEWS: empty level added to the RMdifLR bug-fix bullet.

### DL14 applied (user: yes, 2026-11-11, uncommitted)

New `.cml_threshold_summary(fit, data, se, label)` in utils-theta.R, the shared extraction for three callers: `.fit_cml_thresholds()`, `.rasch_fit_cml()` (RMitemParameters) and `lr_fit_summary()` (RMdifLR).
- Thresholds come from `threshpar()`, sum-zero within the fit.
- **Alignment guard**, moved from `.fit_cml_thresholds()`: one entry per column, names equal to `colnames(data)`, and threshold counts equal to the observed maxima. This is new for RMitemParameters.
- **Uncomputable SEs** (non-finite `sqrt(diag(V))`) become NA with one warning naming `label` ("the CML fit", "group 'B'", "the full sample"). Before, RMitemParameters gave NaN plus R's "NaNs produced".
- Returns `thr_list` named by item (psychotools' names kept when the data have no column names, as before), and with `se` also `se_list` and `vcov_thr`.
- `lr_fit_summary(fit, data, dif_group, is_rm)` takes the group's data instead of `item_names` and `n_thr`. Every group uses every full-data category, so the observed maxima equal the full-data ones. The item-mean SE keeps its own NA handling (warning only when the threshold SEs were all fine).
- `.cml_threshold_summary` added to the four mirai export lists that carry `.fit_cml_thresholds` (local_dependence, reliability, reliability_curve, residual_pca).

**Slip on the way:** the first version `unname()`d the `.fit_cml_thresholds()` result. RMitemCatProb indexes it by item name, so the full suite caught it with 8 failures in test-item_cat_prob.R. Names restored.

Verification, final code:
- item/person baseline 60/60 identical;
- DIF baseline identical;
- step-2 baseline identical except the four `diflr` cells, saved before the eRm replacement, so expected.

New test: `.cml_threshold_summary()` matches `threshpar()` (thresholds, SEs, names), returns NULL `se_list` with `se = FALSE`, and stops on mismatched categories or reordered columns. Full suite: 0 failures, 1760 pass, 7 skips, 1 known mirt warning.

## Third /code-review high (2026-11-11)

User-run `/code-review high` over the unpushed commits plus the uncommitted working tree. DL17 and DL18 were reproduced by the review agent with `pkgload::load_all()`. DL19-DL25 come from reading the code. No fixes applied.

- **DL17 (confirmed), utils-theta.R line 577, raschmodel Hessian fails on incomplete data.** The fallback to `pcmodel()` checks only `ncol(mat) > 2`. The `raschmodel()` Hessian fails with "n > 1L is not TRUE" whenever any respondent answered only 1 or 2 items. Reproduced: 300 x 10 dichotomous data with `d[1:5, 3:10] <- NA` and two groups. RMdifLR stops with "CML estimation failed in the full sample (n > 1L is not TRUE). Small or sparse groups are the usual cause", which names the wrong cause. `pcmodel()` fits the same data. `model = "PCM"` does not help (DL20). Checked with grep: dif_lr.R line 363 is the only `.cml_fit()` call with `model = "auto"` and the Hessian, so only RMdifLR is affected. Proposed fix: use `pcmodel()` whenever `hessian = TRUE`, or base the condition on the smallest number of answered items per row.
- **DL18 (confirmed), utils-theta.R line 500, duplicated item names reported as a category mismatch.** The `identical(names(tp), items)` guard in `.cml_threshold_summary()` catches psychotools' silent mispairing, but the message is about response categories. Reproduced: RMitemParameters() on names I1, I1, I2, ..., I9 stops with "The CML fit returned thresholds that do not match the observed response categories". All `.fit_cml_thresholds()` users (RMreliability bootstrap, residual PCA, Q3, cicc, item_cat_prob, Martin-Löf) behave the same. Only RMdifLR has the "must be unique" check. Proposed fix: move the check to `validate_response_data()` (utils-validation.R unit).
- **DL19 (plausible), line 612, weaker threshold-count guard.** The DL5 guard compared each group fit's threshold counts with the full-data `n_thr`. It now compares them with the group's own observed maxima, so a group whose top category psychotools dropped would pass. The item location would then average over fewer thresholds, `build_wide_table()` would fill NA for the missing cN, and MaxDiff/Flagged would not be comparable, with no error. The comment at lines 600-601 still describes the full-data check. Proposed fix: pass the full-data maxima to the guard.
- **DL20 (plausible), line 363, `resolved_model` ignored.** `fit_lr()` always uses `model = "auto"`. The docs (lines 65-66) say every fit is a `pcmodel()` PCM, but dichotomous data with 3+ items are fitted with `raschmodel()`. The documented estimator and the `lr_test$model` label do not match the fit.
- **DL21 (plausible), line 243, dead branch.** `is.null(colnames(data))` cannot be TRUE after `as.data.frame()` at line 239, and its message would read "Duplicated: .". Reduce to `anyDuplicated(colnames(data)) > 0L`.
- **DL22 (plausible), line 365, dead handler.** `fit_lr()` still passes through `easyRasch2_category_error`, but `.cml_fit(check = FALSE)` can no longer raise it. The handler suggests the check runs inside the fit.
- **DL23 (plausible), lines 342 and 582, utils-theta.R line 667, repeated condition constructor.** The classed `easyRasch2_category_error` condition is built by hand three times. A `.category_error(msg, items)` helper would keep the class and fields in one place.
- **DL24 (plausible), line 391, repeated subsetting.** Each group is subset with `data[dif_factor == g, , drop = FALSE]` for the fit (line 377), again for `lr_fit_summary()`, and again in `.check_group_categories()`. One `split(data, dif_factor)` could serve all three.
- **DL25 (plausible), NEWS.md line 52, miscount.** The RMdifLR bug-fix bullet lists three defects (an "All" level, an empty level, duplicated item names) and then says "Both now stop". Change to "All three now stop".

### Main-session check and fixes of DL17-DL25 (2026-11-11, uncommitted)

- **DL17 confirmed, fixed:** reproduced with 300 x 10 dichotomous data, `d[1:5, 3:10] <- NA` ("CML estimation failed in the full sample (n > 1L is not TRUE)"). `raschmodel(hessian = FALSE)` works on the same data with identical thresholds. The `.cml_fit(model = "auto")` rule is now "raschmodel only without a Hessian", replacing the DL7/DL12 column-count rules. Without a Hessian it works with 2 items and with incomplete data, and it stays the fast route in the bootstrap and Martin-Löf loops. RMdifLR calls `.cml_fit(d, check = FALSE)` (pcmodel) explicitly. Tests: the reproduced case, and `.cml_fit` class by Hessian for 10 items.
- **DL18 confirmed, fixed in part:** `.cml_threshold_summary()` now stops on duplicated names with "Item names (column names of `data`) must be unique. Duplicated: I1." before the category guard. This covers RMitemParameters and every `.fit_cml_thresholds()` user (checked with RMitemParameters). Moving the check up front into `validate_response_data()` stays with the utils-validation.R unit.
- **DL19 partly right, fixed:** the guard compares the fit with the data's observed maxima, so a threshold psychotools dropped from data that has the category would still be caught. It missed only a group whose data lack the top category, which `.check_group_categories()` blocks. `.cml_threshold_summary()` gains `max_cat`, and `lr_fit_summary()` passes the full-data maxima, so the guard no longer relies on the group check. Comment and roxygen updated. Test: a wrong `max_cat` stops.
- **DL20 resolved by DL17:** every RMdifLR fit is now `pcmodel()`, as the docs say. The `lr_test$model` label (RM/PCM) names the model, and the two coincide for dichotomous data.
- **DL21 fixed:** dead `is.null(colnames(data))` branch removed.
- **DL22 fixed:** `fit_lr()` no longer passes `easyRasch2_category_error` through (cannot occur with `check = FALSE`). Comment says so.
- **DL23 fixed:** new `.category_error(msg, items)` constructor in utils-theta.R, used by `.check_response_categories()`, `.check_group_categories()` and the RMdifLR full-data extreme-only stop.
- **DL24 not changed:** three group subsets per group cost little next to the fits and read clearly.
- **DL25 fixed:** NEWS "All three now stop", reflowed.

Baselines vs the DL14 runs: item/person 60/60 identical. The DIF and step-2 sets differ only in the dichotomous RMdifLR cells (max 7.1e-9, back to pcmodel). Tests: test-dif_lr 85/85, test-utils-theta 28/28. Full suite (installed tree): 0 failures, 1766 pass, 7 skips, 1 known mirt warning.

## Remaining test gaps written (2026-11-11, uncommitted)

Three new tests at the end of test-dif_lr.R (polytomous, n = 300, 5 items, two groups):
- **Kable bold:** with `cutoff` at the median MaxDiff, every data cell of a flagged row is bold (`**yes**`, the item name not), unflagged rows have no bold and "no", and the caption gives the bold note. `cutoff = NULL` gives no Flagged column, no bold and no note. `sort = TRUE` orders the kable rows by MaxDiff.
- **Item ggplot:** the group locations and SEs in the plot data equal the dataframe output, the error-bar half-widths equal qnorm(0.975) times the item SE, `conf = 0.80` scales them by qnorm(0.90)/qnorm(0.975) and says "80% CI", and the diamonds sit at the All locations.
- **DL4 fit failure:** `.cml_fit()` mocked with `local_mocked_bindings()`. A failure in the full fit gives "CML estimation failed in the full sample (boom). Small or sparse groups ...", and one in the first group fit names group 'a'.

test-dif_lr.R passes. Unit status: dif_lr done.
