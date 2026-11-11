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
