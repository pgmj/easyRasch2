# Code review: R/utils-validation.R

Started 2026-11-11 (Tier 1 shared helpers). 460 lines, 13 helpers. Call sites (excluding definitions and roxygen):
- `validate_response_data` 40 (29 files, all model-fitting functions; the data plots do not call it)
- `.n_caption` 37
- `.drop_empty_respondents` 19
- `.round_display` 18
- `.check_cutoff_sample` 14
- `.attrition_clause` 10
- `.all_sims_failed_message` 6
- `.extreme_only_categories` 5
- `.extreme_only_message` 4
- `.notify_dropped_resamples` 4
- `.validate_filename` 3
- `.mi_completed_to_numeric` 2
- `.write_output_csv` 2

Status key as in utils-theta.md. Nothing fixed yet.

## Findings

| # | Helper | Finding | Status | Reaches users? |
|---|---|---|---|---|
| UV1 | `validate_response_data` | No check for at least 2 items or unique item names | confirmed | yes, misleading errors |
| UV2 | `validate_response_data` | Character or factor columns give "non-numeric argument to mathematical function" | confirmed | yes, obscure error |
| UV3 | `validate_response_data` | `Inf` passes. RMitemParameters then says "CML estimation failed (attempt to make a table with >= 2^31 elements) ... try estimator = "MML"" | confirmed | yes, wrong advice |
| UV4 | `.mi_completed_to_numeric` | Factor levels that are labels ("never", "often") silently become NA, with only R's coercion warning | confirmed | yes, obscure downstream error |
| UV5 | `.drop_empty_respondents` | Doc: "one-time message" is per call. It names only pcmodel's error, not raschmodel's crash | confirmed (reading) | no, doc only |
| UV6 | several | No direct tests for `.drop_empty_respondents`, `.n_caption`, `.mi_completed_to_numeric`, `.validate_filename`, `.write_output_csv`, `.all_sims_failed_message` | confirmed | no, coverage |

**UV1.** Carried over from the dif_lr unit (DL18).
- **Duplicated item names:** they now stop with a clear message in `.cml_threshold_summary()`, but only after the model is fitted, and only on paths that use it. RMdifLR has its own up-front check.
- **One item:** RMdifLR checks for this. Elsewhere, a single item reaches the fitters, where every respondent is extreme.
- **Proposal:** add both checks to `validate_response_data()`, so all 29 callers stop up front. Every caller fits a model, so no legitimate single-item use is lost. RMdifLR's local copies can then go.

**UV2.** `as.matrix()` on a data.frame with any character or factor column turns everything into character, and `floor()` then fails.
- **Proposal:** stop with "`data` must contain numeric item responses. Non-numeric columns: a, b." For factors, add "convert factors with as.numeric(as.character(x)) if their levels are the response codes". Logical columns pass and are treated as 0/1, which is reasonable, so leave them.

**UV3.** `Inf == floor(Inf)` is TRUE, so `Inf` passes the integer check.
- **Proposal:** stop on non-finite values ("`data` contains infinite values").

**UV4.** RMitemInfitMI and RMitemInfitCutoffMI pass each completed imputation through this helper. `as.numeric(as.character(x))` is right for factors whose levels are the codes ("0", "1", ...). Labelled levels become NA, and the item then fails later as "no responses" or an empty category, which points the user at the wrong problem.
- **Proposal:** stop when a factor's levels are not all numeric, naming the columns, and say that the levels must be the 0-based response codes.

**UV5.** Wording only: "emits a message when any are dropped", and "pcmodel() errors on them and raschmodel() can crash the R session" (as in the `.cml_fit()` doc).

**UV6.** Proposed tests (not written):
- `.drop_empty_respondents` keeps row names and order, and works on matrix and data.frame input.
- `.n_caption` covers all four shapes (with or without "of Y", with or without qualifiers) and the noun argument.
- `.mi_completed_to_numeric` handles code levels, label levels (after UV4: the stop) and mixed factor/numeric columns.
- `.validate_filename` and `.write_output_csv` cover the bad-input cases and a round trip in a temp file.
- `.all_sims_failed_message` gives the extreme-only branch and the generic branch.

## Checked, no issue

- **`.extreme_only_categories`:** the extreme flag uses only answered items. Unused categories are skipped (handled by `.check_response_categories()`). The new `item_max` argument (DL13) defaults to the old behaviour.
- **`.notify_dropped_resamples`:** parallel failures arrive as the string "mirai_error" (reliability.R:742, reliability_curve.R:762), so they are counted as "another reason" correctly.
- **`.n_caption`, `.attrition_clause`, `.round_display`, `.validate_filename`, `.write_output_csv`, `.check_cutoff_sample`:** they do what their docs say. `.round_display` also accepts a list of digits (RMdifLR passes one).
- **`.drop_empty_respondents` in the person-level functions:** RMpersonParameters and RMpersonFit return fewer rows than the input when respondents are empty. This is documented in `@return` ("respondents with no responses are dropped with a message"), and the row names (RMpersonParameters) and `id` (RMpersonFit) keep the original row identity. Checked: 40 rows with one empty, 39 returned, row names 1, 2, 4, ...

## Related, other units

- `.attrition_clause()` says simulated datasets are lost "usually because an item ended up with an unused response category". That ties to the open Cutoff any-empty-category validation (Cutoff units).

## Fixes applied (user: fix and write tests, 2026-11-11, uncommitted)

- **UV1:** `validate_response_data()` stops on fewer than 2 items ("`data` must have at least 2 items (columns). Found 1.") and on duplicated names ("Item names (column names of `data`) must be unique. Duplicated: x."). A matrix without column names still passes.
  - Removed the now-unreachable local checks in RMdifLR, RMitemICCPlot, RMitemCatProb and RMitemHierarchy.
  - Kept RMreliabilityCurve's: it runs after `.align_items()` can reduce the items to those in `item_params`.
  - The data plots (RMplotTile, RMplotBar, ...) do not call the validator and keep their own "at least 2 columns" checks.
  - test-item_cat_prob.R now expects "at least 2 items".
  - The `.cml_threshold_summary()` duplicate check (DL18) stays as a guard for internal callers.
- **UV2:** non-numeric columns stop with "`data` must contain numeric item responses. Non-numeric column(s): a.", plus a conversion hint when any are factors. Logical columns pass, read as 0/1.
- **UV3:** non-finite values stop ("`data` contains infinite values"). NaN counts as missing, as before.
- **UV4:** `.mi_completed_to_numeric()` stops when a factor's levels are not all numeric, naming the columns and suggesting `factor(x, levels = 0:3)` before imputing. Code-level factors, ordered ones included, convert as before.
- **UV5:** `.drop_empty_respondents()` doc corrected (per-call message, raschmodel crash, keeps order and row names).
- **UV6:** new tests/testthat/test-utils-validation.R, 51 expectations. It covers:
  - `validate_response_data`: every check, including the old ones;
  - `.drop_empty_respondents`: matrix and data.frame, row names, the silent no-op;
  - `.n_caption`: all four shapes plus the noun;
  - `.mi_completed_to_numeric`: code levels, ordered, labels stop;
  - `.validate_filename` and `.write_output_csv`: bad inputs, CSV round trip with and without row names;
  - `.all_sims_failed_message`: both branches.

NEWS: one bug-fix bullet. Baselines (item/person, step 2, DIF) identical. Full suite (installed tree): 0 failures, 0 errors, 1817 pass, 7 skips, 1 known mirt warning. Earlier suite totals in these logs counted only failures; each of those runs also printed an empty failed-test section, so none hid errors.

## Review of the fixes (/code-review high, 2026-11-11)

The user ran `/code-review high` on the uncommitted UV1-UV6 changes. It reported 7 findings. FR1-FR3 were reproduced in R with the package loaded. FR4-FR7 come from reading the code and were not run. No separate verification pass.

| # | Where | Finding | Status | Reaches users? |
|---|---|---|---|---|
| FR1 | utils-validation.R:23 | The 2+ items check breaks single-item scoring with fixed `item_params` | confirmed | yes, regression |
| FR2 | utils-validation.R:23 | The item count is checked before `.align_items()`, so the check is inconsistent in item_params paths | confirmed | yes, inconsistent |
| FR3 | targeting_plot.R:192 | RMtargeting calls `rowSums(data)` before validating, so UV2's message never shows | confirmed | yes, obscure error |
| FR4 | conditional_infit_mi.R:230, infit_cutoff_mi.R:192 | MI functions validate inside the per-imputation tryCatch | from reading | yes, obscure error |
| FR5 | NEWS.md:45 | The bug-fix bullet overstates coverage and omits the FR1 regression | from reading | yes, docs |
| FR6 | utils-validation.R:32 | The duplicate-name stop is a copy of the one in `.cml_threshold_summary()` (utils-theta.R:503-510) | from reading | no, maintenance |
| FR7 | item_cat_prob.R:193 | A stale comment explains the removed local 2-item check | from reading | no, code only |

**FR1.** RMpersonParameters, RMretestSD and RMpersonChange call `validate_response_data()` in their `item_params` paths, where no model is fitted. `RMpersonParameters(d[, 1, drop = FALSE], item_params = ip)` used to return WLEs and now stops with "`data` must have at least 2 items (columns). Found 1." The roxygen line "Every caller fits a Rasch model" is wrong for these paths. Single-item scoring with fixed parameters had no test, which is why the full suite passed.

**FR2.** In the item_params paths, `.align_items()` decides the item set, and it runs after the validator. With 5 data columns and `item_params = ip[1, , drop = FALSE]`, RMpersonParameters returns one-item WLEs without an error, while the 1-column call in FR1 is refused. RMreliabilityCurve keeps its own check after alignment for this reason.

**Decision needed (user) for FR1 and FR2.** Should single-item scoring with known parameters be allowed?
- Allow: skip the item-count check in the validator for item_params paths and decide after `.align_items()`.
- Refuse: keep it, also check after `.align_items()` so both cases agree, and state the change in NEWS.

**FR3.** The default `bins` uses `rowSums(data)` before `validate_response_data()`. A character column gives "'x' must be numeric or complex". Proposal: validate before computing `bins`.

**FR4.** RMitemInfitMI validates each completed dataset inside its tryCatch, and RMitemInfitCutoffMI does the same through RMitemInfitCutoff. A single item, duplicated names or a non-numeric column fail identically in every imputation. The user gets m warnings "Model fitting failed for imputation i: ..." and then "Model fitting failed for all m imputed datasets. Check your data." Proposal: validate the first completed dataset once, before the loop.

**FR5.** The bullet says every model-fitting function now stops up front on non-numeric columns or a single item. FR3 and FR4 are exceptions. Depending on the FR1/FR2 decision, the bullet should also say that single-item scoring with fixed parameters now errors. Proposal: rewrite once FR1-FR4 are settled.

**FR6.** The message and logic are now in two places, each enforced by its own tests (test-utils-validation.R, test-utils-theta.R:108, test-dif_lr.R:371). Proposal: a small `.check_unique_item_names(items)` helper called from both.

**FR7.** Three comment lines plus a parenthetical describe a check that no longer exists. Proposal: replace them with one line next to the `validate_response_data()` call, as in cicc_plot.R and item_hierarchy.R.

**Tests to add with the fixes:** single-item scoring with `item_params` for RMpersonParameters, RMretestSD and RMpersonChange (whichever way FR1 is decided), the FR2 alignment case, RMtargeting with a character column, and an MI call with a single item.

### Fixes applied (user: minimum 2 items is fine, fix single-item scoring only if simple, 2026-11-11, uncommitted)

- **FR1, FR2 (refuse):** single-item scoring with `item_params` stays an error, now consistently. `.align_items()` stops after choosing the items: "`item_params` must cover at least 2 items in `data`. Found 1." This covers RMpersonParameters, RMreliabilityCurve, RMretestSD and RMpersonChange (and RMpersonFit internally). RMreliabilityCurve's own post-alignment check was removed as redundant. The validator roxygen no longer says every caller fits a model.
- **FR3:** RMtargeting validates before the default `bins` uses `rowSums()`.
- **FR4:** RMitemInfitMI and RMitemInfitCutoffMI check the first completed dataset once before the loop, with the new `.validate_response_structure()` (class, 2+ items, unique names, numeric columns), which `validate_response_data()` now calls first. The value checks stay per imputation: a first attempt that ran the full validator up front broke three test-infit_mi.R tests, where a negative value in only some imputations must fail only those.
- **FR5:** NEWS bullet rewritten ("fits or scores a model", MI checked once), plus a second bullet stating that single-item scoring with `item_params` now errors.
- **FR6:** new `.check_unique_item_names()` in utils-validation.R, called from `validate_response_data()` and `.cml_threshold_summary()`.
- **FR7:** item_cat_prob.R comment reduced to "# also checks for 2+ items" on the validator line.
- **Tests** (test-utils-validation.R): one data column and one `item_params` item for RMpersonParameters, RMreliabilityCurve, RMretestSD, RMpersonChange, 2 items still score; RMtargeting with a character column; both MI functions with a single-item `mids`, error with no per-imputation warnings. RMpersonFit has no `item_params` argument, so its alignment path is internal only.

Full suite after the FR fixes (installed tree): 0 failures, 0 errors, 1828 pass, 7 skips, 1 known mirt warning.

## Second review of the fixes (/code-review high, 2026-11-11)

The user ran `/code-review high` again, on the UV + FR changes. It reported 9 findings. SR1, SR3, SR4 and SR5 were reproduced in R with load_all() (4-item dichotomous data, n = 400). The others come from reading the code, and SR2 and SR8 were checked against NEWS.md.

| # | Where | Finding | Status | Reaches users? |
|---|---|---|---|---|
| SR1 | reliability_curve.R:288 | Removing RMreliabilityCurve's own 2-item check unguards `items =`: `items = "I1"` returns a curve (reliability ~0.04) with no error | confirmed (run) | yes, regression from FR1 fix |
| SR2 | NEWS.md:45 | First bullet says every function stops on a single item, but SR1 does not | confirmed (reading) | yes, docs |
| SR3 | reliability_curve.R:523 | Pre-existing: `.curve_select_items()` accepts repeats; `items = c(1, 1)` gives sem ~1700, reliability ~1e-5 | confirmed (run) | yes, nonsense output |
| SR4 | person_parameters.R:363 | `.align_items()` 2-item guard bypassed by duplicated names in list `item_params`: `list(I1 = t, I1 = t)` scores I1 twice (n_answered = 2) | confirmed (run) | yes, defeats FR1/FR2 |
| SR5 | person_parameters.R:143 | With `item_params`, the validator checks all columns, so duplicated names in unused columns now stop (`cbind(d, x = 1, x = 2)`) | confirmed (run) | yes, new refusal |
| SR6 | conditional_infit_mi.R:226, infit_cutoff_mi.R:181 | Up-front MI structure check completes and converts imputation 1 twice; 5-line block copied in two files | from reading | no, maintenance |
| SR7 | utils-validation.R:410 | `.mi_completed_to_numeric()` walks the columns twice; `labelled` could run over `data[factor_cols]` only | from reading | no, minor |
| SR8 | NEWS.md:53 | Second bullet names RMreliabilityCurve, which already refused one aligned item (old post-alignment check); it also never returned person locations | confirmed (diff) | yes, docs |
| SR9 | NEWS.md:54 | The two new bullets lack the blank line every other bullet has; the item_params refusal removes behaviour and arguably belongs under Changes | confirmed (reading) | yes, docs |

All nine confirmed on reading the code (SR1, SR3-SR5 also by the reviewer's runs).

### Fixes applied (2026-11-11, uncommitted)

- **SR1, SR3:** `.curve_select_items()` stops when `items` repeats an item ("`items` must not repeat an item.") or selects fewer than 2 ("`items` must select at least 2 items. Found 1."). The `items` roxygen says "at least 2 distinct items" (Rd regenerated).
- **SR4:** `.align_items()` calls `.check_unique_item_names(names(thr_list), "names of `item_params`")` first. The helper gained a `where` argument for the message.
- **SR5: left as is.** With `item_params`, duplicated names in unused data columns now stop. Duplicated column names make the data ambiguous, the old non-numeric refusal in these paths already applied to unused columns, and it is a one-line change if the user wants unused columns ignored. Flagged to the user.
- **SR6:** new `.mi_validate_structure(mids_object)` in utils-validation.R. Both MI functions call it in one line. It still completes imputation 1 twice in RMitemInfitMI (negligible).
- **SR7:** `.mi_completed_to_numeric()` checks levels only over the factor columns.
- **SR2, SR8, SR9 (NEWS):** blank lines restored. The `item_params` refusal moved to Changes and names RMpersonParameters, RMretestSD and RMpersonChange (not RMreliabilityCurve, which already refused), "used to run" instead of "returned person locations", and it mentions duplicated `item_params` names. A separate Bug-fixes bullet covers the RMreliabilityCurve `items` fixes. The general bullet is accurate again now that `items =` is covered.
- **Tests** (test-utils-validation.R): `items = "I1"`, `c(1, 1)`, a repeated name, a bad type; duplicated `item_params` names for RMpersonParameters and RMreliabilityCurve.

Full suite after the SR fixes (installed tree): 0 failures, 0 errors, 1834 pass, 7 skips, 1 known mirt warning.

## Third review of the fixes (/code-review high, 2026-11-11)

The user ran `/code-review high` on the working tree after the SR fixes. It reported 6 findings. TR1-TR4 were reproduced in R with load_all() (3-item dichotomous data, columns a, b, c). TR5 was checked with grep and TR6 against the roxygen.

| # | Where | Finding | Status | Reaches users? |
|---|---|---|---|---|
| TR1 | person_parameters.R:364 | Pre-existing (same in HEAD): names in `item_params` that partly match the data, with an equal item count, fall back to positional pairing. `list(b, a, z)` on columns a, b, c scores a with b's thresholds, silently | confirmed (run) | yes, wrong person locations |
| TR2 | reliability_curve.R:535 | Repeat check compares raw numerics; `items = c(1, 1.5)` passes and selects a and a.1 | confirmed (run) | yes, double-counted item |
| TR3 | reliability_curve.R:524 | `items = c(1, NA)` gives R's "missing value where TRUE/FALSE needed" | confirmed (run) | yes, obscure error |
| TR4 | person_parameters.R:363 | `list(a = t, t2, t3)` (some names empty) now stops with "Duplicated: ." Before, it fell through to positional pairing (TR1). An all-unnamed list is refused by `.coerce_item_params()` already | confirmed (run) | yes, misleading message |
| TR5 | conditional_infit_mi.R:360 | RMitemInfitMI completes imputation 1 three times (structure check, loop i = 1, `item_names` at :360) | confirmed (grep) | no, efficiency |
| TR6 | utils-validation.R:410 | A factor covariate left in the mids (e.g. sex) stops with "levels are not response codes". The reviewer said the roxygen does not say item columns only; it does (both MI functions, @param mids_object). Only the message could add "or is not an item column" | partly confirmed | yes, minor wording |

All six confirmed on reading the code (TR6: the roxygen already says item columns only, so only the message changed).

### Fixes applied (2026-11-11, uncommitted)

- **TR1 (pre-existing, results move for affected calls):** `.align_items()` stops when some, but not all, `item_params` names are columns of `data` ("Names in `item_params` not found in `data`: z. The names must match the item columns of `data`."). Names that match no column still pair by position when the counts are equal (unchanged). RMpersonParameters `@param item_params` now states the rule (the other functions point to it).
- **TR4:** empty or NA names in a list `item_params` stop first: "Every element of `item_params` must be named by item. Unnamed: element(s) 2, 3."
- **TR2, TR3:** `.curve_select_items()` stops on `NA` ("`items` must not contain NA.") and on fractional indices ("`items` must be whole column indices.").
- **TR5:** `.mi_validate_structure()` returns the item names invisibly. RMitemInfitMI uses them instead of completing imputation 1 a third time.
- **TR6:** the label-levels message adds "Remove columns that are not items, such as covariates, from the mids object."
- **NEWS:** the RMreliabilityCurve `items` bullet adds NA and fractional indices. A new Bug-fixes bullet covers TR1.
- **Tests** (test-utils-validation.R): pairing by name in reverse order, by position with no matching names (same locations), partial overlap stops, empty names stop; `items` with NA (numeric and character) and 1.5; the covariate hint.

Full suite after the TR fixes (installed tree): 0 failures, 0 errors, 1842 pass, 7 skips, 1 known mirt warning.

**TR1 follow-up (user: "Using position seems risky", 2026-11-11):** the positional fallback is removed. `.align_items()` now pairs by name only: any `item_params` name that is not a column of `data` stops ("Names in `item_params` not found in `data`: x, y, z. ..."), and a list with no names stops ("`item_params` must be named by item."). The old "Number of items ... does not match the number of columns" branch is gone. test-reliability_curve.R:182 now expects the new message, the TR1 test expects the all-mismatch stop, and the roxygen and NEWS say name-only.

Full suite after the TR1 follow-up (installed tree): 0 failures, 0 errors, 1842 pass, 7 skips, 1 known mirt warning.
