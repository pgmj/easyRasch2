# Code review: RMdifAnchor() and the LR-only RMdifLR() (uncommitted 1.4.0.9004 diff)

Reviewed 2026-11-16 with `/code-review high` on the uncommitted diff. The diff includes the new `R/dif_anchor.R` (`RMdifAnchor()`), the reduction of `RMdifLR()` to the LR test (`R/dif_lr.R`), and the equal-tailed WY work in `R/utils-multiplicity.R` already reviewed in `restscore-equal-tails.md`. Background for the DIF work is in `dif_lr.md`, from "Redesign decided 2026-11-14" onward.

The findings below come from the review agent. All ten were verified and fixed on 2026-11-16, see "Fixes applied" at the end. The agent reproduced DA1 and DA3 itself. Verify each one before fixing it, as for ET1-ET9. Line numbers refer to the working tree on 2026-11-16.

Tests run by the agent: test-dif_anchor.R, test-dif_lr.R, test-item_restscore_cutoff.R, test-infit_pvalue.R. All pass under C.UTF-8. Under a C locale, six RMdifLR kable expectations fail (DA1).

## Findings

| # | Where | Finding | Status | Reaches users? |
|---|---|---|---|---|
| DA1 | `RMdifLR()` kable, dif_lr.R:266 | Column renamed to `"LR χ²"`. kable fails with "unable to translate ... to native encoding" in a non-UTF-8 locale | reproduced by agent | yes, default output breaks in C/POSIX locales |
| DA2 | threshold-level plot, dif_anchor.R:830 | `rev(sort(unique(Threshold)))` sorts labels as text, so c10 comes between c1 and c2 | plausible | yes, items with 11+ categories |
| DA3 | `seed` validation, dif_anchor.R:280 | Accepts whole numbers above `.Machine$integer.max`. `as.integer(seed)` gives NA and `set.seed()` errors. The settings row has the same `as.integer()` | reproduced by agent (`seed = 3e9`) | edge only |
| DA4 | spacing test, dif_anchor.R:561 | `solve()` on each item's block of S with no guard. A singular block gives a bare base-R error instead of the "small or sparse groups" message | plausible | sparse groups |
| DA5 | `seed` default, dif_anchor.R:237 | Default `seed = 20260113L` and restored RNG state. The reproducibility help topic says `seed = NULL` everywhere and does not list RMdifAnchor(). `.with_local_seed()` does not pin `RNGkind` | partly a known decision, see below | docs, and non-default RNG kinds |
| DA6 | `.wy_attrition_message()`, utils-multiplicity.R:158 | With `floor_mult = 2` and few or no complete iterations, the message can report a floor above 1 ("2/(0+1) = 2") | plausible | edge only |
| DA7 | `.check_group_categories()` message, dif_lr.R:361 | Still opens "Item locations cannot be compared across DIF groups", but RMdifLR() no longer compares locations. The roxygen above the helper says the same | confirmed by reading | yes, misleading error text |
| DA8 | preprocessing, dif_anchor.R:270-347 vs dif_lr.R:141-225 | About 50 lines duplicated (dif_var length check, NA drop, empty-respondent drop, extreme-only category check, `.check_group_categories()`, tryCatch around `.cml_fit()`). The empty-respondent drop is a third copy of `.drop_empty_respondents()` | plausible | no, maintenance |
| DA9 | `isTRUE_vec()`, `.with_local_seed()`, dif_lr.R:379 | Used only by RMdifAnchor() but live in dif_lr.R. `isTRUE_vec()` is not needed, since `Flagged` and `Anchor` are always non-NA logicals | plausible | no, maintenance |
| DA10 | `loc` computation, dif_anchor.R:447 | A second full `.dif_anchor_tests()` call runs every time, but only `output = "list"` and `plot_type = "criterion"` use `loc` | plausible | speed only |

**DA1.** No other table in the package uses non-ASCII column names. `cicc_plot.R` uses `\u` escapes, but only in plot text. Proposal: an ASCII column name such as `LR chisq`, matching the other tables. Check the RMdifLR() roxygen and NEWS for the column name afterwards.

**DA2.** Proposal: order by the numeric threshold index (`res$k`) instead of the label. Add a test with an 11-category item, or check one by hand.

**DA3.** Proposal: check `abs(seed) <= .Machine$integer.max` in the validation and say so in the error message. Check whether other seeded functions have the same gap before choosing the wording.

**DA4.** Proposal: wrap the `solve()` call in `tryCatch()` and stop with the same sparse-groups message the fit and vcov steps give. A near-singular block that inverts without error but gives huge values is a separate question. Decide whether to check the condition number too.

**DA5.** The fixed default is a deliberate exception, approved by the user on 2026-11-15 (`dif_lr.md`, "RMdifAnchor() design approved": a test result should not change between runs). The two remaining issues are:
- The reproducibility help topic says `seed = NULL` everywhere and does not list RMdifAnchor(). Add the function and state the exception.
- `.with_local_seed()` does not pin `RNGkind`, unlike the other seeded functions. Under a non-default generator (for example L'Ecuyer in a mirai daemon) the same `seed` gives different tie-break orderings and possibly different anchors and flags. Proposal: pin the default kind inside the helper and restore it on exit. Ask before changing.

**DA6.** Proposal: `min(1, floor_mult / (b_used + 1))`, or a separate message when no complete iteration is left. This overlaps ET6, which added the helper.

**DA7.** Proposal: reword to say that the group fits would have different thresholds, so the LR test is undefined. Keep the class `easyRasch2_category_error` and the `items` field. Update the roxygen and any test that matches the message text.

**DA8.** Proposal: a shared `.prepare_dif_data(data, dif_var, fn)` returning `data`, `dif_factor`, `n_total` and `has_na`. Results must not change. Run test-dif_lr.R and test-dif_anchor.R afterwards.

**DA9.** Proposal: move `.with_local_seed()` to dif_anchor.R or a utils file, and replace `isTRUE_vec()` with plain indexing. `withr::with_seed()` is an option, but withr is only in Suggests, so keep the helper unless withr moves to Imports.

**DA10.** Proposal: compute `loc` only on the list and criterion paths, or return it from the first `.dif_anchor_tests()` call when `level = "item"`.

## Not in the table

- The pre-knitted `vignettes/easyRasch2.Rmd` still calls `RMdifLR(level = "threshold", output = "ggplot")`, which now errors. This is the known re-knit (see `dif_lr.md` and ET7).
- The agent checked the equal-tailed `.split_pvalues()` logic and how RMitemRestscore() takes the flag direction from it, and found no problem.

## Suggested order

1. DA1, DA7 (user-facing text, small).
2. DA2, DA3, DA4, DA6 (small local fixes).
3. DA5 (docs now, `RNGkind` after asking).
4. DA8, DA9, DA10 (cleanups, optional).
5. Vignette re-knit at release.

## Fixes applied (2026-11-16)

Verification first: DA1, DA2 and DA3 reproduced, DA5, DA7, DA8, DA9 and DA10 confirmed by reading, DA4 and DA6 not reproduced but guarded anyway. DA1 was milder than reported: in a C locale the kable is still produced, with a warning and the header printed as `LR <U+03C7><U+00B2>`.

| # | Fix |
|---|---|
| DA1 | The rename is gone. The kable header is `Chisq`, as in the data.frame and the RMdifAnchor() spacing table. The caption already names the likelihood-ratio test. |
| DA2 | Threshold rows are ordered by the number after `c`. New test with an 11-category item. |
| DA3 | `seed` must lie within `.Machine$integer.max` in absolute value, and the error says so. No other function takes a user seed this way, so nothing else changed. |
| DA4 | `solve()` on each spacing block sits in `tryCatch()` and stops with the sparse-groups message, naming the item by its column number. No condition-number check. |
| DA5 | The reproducibility topic now names RMdifAnchor() as the one exception to `seed = NULL` and says what it pins and restores. `.with_local_seed()` sets `kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection"` and restores the caller's `RNGkind()` and `.Random.seed` on exit (approved by the user, "do all 10"). New test: same result under L'Ecuyer, the kind and state left as they were. |
| DA6 | The floor is capped at 1, and when no iteration is complete the message says the corrected p-values are missing (`.wy_stepdown()` returns NA then). New file test-wy-attrition.R. |
| DA7 | The message opens "Item parameters cannot be compared across DIF groups", which holds for both callers. The helper's roxygen says why. test-dif_lr.R matches the new text. |
| DA8 | New `R/utils-dif.R` with `.prepare_dif_data()` (length check, NA drop, empty respondents through `.drop_empty_respondents()`), `.check_dif_categories()` and `.cml_fit_or_stop()`. Both functions call them. |
| DA9 | `.check_group_categories()` and `.with_local_seed()` moved to `R/utils-dif.R`. `isTRUE_vec()` removed, the four calls use the logical columns directly. |
| DA10 | `loc` is computed only for `output = "list"` and the criterion plot. |

Results unchanged: on PHQ-9 by gender, the list output for item, threshold, spacing, Gini and threshold/Holm, and the RMdifLR() data.frame, are `identical()` to the pre-review code (backups in the session scratchpad).

Tests run: test-dif_anchor (19), test-dif_lr (20), test-wy-attrition (1), test-category-check (14), test-rounding (3), test-item_restscore (9, 1 skipped), test-documentation (1). All pass. Full suite not run.

## Second review round (2026-11-16)

A second `/code-review high` ran on the full working-tree diff after the DA1-DA10 fixes, including the untracked `R/dif_anchor.R` and `R/utils-dif.R`. It found nine issues. DA11-DA18 were verified and fixed on 2026-11-16, see "Second-round fixes applied" at the end. DA19 waits for the re-knit. Line numbers refer to the working tree on 2026-11-16.

The agent reproduced DA11, DA13 and DA14 with a script on `pcmdat2` with planted DIF. It found nothing in the equal-tailed restscore code (`.split_pvalues()`, the doubled floor in `.fdr_min_iterations()` and `.wy_min_iterations()`, the flag direction from the rejecting tail).

Tests run by the agent: test-dif_anchor, test-dif_lr, test-wy-attrition, test-category-check, test-item_restscore_cutoff (`NOT_CRAN=true`). All pass.

| # | Where | Finding | Status | Reaches users? |
|---|---|---|---|---|
| DA11 | `level = "spacing"`, dif_anchor.R:339 | The forward anchor search still runs (100 `psychotools::anchor()` fits), although the spacing test uses no anchor. It can print the tie-break message, an `anchor()` error stops the spacing test, and `settings$Anchor_items` lists anchor items that were never used | reproduced by agent | yes, misleading settings and message, wasted time |
| DA12 | single-step `padj`, dif_anchor.R:570 | The Monte-Carlo `padj` is not floored at the raw `p`. With one tested unit (Gini on 2 items, or one polytomous item at spacing level) `padj < p` about half the time. `.bootstrap_pvalues()` already has this guard (UM12, `pmax(padj, p_marg)`) | plausible | edge only, but looks like an error |
| DA13 | positional calls, dif_lr.R:92 | Removing `level` moved `output` to the 4th position. `RMdifLR(d, g, "auto", "threshold")` fails with the bare `match.arg` error, and `RMdifLR(d, g, "PCM", "item", "dataframe")` with "Unused argument(s) ... <unnamed>". Neither reaches the "use RMdifAnchor()" pointer | reproduced by agent | yes, 1.4.0 code written with positional arguments |
| DA14 | kable NA cells, dif_anchor.R:745 | Numeric columns are turned into strings with NA -> "" only inside `if (any(flag))`. With nothing flagged, the Gini anchor row shows NA in DIF, SE, lower, upper and z. Once a row is flagged, the same cells are blank | reproduced by agent | yes, default output |
| DA15 | p formatting, dif_anchor.R:739 | `formatC(digits = 4)` with no lower bound, so z = -5 prints as p = 0.0000. RMdifLR() prints "< 0.001" | confirmed by reading | yes, default output |
| DA16 | consensus fallback, dif_anchor.R:530 | When no item reaches 50% inclusion, `which(freq == max(freq))` is used. If every run returns an empty anchor (all freq 0) it selects every item, so the linking becomes centring on all items, with no message | plausible | edge only, silent |
| DA17 | `info` list, dif_anchor.R:408 | `info$inclusion`, `info$untested` and `info$seed` are never read. The settings data.frame already reports them | plausible | no, maintenance |
| DA18 | `.round_display()`, dif_lr.R:209 | `p = 4L` has no effect, since `formatC(digits = 3)` overwrites p right after | plausible | no, maintenance |
| DA19 | vignettes/easyRasch2.Rmd:550 | The pre-knitted vignette still calls `RMdifLR(level = "threshold", output = "ggplot")`, which now errors | known | yes, at release |

**DA11.** Proposal: skip `.dif_select_anchor()` and the tie-break message when `level == "spacing"`, and leave the anchor rows out of `settings` (or set them to "none, not used"). Check that the spacing results stay `identical()`.

**DA12.** Proposal: `padj <- pmax(padj, p)` after the adjustment for `correction = "fwer"`, as in `.bootstrap_pvalues()`.

**DA13.** Proposal: when `output` is "item" or "threshold" (or another value of the removed `level`), stop with the same removed-argument message that points to RMdifAnchor(). Add a test with the 1.4.0 positional call.

**DA14.** Proposal: apply the NA -> "" formatting whether or not a row is flagged, so the table looks the same in both cases.

**DA15.** Open: pick one style for both DIF functions ("< 0.001" as in RMdifLR(), or "< 0.0001" to match the 4 digits) and check what the other kable outputs in the package do first (within-package consistency). A shared helper is an option.

**DA16.** Proposal: stop, or warn, when `max(freq) == 0`, and name the fallback in the caption whenever it is used.

**DA17, DA18.** Proposal: drop the unused fields and the `p = 4L` entry. Results must not change.

**DA19.** This is the known re-knit (DIF section above, `dif_lr.md`, ET7). Re-knit before release.

### Suggested order

1. DA13, DA14, DA15 (default output and the migration path).
2. DA11, DA12, DA16 (small local fixes).
3. DA17, DA18 (cleanups).
4. DA19 at release.

## Second-round fixes applied (2026-11-16)

Verification: DA11, DA13, DA14, DA15, DA17 and DA18 confirmed by reading. DA12 is real and not only for one tested unit, since Monte-Carlo noise can put padj below p with several units too. DA16 was not reproduced (psychotools probably never returns an empty anchor) and is guarded anyway.

| # | Fix |
|---|---|
| DA11 | At spacing level no anchor is searched for, unless `plot_type = "criterion"` asks for one. No tie-break message. `settings`: `Anchor = "none (spacing needs no anchor)"`, `Anchor_items` and `N_seeds` NA (Seed kept, it seeds the single-step draws). `criterion` is `NULL` at spacing level. Roxygen updated. |
| DA12 | The `fwer` branch of `adjust()` returns `pmax(p, MC estimate)`, covering item, threshold and spacing. |
| DA13 | `output = "item"` or `"threshold"` (a positional 1.4.0 call) stops before the `...` check with a message naming the removed `level` and RMdifAnchor(). Tests for both positional calls in the review. |
| DA14 | NA cells are blank in every kable, flagged or not. Non-label columns are now right-aligned, since they are strings in both cases. |
| DA15 | User decision: 4 decimals in tables. New `.format_p_table()` in utils-multiplicity.R ("< 0.0001" below that, "" for NA), used by both DIF kables. RMdifLR's table moves from 3 decimals and "< 0.001" to this. `.format_p()` for prose is unchanged. |
| DA16 | Stops with a message pointing to `anchor = "gini"` when every inclusion frequency is 0. Not tested, no data triggers it. |
| DA17 | `info$inclusion`, `info$untested` and `info$seed` removed. |
| DA18 | `p` is formatted before `.round_display()`, so the dead `p = 4L` is gone and test-rounding still sees every numeric column handled. |

Results: on PHQ-9 by gender, the tables for item, threshold, Gini, threshold/Holm and spacing are `identical()` to the round-1 code, and no padj was raised by the floor there. The list output is identical except at spacing level (DA11, intended). RMdifLR()'s data.frame is identical.

Tests: test-dif_anchor (22), test-dif_lr (21), test-wy-attrition (1), test-category-check (14), test-rounding (3), test-documentation (1). All pass. Full suite not run.
