# Code review: R/utils-multiplicity.R

Started 2026-11-11 (Tier 1 shared helpers). 568 lines. The first pass was `/code-review high` on the whole file, since there was no diff. Its ten findings (UM1-UM10) were verified in the main session by sourcing the file in R (`LANG=C.UTF-8`), using the probe script in the session scratchpad, and by reading the callers.

`.bootstrap_pvalues()` callers: conditional_infit, item_restscore, local_dependence (Q3), ld_partgam, dif_partgam, cfa_cutoff (fit indices and loadings), residual_pca (one column, `correction = "none"`).

Status key as in utils-theta.md. Nothing fixed yet.

## Findings

| # | Helper | Finding | Status | Reaches users? |
|---|---|---|---|---|
| UM1 | `.bootstrap_pvalues` | NA simulated values are left out of the exceedance count but stay in the `B + 1` denominator, so the p-value is too small | confirmed | yes: CFA (non-finite RMSEA/CFI set to NA), partial gamma (NA gamma when a pair has no concordant or discordant pairs) |
| UM2 | `.wy_stepdown` | Same NA mechanism in the step-down: `max(na.rm = TRUE)` shrinks the max-T, and all-NA rows give `-Inf` and a warning on every row | confirmed | yes, same paths as UM1, plus warning spam |
| UM3 | `.bootstrap_pvalues` | A column whose simulated values are constant gets the minimum p (flagged) under `"joint"`, while `"simulated"` and the docs give NA | confirmed | rarely; it needs a constant simulated column |
| UM4 | `.warn_fdr_floor` | The BY message quotes `alpha/m` and `m/(s(B+1))` without the `c(m)` factor, which contradicts its own `need` | confirmed | yes, wrong number in a warning (fdr_by only) |
| UM5 | `.warn_fdr_floor` + roxygen of `.fdr_min_iterations` | Westfall-Young is said to need `1/alpha - 1 = 19` iterations. Flagging is a strict `padj < alpha`, and at B = 19 `padj >= .05`, so 20 are needed | confirmed | yes, wrong number in a warning |
| UM6 | `.wy_stepdown` | O(k² B): a matrix copy and a row-wise `apply(max)` at every step | confirmed | speed only: k = 1225 (50-item Q3), B = 1000 takes 4.3 s here |
| UM7 | `.fmt_width` | 5 decimals show the Sidak width as "1" at k = 19900 (200 items, Q3) | confirmed | edge only; k = 1225 shows .99996 |
| UM8 | `.p_stars` | Strict `<`, while iarm's `symnum()` is right-closed, so p = .05 gives "." here and "*" in iarm | confirmed | cosmetic: input is the asymptotic `padj_bh`, which almost never lands on a cutpoint |
| UM9 | `.notify_low_iterations`, `.iteration_note` | Both say error rates are controlled at any number of iterations whatever the correction, including `"none"` and FDR | confirmed (reading) | yes, wording |
| UM10 | `.bootstrap_pvalues` | `names(observed)` is not validated. Unnamed input returns an empty data.frame silently, and duplicate names reuse one column | confirmed (reading) | no, every current caller names `observed` |

**UM1.** In the probe, setting 100 of 200 simulated rows of item `a` to NA moved p from .0995 to .0498. Dropping those rows instead (B = 100) gives .0990, which is the correct answer.
- **Where it reaches users:** RMcfaCutoff deliberately sets Inf RMSEA to NA (cfa_cutoff.R:688). The partial gamma helpers return `NA_real_` for an undefined gamma (ld_partgam.R:1248, dif_partgam.R:1063), and those NA values feed `sim_mat` through `tapply`. The Q3 and infit paths should not produce NA, because the Cutoff validation discards bad iterations, but that has not been checked.
- **The Inf case is a second question.** An Inf simulated RMSEA is more extreme than any observed value, so counting it as an exceedance is arguably more correct than dropping it.
- **Proposal:** use a per-column denominator, `(1 + sum(x >= obs, na.rm = TRUE)) / (1 + sum(!is.na(x)))`, and compute the joint moments on the same non-NA values. For CFA, keep Inf as Inf in the upper tail rather than turning it into NA (needs a decision). Alternatively, if NA cells indicate a failed iteration, drop whole rows. With per-column NA the WY maximum is not defined consistently, so dropping rows is the cleaner choice for the FWER path.

**UM2.** The probe (two columns, alternating NA) gave 100 "no non-missing arguments to max" warnings. The fix follows from UM1. If rows containing NA are dropped before the step-down, `max_star` is always finite and B is the number of rows actually used.

**UM3.** Simulated column constant at 1, observed 1.01: `"joint"` gives p = padj = .005, `"simulated"` gives NA. Under `"joint"`, the SD guard on line 52 sees the observed value's contribution, so it no longer catches the degenerate case.
- **Proposal:** apply the zero-variance guard to the simulated values alone (`sd(sim[, i])`), and use the joint moments only for studentising. That matches the docs (lines 31-32).

**UM4.** For m = 10 with BY: the message says a raw p below .005 is enough, but BY needs .005 / 2.93 = .0017, which is what `need = 585` reflects.
- **Proposal:** multiply both quoted numbers by `penalty` (the same `c(m)` as `.fdr_min_iterations`).

**UM5.** At B = 19 the probe gives padj = .05 for an extreme item, which is never `< .05`.
- **Proposal:** change line 499 to `ceiling(1 / alpha)` (= 20), and change the roxygen text on line 411 to "`1 / alpha` iterations".

**UM6.** The successive maxima are a reverse cumulative max over `ord`: run `running <- pmax(running, stat_sim[, ord[j]])` from the last-ranked comparison back to the first, which is O(kB) and gives identical results. It only matters for large Q3 / LD-gamma families, so it is optional. If done, combine it with the UM1/UM2 NA handling and check the results are identical against the current version.

**UM7.** **Proposal:** format with enough significant digits to show the first non-9 digit, for example `sprintf("%.*f", max(5, ceiling(-log10(1 - w)) + 1), w)`.

**UM8.** Making `.p_stars` right-closed (`<=`) would match iarm, but a p of exactly .05 would then get "*" without being flagged under the strict `padj < alpha` rule. Since the input is an asymptotic p, this is cosmetic.
- **Proposal:** keep the strict rule and change the comment from "reproduces" to "uses the same cutpoints as iarm".

**UM9.** RMlocdepQ3 etc. call `.notify_low_iterations` whatever the `correction`. With `"none"` there is no family-wise control. With FDR, the rate is controlled only once B clears `.fdr_min_iterations`, and below that nothing can be flagged.
- **Proposal:** pass `correction` and word the message for each case, or drop the error-rate clause when the correction is not `"fwer"`.

**UM10.** **Proposal:** `stopifnot(!is.null(nm), !anyDuplicated(nm), all(nm %in% colnames(sim)))`. This is a guard for future callers.

## Suggested order

UM1 + UM2 + UM3 together (the p-value logic, all in `.bootstrap_pvalues` / `.wy_stepdown`; baselines on the complete-data paths should be identical), then the messages UM4, UM5, UM9, then UM10, UM7, UM8 docs. UM6 optional.

## Fixes applied: UM1-UM3 (2026-11-11, uncommitted)

- **UM1:** the marginal p-value uses each comparison's non-missing simulated values, `(1 + count) / (1 + B_i)`. The studentising moments already used `na.rm`. The CFA Inf-to-NA conversion (cfa_cutoff.R:688, :720) is unchanged, so an Inf RMSEA is now left out rather than counted as a non-exceedance. Counting it as an exceedance (keeping Inf in the upper tail) is still open for the user.
- **UM2:** `.wy_stepdown(stat_obs, stat_sim)` (the `B` argument is gone) drops iterations with a missing value among the valid comparisons, so every maximum is over the same set and B is the number of complete iterations. No complete iteration gives NA for all padj (the marginal p stays). `max(na.rm = TRUE)` and its warnings are gone.
- **UM3:** the zero-variance guard also checks the SD of the simulated values alone, so a constant simulated column gives NA under `"joint"` as under `"simulated"` and as documented.
- Roxygen of both helpers states the missing-value rule.
- **Checks:** the old version (copy in the session scratchpad) and the new one give identical output on complete data, over 300 random cases x 4 corrections x 2 tails x 2 moments. UM1 probe: .055 before, .109 after, which equals dropping the rows. UM2 probe: 102 warnings before, 0 after. UM3 probe: p = .005 before, NA after.
- **Tests** (test-infit_pvalue.R): UM1 against dropping the rows, per-column denominators, both tails; UM2 step-down equals the complete-iteration run, no warnings, all-incomplete gives NA padj and finite p; UM3 under both moments. The file passes.
- NEWS bug-fix bullet. Full suite not run (user: no full runs without approval).
- Side effect to note: in RMdifGamma/RMlocdepGamma a pair whose gamma is often undefined now removes those iterations from the whole family's step-down. Captions still report the cutoff's iteration count, not the number the step-down used.
- Caller test files run singly (NOT_CRAN=true): test-cfa_cutoff (85 pass), test-ld_partgam_gamma (78), test-partgam_dif (38), test-partgam_ld (67), 0 failures or errors.

## CFA scaled indices (user: all three suggestions, scaled metrics only, 2026-11-11, uncommitted)

Follow-up to the UM1 Inf question. Reading cfa_cutoff.R showed `extract_cfa_fit()` already turned non-finite values into NA, so an Inf never reached the p-values (the comment said otherwise). The real issue was `pick_scaled()`: per dataset it fell back from `cfi.scaled`/`rmsea.scaled` to the uncorrected index, so one null distribution could mix two metrics. lavaan 0.7.3 check (container): WLSMV, ULSMV, WLSM, ULSM, WLSMVS and ULSMVS give `.scaled` indices; DWLS, ULS and WLS give none (they had silently used the uncorrected ones throughout). With WLSMV the scaled indices were never missing in 150 datasets each at n = 100 (targeted, mu = -1.5) and n = 60 (mu = -1.5). User decision: ordinal data only, so scaled metrics consistently.

- `extract_cfa_fit()` returns the scaled CFI/RMSEA only (NA when missing or not finite), and asks lavaan only for `cfi.scaled`, `rmsea.scaled`, `srmr`.
- RMdimCFACutoff refuses estimators without a scaled test ("gives no scaled test statistic ... Use one of: WLSMV, ULSMV, WLSM, ULSM, WLSMVS, ULSMVS"). ML-family message no longer suggests DWLS.
- RMdimCFACutoff messages the count of simulated datasets without a scaled CFI, RMSEA or SRMR ("3 of 6 simulated datasets had no scaled RMSEA. Their cutoffs and p-values use the remaining datasets.").
- Stale Inf comment in RMdimCFA replaced. Roxygen (`estimator`, estimation-model paragraph) updated, Rd regenerated.
- Tests (test-cfa_cutoff.R, NOT_CRAN=true, 94 pass): DWLS/ULS/WLS refused; `extract_cfa_fit()` gives NA for DWLS and the scaled values for WLSMV; the message via a mocked `extract_cfa_fit()`.
- NEWS bullet under Changes.

## Fixes applied: UM4, UM5, UM9 (2026-11-12, uncommitted)

- **UM5:** new `.wy_min_iterations(alpha)`, the first B with `1/(B+1) < alpha` (20 at .05, 100 at .01, 33 at .03). `.warn_fdr_floor()` quotes it; the `.fdr_min_iterations()` roxygen points to it. Checked: one extreme comparison has padj = .05 (not flagged) at B = 19 and .048 at B = 20.
- **UM4:** `.warn_fdr_floor()` multiplies the quoted raw-p threshold by c(m) for BY and gives the formula as c(m)m/(s(B+1)). m = 10: .00171, matching need = 585 (p.adjust confirms 585 reaches it, 584 does not). BH unchanged (.005, 200).
- **UM9:** new `.error_rate_sentence(correction, when)`: fwer "The family-wise error rate is controlled at any number of / this many iterations."; fdr "The false discovery rate is controlled once there are enough iterations for the correction to flag anything."; none "Each p-value is valid on its own, but none is corrected for multiple comparisons." `.notify_low_iterations()` and `.iteration_note()` take `correction` (default "fwer"), and all five callers pass it (RMitemInfit, RMitemRestscore, RMlocdepQ3, RMlocdepGamma, RMdifGamma). RMitemInfit's own caption copy uses the same sentence. RMitemRestscore's caption made no control claim and is unchanged.
- **Tests** (test-locdep-flagging-defaults.R): `.wy_min_iterations` at four alphas plus the 19/20 step-down check; BY and BH warning numbers; notice and caption wording per correction; RMlocdepQ3 with `correction = "none"` end to end. The file passes (101), as does test-infit_pvalue.R (55).
- NEWS bug-fix bullet. Full suite not run.
- Caller test files (NOT_CRAN=true): conditional_infit 56, infit-flagging-defaults 42, item_restscore 18, local_dependence 86, partgam_dif 38, partgam_ld 67, ld_partgam_gamma 78 passed, 0 failures or errors.

## Plain-decimal p-value floor (user request, 2026-11-12, uncommitted)

The captions printed the smallest attainable p-value with `round(1 / (B + 1), 4)` inside `paste()`, which R shows as "7e-04" at B = 1500. New `.fmt_small(x, digits = 2)` (`format(signif(x, digits), scientific = FALSE)`) gives "0.00067". Used in the six kable captions that state the floor (RMitemInfit, RMitemRestscore, RMdifGamma, RMlocdepGamma, RMdimCFA, RMdimResidualPCA) and the two numbers in `.warn_fdr_floor()` (3 digits). Two significant digits instead of four decimals: 1/41 now shows 0.024 rather than 0.0244, 1/401 is 0.0025 either way. Tests in test-infit_pvalue.R (helper values, no e-notation in the BY warning, RMitemInfit caption at a pretended 1500 iterations); it and test-locdep-flagging-defaults.R pass. NEWS folded into the messages bullet.

## Fixes applied: UM10, UM7, UM8 (2026-11-12, uncommitted)

- **UM10:** `.bootstrap_pvalues()` stops with "Internal error in .bootstrap_pvalues(): `observed` needs unique names that are all column names of `sim`" when the names are missing, empty, duplicated or not all columns of `sim`. A named subset in any order still works. No user-facing change (every caller names `observed`).
- **UM7:** `.fmt_width()` uses at least 5 decimals, and `ceiling(-log10(1 - w)) + 1` when that is more, so a width below 1 never prints as "1". Sidak widths: k = 10 ".99488" (unchanged), k = 1225 ".999958" (was ".99996"), k = 19900 ".9999974" (was "1"). w = 1 still prints "1".
- **UM8:** `.p_stars()` code unchanged (strict `<`, matching `padj < alpha`). The roxygen now says it uses iarm's cutpoints but closes the intervals on the left, so exactly .05 gets "." here and "*" in iarm.
- **Tests:** UM10 in test-infit_pvalue.R; UM7 and UM8 in test-locdep-flagging-defaults.R.
- No NEWS entry: UM10 and UM8 change nothing users see, and UM7 only shows at 50+ items. Say if UM7 should get one.
- Test files run singly (NOT_CRAN=true), 0 failures or errors: infit_pvalue 68, locdep-flagging-defaults 109, infit-flagging-defaults 42, conditional_infit 56, item_restscore 18, local_dependence 86, partgam_dif 38, partgam_ld 67, ld_partgam_gamma 78, cfa_cutoff 94, residual_pca 33.

## Fix applied: UM6 (2026-11-12, uncommitted)

`.wy_stepdown()` computes the successive maxima as a reverse running maximum (`pmax()` from the last-ranked comparison back to the first, then `cummax()` for monotonicity): O(K B) instead of O(K^2 B).
- Identical output to the pre-UM6 version (copy in the session scratchpad) over 1600 random `.bootstrap_pvalues(correction = "fwer")` cases with ties, tied observed values, missing values and constant columns, both tails and both moments.
- k = 1225 pairs (50-item Q3), B = 1000: 6.21 s before, 0.21 s after (container), identical results.
- Test (test-infit_pvalue.R) against a direct reference implementation of the step-down, 30 random cases with ties. test-infit_pvalue (98), test-locdep-flagging-defaults (109) and test-local_dependence (86) pass.
- No NEWS entry (results unchanged). Could add "faster for large item sets" if wanted.

**Unit status:** utils-multiplicity.R done (UM1-UM10, plus the CFA scaled-index and plain-decimal follow-ups). Full suite not run in this unit.

Full suite (installed tree, 2026-11-12, after UM1-UM10, CFA scaled indices, plain-decimal floor, median SEM, dif_lr tests): 0 failures, 0 errors, 1997 pass, 7 skips, 1 warning (the known mirt EM warning in test-reliability.R, RMreliability boot sequential).

## Second pass: `/code-review high` on the diff (2026-11-12, not verified)

Run by the user in a fresh session, after all of the above. It reviewed the uncommitted working tree at 1.4.0.9003 (no git repo here, so the "diff" is the set of changes in this unit): `utils-multiplicity.R`, the scaled-only CFA indices, the median SEM in RMreliabilityCurve, the iteration notices and captions, `.fmt_small()` and `.fmt_width()`. Ten findings, numbered on from the first pass. **None has been verified.** The reviewer ran no code, and neither has the main session yet. Line numbers refer to the working tree on 2026-11-12. Nothing fixed.

| # | File:line | Finding | Status | Reaches users? |
|---|---|---|---|---|
| UM11 | utils-multiplicity.R:646 | The step-down keeps only complete iterations, so one comparison that is often NA shrinks B for the whole family, silently. Captions and `.notify_low_iterations()` still use the cutoff's iteration count | reported | yes: RMlocdepGamma, RMdifGamma (sparse pairs or items) |
| UM12 | utils-multiplicity.R:646 | The marginal p uses each column's non-NA rows (line 96), the step-down only complete rows. padj can then fall below p, and conditioning on complete iterations may not be missing at random | reported | yes, same paths as UM11 |
| UM13 | cfa_cutoff.R:686 | Without the fallback, an observed scaled CFI/RMSEA that lavaan cannot compute is NA, so `compute_cfa_flagged()` gives not flagged and p = NA with no message. Only simulated NAs are reported | reported | yes: RMdimCFA |
| UM14 | cfa_cutoff.R:686 | RMdimCFA and RMdimCFAPlot refit with `cutoff$estimator` without checking it. A saved 1.4.0 cutoff object built with DWLS/ULS gives NA observed indices against a null of uncorrected values | reported | yes: cached cutoff objects from 1.4.0 |
| UM15 | cfa_cutoff.R:890 | The RMdimCFA caption floor uses `actual_iterations`, but each index's p uses its own non-NA count and the loadings step-down only complete iterations. RMdimCFACutoff reports no count for missing loadings | reported | yes, wrong number in a caption |
| UM16 | utils-multiplicity.R:401 | The once-per-session `.frequency_id` of `.notify_low_iterations()` does not include `correction`, so after one notice, a call with another correction shows none | reported | yes, missing notice (follows from UM9) |
| UM17 | conditional_infit.R:583 | RMitemInfit builds its two iteration tiers inline instead of calling `.iteration_note(n_iter, correction)`, and the two copies already differ (the extra "Two analysts ... 10% ... 4%" sentence) | reported | no, maintenance |
| UM18 | utils-multiplicity.R:72 | `s_sim` is always a second column-wise `sd` pass. Under `moments = "simulated"` it equals `s` | reported | speed only |
| UM19 | utils-theta.R:431 | `.marginal_summaries()` always sorts the nodes and interpolates the median SEM, which RMreliability (reliability.R:644 and its bootstrap workers) discards | reported | speed only, per resample |
| UM20 | utils-multiplicity.R:96 | The marginal p drops NA but not Inf. A caller that skips the non-finite-to-NA step gets s = NaN and an NA p-value | reported | no, current callers convert first (RMdimCFA does) |

**UM11.** Already noted as a side effect under UM2 above. The reviewer's examples: 20 items (190 pairs) with one pair undefined in 95% of 1000 iterations leaves B of about 50 and a padj floor of 1/51. If each of the 190 pairs is undefined 1% of the time, independently, about 15% of iterations survive (0.99^190 = .15). Below B = 20 nothing can be flagged at .05 (UM5), and the caption still quotes 1/1001.
- **Options:** (a) keep complete-case rows, report the B the step-down used in the caption, and message when it is well below the cutoff's count, or (b) handle NAs per comparison in the step-down. Either changes behaviour, so it needs a decision.

**UM12.** Check first whether `.wy_stepdown()` or its caller already enforces `padj >= p`. If not, a probe with NAs placed in the most extreme rows of one column would show whether padj < p can occur. The missing-at-random question needs a simulation where undefined gammas are tied to sparse resamples, so ask before running it.

**UM13.** **Proposal:** message when an observed scaled index is NA, naming the index, and show it as not tested rather than not flagged.

**UM14.** **Proposal:** in RMdimCFA and RMdimCFAPlot, run the same estimator check as RMdimCFACutoff on `cutoff$estimator`, so an old object stops with the new error.

**UM15.** **Proposal:** quote the floor per index from its non-NA count, or give the smallest of them, and have RMdimCFACutoff count missing loadings as it does the fit indices. Resolve together with UM11.

**UM16.** **Proposal:** add `correction` to the `.frequency_id`.

**UM17.** **Proposal:** call `.iteration_note()` from RMitemInfit and keep the extra sentence as an argument or a separate line. Check the wording against the RMitemInfit tests.

**UM18.** **Proposal:** compute `s_sim` only when `moments == "joint"` and reuse `s` otherwise.

**UM19.** **Proposal:** compute `sem_median` only when asked (an argument, or in RMreliabilityCurve itself).

**UM20.** **Proposal:** filter on `is.finite()` in the marginal p and the moments, or say in the roxygen that callers must pass NA, not Inf. The upper-tail Inf question from UM1 still applies.

## Suggested order (second pass)

Verify first: UM12 (the padj < p check), UM13 and UM14 with small probes, UM11 by counting complete iterations in a sparse RMlocdepGamma example. Then the UM11 + UM12 + UM15 decision together, since they share the complete-iteration rule. Then UM13, UM14, UM16, UM20. UM17-UM19 are optional cleanups.

## Verification of UM11-UM20 (main session, 2026-11-12)

Probe scripts in the session scratchpad (`cr/probe_um11.R`, `cr/probe_um11b.R`, `cr/probe_um12.R`, `cr/probe_um12b.R`). Nothing fixed yet.

| # | Verdict | Evidence |
|---|---|---|
| UM11 | confirmed in the code, rare in practice | RMlocdepGammaCutoff, 200 iterations: n = 60, 8 items, 4 categories, mean -1.5 gave 1 of 28 pairs NA in 1 of 93 iterations (92 complete). n = 100, 10 items: no NA. A sparser design stopped on an empty category before any gamma could be NA. RMdifGammaCutoff with a 10-20% small group (n 80-200, 2 or 4 categories): no NA in 200 iterations. Sparse data mostly fail the empty-category rule first, so the complete-iteration loss stays small |
| UM12 | confirmed | 2000 random cases with 0-20% missing values per cell (B = 200, 2-8 comparisons): padj < p in 518. Mostly tiny (median ratio 0.99), the worst 0.45, and in 1 case p >= .05 with padj < .05. My targeted probe (NA in one column where another is extreme) did not show it. It comes from the two p-values using different subsets of iterations |
| UM13 | confirmed by reading | `compute_cfa_flagged()` treats NA as not flagged, and the p is NA with no message. How often lavaan returns NA for an observed scaled index under WLSMV is unknown. It was never NA in the earlier simulated checks down to n = 60 |
| UM14 | confirmed by reading | RMdimCFA (cfa_cutoff.R:680) and RMdimCFAPlot (line 1116) pass `cutoff$estimator` to lavaan unchecked |
| UM15 | confirmed by reading | Caption floor at line 893 uses `actual_iterations`. RMdimCFACutoff counts missing fit indices (line 347) but not missing loadings. Same root as UM11, same rarity |
| UM16 | confirmed | `.frequency_id = id`, with no `correction` in it |
| UM17 | duplication confirmed, difference intentional | The `.iteration_note()` roxygen says the item seed-disagreement percentages are deliberately left out for pairs. A refactor would need an `extra` argument. Maintenance only |
| UM18 | confirmed | One extra `sd` pass per column. Negligible next to the simulations |
| UM19 | confirmed | One sort over the quadrature nodes per call. Negligible next to a mirt fit per bootstrap resample |
| UM20 | confirmed by reading | All current callers convert non-finite values to NA first |

## CFA minimum raised to 4 items (user decision, 2026-11-12, uncommitted)

The user ran a 3-item model. The fit part was empty while loadings worked. A one-factor model on 3 ordinal items is just-identified (df = 0). lavaan WLSMV gives CFI 1, RMSEA 0 and SRMR 0 for the observed data and for every simulated dataset, so the null is a single point and the UM3 guard sets the p-values to NA. The user chose to require 4 items, since the fit indices are the key part of the function.
- `RMdimCFACutoff()` stops below 4 items, with a df = 0 explanation in the message. `@param data` and NEWS were updated.
- A cached 3-item cutoff object would still run in RMdimCFA/RMdimCFAPlot. Add an item-count check there together with the UM14 estimator check.
- test-cfa_cutoff.R: 95 pass, 1 skip.

## Fixes applied: UM11, UM12, UM13, UM14, UM16 (2026-11-12, uncommitted)

- **UM12.** `.bootstrap_pvalues()` sets `padj <- pmax(padj, p)` under `"fwer"`. It changes results only when simulated values are missing.
- **UM11 (and the iteration part of UM15).** `.wy_stepdown()` returns the number of complete iterations as attribute `B`, which is stripped before return. `.bootstrap_pvalues()` messages when it is below 95% of the iterations, with the step-down's own floor 1/(B+1). This covers the CFA fit and loadings families too, since they go through the same helper. The caption still quotes the cutoff's count, and the message covers the rare case.
- **UM13.** RMdimCFA: an observed fit index that is NA gets `Flagged = "not tested"` and a message. `@return` updated.
- **UM14.** New `.check_cfa_cutoff_object()`, called by RMdimCFA and RMdimCFAPlot, stops on a non-scaled estimator or fewer than 4 items. `.cfa_scaled_estimators` is shared with RMdimCFACutoff.
- **UM16.** `.frequency_id = paste0(id, "_", correction)`. The rlang resets in six test files were updated to the per-correction ids.
- Not done, as recommended: UM17 (difference intentional), UM18, UM19 (negligible), UM20 (callers convert first).
- NEWS: four bug-fix bullets.
- Tests: test-infit_pvalue (506), test-cfa_cutoff (102, 1 skip), test-locdep-flagging-defaults (109), test-infit-flagging-defaults (42), test-dif-flagging-defaults (48), test-partgam_dif (38, 1 skip), test-partgam_ld (67, 1 skip), test-local_dependence (86): 0 failures. Full suite not run.
- Full suite (installed tree, 2026-11-12, after this batch, the 4-item minimum and the condensed NEWS): 0 failures, 0 errors, 2413 pass, 7 skips, 2 warnings. One was the known mirt EM warning in test-reliability.R. The other was the expected 20-iteration warning in the new UM13 test, which now expects it (test-cfa_cutoff rerun: 103 pass, 0 warnings).
