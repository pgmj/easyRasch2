# Code review: equal-tailed restscore p-values (uncommitted 1.4.0.9004 diff)

Reviewed 2026-11-13 with `/code-review high` on the uncommitted diff. The diff adds `tail = "split"` to `.bootstrap_pvalues()` (new `.split_pvalues()` in `R/utils-multiplicity.R`), wires it into `RMitemRestscore()`, and updates the tests in `test-item_restscore_cutoff.R` and `test-infit_pvalue.R`. Background is in `dev/HANDOFF-restscore-equal-tails.md`.

The findings below come straight from the review agent. **None have been verified in R yet.** Verify each one before fixing it, as for UM1-UM12 in `utils-multiplicity.md`. Line numbers refer to the working tree on 2026-11-13.

Status key as in utils-theta.md. Verified and fixed 2026-11-13, see "Verification and fixes" at the end.

## Findings

| # | Where | Finding | Status | Reaches users? |
|---|---|---|---|---|
| ET1 | `RMitemRestscore()`, item_restscore.R:~470 | No `.warn_fdr_floor()` call. Under fdr_bh/fdr_by, a run with too few iterations to flag anything stays silent. The split tail makes this worse because it doubles the floor | confirmed, fixed | yes, silent no-power runs |
| ET2 | `.fdr_min_iterations()`, `.warn_fdr_floor()`, `.wy_min_iterations()`, utils-multiplicity.R:~600-640 | All three hard-code a `1/(B+1)` floor. Under `tail = "split"` the floor is `2/(B+1)` | confirmed, fixed | only once ET1 is wired in, but it blocks ET1 |
| ET3 | test-item_restscore_cutoff.R:183 | The planted-underfit test uses `iterations = 40`, so the floor is 2/41 = .0488, only just under .05 | plausible, fixed | no, test fragility |
| ET4 | test-infit_pvalue.R:447 | `expect_equal(rate[1], 2 / (B + 1))` depends on the seed | confirmed, fixed | no, test fragility |
| ET5 | `.split_pvalues()` `one_sided()`, utils-multiplicity.R:164 | Duplicates `.mc_pvalue(tail = "two.sided")` (person_fit.R:677). The two disagree when every simulated value is NA | no issue (NA case unreachable); ET9 part fixed | edge only |
| ET6 | `.split_pvalues()`, utils-multiplicity.R:~183 | The Westfall-Young attrition `message()` is copied from `.bootstrap_pvalues()`, differing only in 1 vs 2 | plausible, fixed | no, maintenance |
| ET7 | vignettes/easyRasch2.Rmd:244 | The precomputed caption still reads "Two-sided ... 1/(100+1) = 0.0099". The function now prints "Equal-tailed ... 2/(100+1) = 0.0198" | known | yes, until re-knit |
| ET8 | `.bootstrap_pvalues()` comment, utils-multiplicity.R:84 | "the test statistic is \|t\| (two-sided) or t (upper tail)" does not mention `"split"`, which returns on the next lines | plausible, fixed | no |
| ET9 | `.split_pvalues()`, utils-multiplicity.R:172-173 | Two `vapply` passes, so each column is NA-filtered twice | plausible, fixed | speed only, negligible |

**ET1.** The agent's example: 20 items, `iterations = 400`, `correction = "fdr_bh"`. The smallest BH-adjusted p is 20 × 2/401 = .0998, so nothing can reach `padj < .05`. The caption only quotes the raw floor, 2/(400+1) = .005, which suggests flagging is possible. The gap existed under `|t|` (floor 1/(B+1)), but the split tail roughly doubles the iterations needed (about 800 → 1600 for 20 items).
- **How others do it:** `dif_partgam.R:288`, `local_dependence.R:276` and `ld_partgam.R:284` already call `.warn_fdr_floor()`. Copy their pattern.
- **Depends on ET2.** Wiring in the current helper unchanged would quote the wrong floor and half the iterations needed.

**ET2.** Proposal: add a `floor_mult = 1` argument (or `tail = "two.sided"`) to all three helpers, and pass 2 under `"split"`.
- **`.fdr_min_iterations()`:** use `reachable <- function(b) (penalty * m) * (floor_mult / (b + 1)) < alpha`, and adjust the starting value of `b` to match.
- **`.wy_min_iterations()`:** it gives 20 at alpha = .05. Under split, `2/(B+1) < .05` needs B = 40.
- **`.warn_fdr_floor()`:** the message text and `need` must use the multiplier.
- **Check:** keep the existing `p.adjust()` bit-agreement comment true, and test at `m / alpha` whole numbers.
- **Unchanged callers:** Q3, LD gamma and DIF gamma use the default `tail`, so their results must not change.

**ET3.** If one simulated dataset fails to refit (`actual_iterations = 39`), the floor becomes 2/40 = .05, which is never `< .05`, so I6 cannot be flagged and `expect_identical(..., "underfit")` fails. The same happens if any single simulated row has a min t below I6's observed t. Under `|t|` the floor was 1/41, so the test was robust. Proposal: raise to `iterations = 100` and check the runtime stays acceptable. Run only this test file afterwards.

**ET4.** Rows are exchangeable, so the row holding the global max of the studentised values gets `padj_hi = 1/(B+1)`, and the row holding the global min gets `padj_lo = 1/(B+1)`. Each of those rows gets padj = 2/(B+1).
- **Usual case:** the max and min are in different rows, so two rows reach 2/(B+1) and `rate[1] = 2/(B+1)`.
- **When the test fails:** if one row holds both the max and the min (roughly 1 in 60 seeds at B = 59), only one row reaches it, and `rate[1] = 1/(B+1)`.
- **Current state:** the test passes now with its fixed seed (verify this), but a change in RNG or seed can break it without anything being wrong.
- **Proposal:** keep `expect_true(all(rate <= lv + 1e-12))`. Replace the exact equality with a check that the bound is reached across seeds (e.g. `rate[1] >= 1/(B+1)`, or average over a few seeds), or drop it.

**ET5.** Proposal: `p_marg` from `.mc_pvalue(t_sim[, i], t_obs[i])`, keeping per-tail counts only for `direction`.
- **NA mismatch:** when every simulated value is NA, `.mc_pvalue` returns NA but `one_sided()` returns 1. NA is the right answer and matches UM1's NA policy. Decide this before merging.
- **Combine with ET9:** a single helper returning both tail counts fixes both.

**ET6.** Proposal: a small `.wy_attrition_message(b_used, b_total, floor_mult)` used by both functions. This only touches wording, so the results do not change.

**ET7.** Re-knit the vignette before release (already noted in the handoff and in the restscore memory note).

**ET8.** Change the comment to "|t| (two-sided), t (upper tail), or both tails separately (split, see `.split_pvalues()`)".

**ET9.** Fold into ET5.

## Suggested order

1. ET2, then ET1 (same edit session).
2. ET3, ET4 (tests).
3. ET5 + ET9, ET6, ET8 (cleanups, optional).
4. ET7 at release.

Test scope: `test-item_restscore*.R`, `test-infit_pvalue.R`, plus the test files for any helper changed in ET2. Those helpers are shared, so also run the DIF gamma, LD gamma and Q3 tests to confirm their output is unchanged. Do not commit; the user commits.

## Verification and fixes (2026-11-13)

Checked in R first (script in the session scratchpad, not kept):

- **ET1 confirmed.** 20 items, B = 400, `fdr_bh`, one item at t = 10: the smallest BH padj was .0998 = 20 x 2/401, with no warning. Three standouts reach .033, so the warning wording ("several extreme items lower that requirement") stays right.
- **ET4 confirmed.** Over 200 seeds, `rate[1]` was 1/(B+1) in 3 of them. The level bound `rate <= lv` held in all 200.
- **ET5: the NA mismatch can't happen.** A column with fewer than two non-NA simulated values has a missing SD, so the zero-variance guard sets `t_obs` to NA before `one_sided()` runs, and p is NA. Checked with an all-NA column and with a single non-NA column. `.mc_pvalue()` is not reused because `direction` needs the per-tail counts anyway.

| # | Fix |
|---|-----|
| ET2 | `floor_mult = 1L` argument on `.fdr_min_iterations()`, `.wy_min_iterations()` and `.warn_fdr_floor()`. The default leaves Q3, LD gamma and DIF gamma unchanged. `p.adjust()` bit agreement holds because doubling is exact in binary; tested at m = 9, 20, 21, 36, 190 and alpha = .05, .01 |
| ET1 | `RMitemRestscore()` calls `.warn_fdr_floor(..., floor_mult = 2L)` after the low-iteration notice, as DIF gamma does. Test at 20 iterations: warns under `fdr_bh` and quotes 2/(20+1), silent under `fwer` |
| ET3 | Planted-underfit test raised to 100 iterations |
| ET4 | Exact equality replaced by `rate[1] >= 1/(B+1)`, with the reason in a comment |
| ET6 | `.wy_attrition_message(b_used, b_total, floor_mult)`, used by both functions |
| ET8 | Comment now names the split statistic |
| ET9 | One `vapply` returns both tail p-values |
| ET7 | Open, re-knit at release |

New tests: the floor helpers with `floor_mult = 2` (40 WY iterations at .05, BH and BY boundaries against `p.adjust()`, warning text), and two standouts in opposite directions each getting padj = 2/(B+1).

Checks after the fixes:
- **Tests:** these files pass with NOT_CRAN = true: `test-infit_pvalue.R`, `test-item_restscore_cutoff.R`, `test-item_restscore.R`, `test-bootstrap_restscore.R`, `test-locdep-flagging-defaults.R`, `test-partgam_dif.R`, `test-partgam_ld.R`, `test-ld_partgam_gamma.R`.
- **Reproduction:** the package still matches `split_t_pooled` on all 2000 datasets in `dev/restscore_split_null/` (0 disagreements; FWER 5.00, 2.40 under, 2.75 over).

Not in scope, noted: `RMitemInfit()` (conditional_infit.R) also has no `.warn_fdr_floor()` call. That gap predates this change, and its floor is 1/(B+1).
