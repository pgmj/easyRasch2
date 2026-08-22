# Plan: local dependence flagging defaults for 1.2.0

Status: proposed, nothing implemented. Written 2026-08-22.

Moves `RMlocdepQ3()` and `RMlocdepGamma()` from interval flagging to corrected
bootstrap p-values, completing the change `RMitemInfit()` made in 1.1.1.9000.
Everything below is one release.

## Why the change is defensible now

The Q3 and partial gamma simulation study (`rasch_q3fwer/`) is not finished, but
it does not have to be. The study asks which statistic and which familywise
procedure detect local dependence best. The default asks whether the shipped
rule controls the familywise error rate at all. Only the second question bears
on this change, and it rests on three things that are already settled.

**Arithmetic.** A per-pair interval is a per-comparison rule, so its width sets
a familywise error rate of `1 - width^m`. Pairs grow quadratically in items, so
`hdci_width = 0.99` implies about 30% over the 36 pairs of a 9-item scale, 65%
over the 105 pairs of 15 items, and 85% over the 190 pairs of 20 items. The
same argument was accepted for `RMitemInfit()`, where the family grows only
linearly.

**The source literature.** Christensen et al. (2017) propose maximum-statistic
rules, which are familywise by construction, and offer per-pair critical values
only for a pair specified a priori. Debelak et al. (2020) advise nominating
pairs in advance because a Bonferroni-type correction would cost too much
power, which argues for a better correction rather than none (Ferreira, 2024).
Kreiner and Christensen (2011) run item screening with partial gamma corrected
by Benjamini-Hochberg at 5%, and invite substitution of a better multiplicity
adjustment. Every source the two functions rest on controls multiplicity. The
current default is the only link in the chain that does not.

**The spurious negative halo.** The interval is two-sided, so it also flags
pairs below the lower bound, which reads as multidimensionality. Kreiner and
Christensen (2011, Example 3) show one genuinely dependent pair driving every
pair that shares an item with it significantly negative, in partial gamma. The
`rasch_q3fwer` run found the same suppression (upper-tail false positive rates
.002 to .004 for sharing pairs against .078 to .116 for disjoint pairs). The
below-bound flags are mostly that halo, so losing them is the point of the
change rather than a cost of it.

## Measured, 2026-08-22

The study's rule set contains no rule that is the package's interval, so the
shipped defaults were measured directly: phq9 CML parameters, 9 items, 36
pairs, n = 400, Gaussian copula, each function at its own default iteration
count, both rules read from the same cutoff object per replication.
Scripts in the session scratchpad, results reproducible from
`rasch_q3fwer/R/helpers.R`.

Complete null, at least one pair flagged:

| Function | interval | Westfall-Young | replications |
|---|---|---|---|
| `RMlocdepQ3()`, B = 500 | .344 (SE .030) | .060 (SE .015) | 250 |
| `RMlocdepGamma()`, B = 250 | .410 (SE .022) | .042 (SE .009) | 500 |

Against `1 - .99^36 = .304` predicted. Both intervals exceed the prediction,
and gamma exceeds it further, which is the .99 interval being estimated from
250 draws.

One weak dependent pair (q1+q2, rho calibrated to Q3 = .08):

| | Q3 interval | Q3 WY | gamma interval | gamma WY |
|---|---|---|---|---|
| detects the dependent pair | .928 | .832 | .898 | .754 |
| >=1 false flag, 21 disjoint pairs | .272 | .032 | .264 | .042 |
| >=1 flag, 14 pairs sharing an item | .260 | .000 | .338 | .006 |

The interval buys 10 points of detection for Q3 and 14 for gamma, and pays for
them with a false flag somewhere in the matrix in roughly a quarter of
analyses. The Q3 Westfall-Young figure of .832 matches the `rasch_q3fwer` v1
run at the same cell (.834), so the study's numbers transfer to the package.

Scope: one sample size, one item set, one mechanism, one magnitude. The
familywise half is the robust half, being arithmetic that the measurement
confirms in both functions. The v1 run carries the four sample sizes and three
mechanisms, for Q3 only.

## Changes

### 1. `p_value = NULL` in both LD functions

- `R/local_dependence.R:158`, `R/ld_partgam.R:156`: `p_value = FALSE` becomes
  `p_value = NULL`.
- Resolution mirrors `R/conditional_infit.R:295-320`. `NULL` means TRUE when
  the full cutoff object carrying the simulated distributions was supplied, and
  FALSE otherwise. The field is `cutoff_full$pair_results` for Q3 and
  `cutoff_full$results` for gamma.
- A bare numeric cutoff (Q3) or a bare `$pair_cutoffs` data.frame (gamma)
  resolves to FALSE, so those calls are unchanged.
- No cutoff at all resolves to FALSE. `RMlocdepGamma(data)` keeps the
  asymptotic `iarm` path exactly as it is today. This change cannot reach it.
- Explicit `p_value = TRUE` without simulations stays an error
  (`R/local_dependence.R:216`, `R/ld_partgam.R:231`).

**Behaviour delta.** Scripts that pass the full cutoff object and do not name
`p_value` get two new columns and a `Flagged`/`flagged` column that may name
different pairs. Fewer pairs will be flagged, and no pair will be flagged for
sitting below the lower bound. `p_value = FALSE` restores the old behaviour.

### 2. `hdci_width` 0.99 to 0.95 in both cutoff functions

- `R/local_dependence.R:846`, `R/ld_partgam.R:772`.
- The interval stops being a decision rule and becomes a description of where a
  pair's statistic is expected to fall, which is what `RMitemInfitCutoff()`
  already did in moving from .999 to .95.
- The width and the iteration count have to agree. `.width_iterations()` puts a
  .99 interval at 2000 draws and a .95 interval at 400. Both LD cutoff
  functions currently ship a .99 interval estimated from 250 or 500 draws,
  which is why the measured familywise rates (.344 and .410) exceed the .304
  the width implies. Moving to .95 makes 400 the matched count, which is item 3.
- Affects `$pair_cutoffs`, the plots that overlay them, and any script reading
  `gamma_low`/`gamma_high` or `Q3_low`/`Q3_high`.

### 3. Both cutoff functions default to 400 iterations

- `R/ld_partgam.R:766` (250 to 400) and `R/local_dependence.R:840`
  (500 to 400). Decided 2026-08-22 in favour of one number across the package.
- 400 is the floor measured by Johansson (2026) for the Westfall-Young
  correction, the level of the `rasch_q3fwer` B = 400 arm, which is
  indistinguishable from B = 1000 for Westfall-Young at 36 pairs, the
  `RMitemInfitCutoff()` default since 1.1.1, and the count `.width_iterations()`
  matches to the new .95 width. The measurement above shows the correction
  already holding at 250 (.042), so this is alignment rather than a fix.
- Lowering Q3 from 500 to 400 also narrows the interval slightly, which matters
  less once the interval is descriptive, and the caption still recommends
  1000 to 2000 for a final analysis.

### 4. `requested_iterations` in both cutoff objects

- Add beside `actual_iterations`, following `R/infit_cutoff.R:317`.
- `.notify_low_iterations()` needs it to distinguish a low request from
  iterations lost to failed fits.

### 5. Replace the two hard warnings with the two-tier notice

- `R/local_dependence.R:224-237` and `R/ld_partgam.R:241-251` warn whenever
  p-values are computed from fewer than 1000 iterations. Both cutoff functions
  default below 1000, so flipping the default without touching these makes
  every default call warn. This is the pattern the infit work removed, and it
  would land on two more functions at once.
- Use `.notify_low_iterations()` below 400 (once per session) and a caption
  sentence between 400 and 1000, as `RMitemInfit()` does.

### 6. Generalise `.notify_band_flagging()` and `.band_error_clause()`

- `R/utils-multiplicity.R:180` and `:207`. Both compute `1 - width^k` correctly
  for any count, but their wording says "items" and the notice names
  `RMitemInfitCutoff()`.
- Add a noun and a function name, or take a pre-built label. Keep the Sidak
  helpers untouched.
- Then call the notice from the `p_value = FALSE` branch of both LD functions
  with the pair count, so a user who opts back into the interval is told what
  it costs over 36 or 190 comparisons.
- The interval branch is otherwise unchanged and keeps flagging pairs below the
  lower bound (decided 2026-08-22). The p-value branch stays one-sided upper,
  so the two branches deliberately test different things and the captions have
  to say which is in force.

### 7. Captions

- Both p-value branches already append the correction label and the alpha. Add
  the reproducibility sentence between 400 and 1000 iterations, and the
  familywise clause to the interval branch, matching `RMitemInfit()`.
- Say that the p-value test is one-sided upper while `Low`/`High` remain
  two-sided description, so the columns and the flag are not read as the same
  test.

### 8. Reachability warning for `fdr_bh` and `fdr_by`

Included on the user's decision, 2026-08-22. Fires when `correction` is
`"fdr_bh"` or `"fdr_by"` and the cutoff object carries too few iterations for
the procedure to reject anything.

**The arithmetic, corrected.** A Monte Carlo p-value cannot fall below
`1/(B+1)` (`R/utils-multiplicity.R:54`). When `s` of the `m` p-values sit at
that floor, the smallest attainable adjusted value is

```
m / (s * (B + 1))            BH
c(m) * m / (s * (B + 1))     BY,  c(m) = sum(1/i, i = 1..m)
```

verified against `stats::p.adjust()` for m in {36, 190}, s in {1, 2, 3, 5} and
B in {250, 400, 1000}, exact in every cell. The familiar `B >= m/alpha - 1` is
the `s = 1` case, and the earlier note in this project stated it as though it
were general. It is not. The requirement is

```
B >= m / (s * alpha) - 1
```

so it relaxes as more pairs reach the floor. At m = 36 and alpha = .05 the
required B is 719 for s = 1, 359 for s = 2, 239 for s = 3 and 179 for s = 4.
The same correction was already made in the infit paper, where k = 1 detection
jumps 0, 8 and 98 percent across B = 100, 200 and 400 exactly at the floor,
while k = 3 already reaches 71 percent at B = 200.

**What the warning should therefore say.** `s = 1` is the worst case and the
one to warn on, since a single genuinely dependent pair is the situation the
user cannot rule out in advance. The wording must not claim the procedure
"cannot flag anything", because with several extreme pairs it can.

Required B at alpha = .05, single standout:

| items | pairs | BH | BY |
|---|---|---|---|
| 9 | 36 | 719 | 3 005 |
| 15 | 105 | 2 099 | 10 995 |
| 20 | 190 | 3 799 | 22 142 |

Westfall-Young needs `B >= 1/alpha - 1`, which is 19 and does not grow with the
number of comparisons. That contrast is the reason `"fwer"` is the default and
belongs in the warning.

**Implementation notes.**

- Compute from the function's own `alpha` argument, not a hardcoded .05.
- Key on `actual_iterations`, not the requested count. `.bootstrap_pvalues()`
  takes `B <- nrow(sim)`, which is the number of successful iterations carried
  in `$pair_results` or `$results`.
- `c(m)` is `sum(1/seq_len(m))`. c(36) = 4.175, c(190) = 5.827.
- Best home is a helper in `R/utils-multiplicity.R` beside the Sidak ones, so
  `RMitemInfit()` and the DIF functions can use it later. Their family is items
  or item-by-group cells rather than pairs, so the counts differ but the
  formula does not.

## Tests

`tests/testthat/test-infit-flagging-defaults.R` is the template. A sibling file
covering both LD functions should assert:

- the resolution of `p_value = NULL` in all four cases (full object, bare
  cutoff, no cutoff, explicit TRUE without simulations)
- the new `hdci_width` and iteration defaults
- the notice below 400 iterations and silence at or above it
- the familywise notice on the interval branch, with the pair count and not the
  item count
- both caption branches
- `requested_iterations` recorded
- the reachability warning fires for `fdr_bh` and `fdr_by` below the threshold
  and is silent above it, with the threshold computed from `alpha`
- the closed form `m/(s(B+1))` against `stats::p.adjust()`, so the corrected
  arithmetic is pinned rather than restated in a comment

Existing tests that call with a full cutoff object and no `p_value` change
meaning and must pin `p_value = FALSE` to keep testing the interval:
`test-local_dependence.R:262-277` and `test-partgam_ld.R:160-268`, in
particular the gamma test at `:225-232` that asserts flagging on `gamma_pair`.

Shared helpers change, so the full suite is warranted rather than the two files.

## Documentation

- `NEWS.md`: one block under Breaking changes covering items 1 to 3, with the
  meaning change to `Flagged`/`flagged` in bold and `p_value = FALSE` named as
  the way back. Items 4 to 8 get one bullet each under Other changes.
- `README.md:70` and `:82` describe `p_value = TRUE` as opt-in. Reword.
- `?easyRasch2-reproducibility` is unaffected.
- `_pkgdown.yml` is unaffected, no exported functions change.
- Roxygen for the four functions: `@param p_value`, `@param hdci_width`,
  `@param iterations`, and the `@details` paragraphs that describe flagging.

## Out of scope

- The partial gamma asymptotic path and its family of 72
  (`ld_partgam.R:277`, `:513`, `:521`, `dif_partgam.R:258`). The mislabelled
  "BH" column, the one-sided recomputation and the pairwise family are a
  separate decision that the study bears on directly.
- `RMdifGamma()`, `RMdimCFA()` and `RMdimResidualPCA()`, which keep
  `p_value = FALSE`. The rationale for moving two of the six is that infit, Q3
  and partial gamma LD are what the two studies examined.
- The argument-order difference between `RMlocdepQ3()` and `RMlocdepGamma()`.
- Extending the item 8 warning to `RMitemInfit()` and the DIF functions, which
  share the correction argument and would use the same helper.

## Decisions taken 2026-08-22

1. The interval branch is unchanged and keeps flagging below the lower bound.
   The p-value branch stays one-sided upper.
2. `RMlocdepQ3Cutoff()` moves from 500 to 400, so all three cutoff functions
   share one iteration default.
3. The reachability warning ships in this release, on the corrected arithmetic
   in item 8.

Nothing is left open. The plan is ready to implement on request.

## Verification before finishing

- `devtools::document()` with zero warnings, then `test-documentation.R`.
- Full `devtools::test()`.
- Re-run the two measurement scripts against the changed defaults and confirm
  the familywise rate lands near .05 for both functions with no arguments
  beyond the data and the cutoff object.
- `R CMD check --as-cran` under `LC_ALL=en_US.UTF-8`.
