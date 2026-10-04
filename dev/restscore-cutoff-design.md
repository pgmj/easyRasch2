# Design: RMitemRestscoreCutoff()

Status: BUILT 2026-09-30 (function, consumer path, plot, tests, docs). Decisions taken: dgp default "resample", `cutoff` second, sample-size check added to all cutoff consumers, plot now. Decision 3 done (doc note + caption on the asymptotic path). Decision 6 dropped: the Boot sentence is about practical vs statistical significance, not calibration, so it stands. Validation run done 2026-10-01 (`restscore_cutoff_validation.qmd`).

**Revised 2026-10-01 after the validation run.** (1) The `dgp` default is now
`"conditional"`: pooled WY family-wise rate 5.0 percent against 6.35 for
`"resample"` on the same datasets, at 1.2 to 1.6 times the run time. This
breaks the signature match with `RMitemInfitCutoff()` on purpose. (2) The
p-values stay two-sided on |t|. An equal-tailed version was built and then
withdrawn on 2026-10-02: the rank-based Westfall-Young step-down it needed
cannot reject with 20 items at B = 400 (about 2k simulated rows tie at the
largest |z|), and continuous alternatives either held the level but left the
family-wise flags lopsided or reintroduced the imbalance. See
`restscore_eti_null.qmd`. The tail imbalance of |t| (underfit 2.6 to 3.4
percent per item, overfit 1.8 to 2.7) is documented in the help page.

**Revised 2026-10-03: an equal-tailed test that works, not yet shipped.**
Two one-sided Westfall-Young step-downs, on t for overfit and on −t for
underfit, each at α/2, with t studentised by the mean and SD of all B + 1
values (`split_t_pooled`). Validated on synthetic exchangeable nulls, the
real restscore null, power (300 datasets per cell) and seed stability:

- Family-wise 5.0 percent on the real null, split 2.40 underfit / 2.75
  overfit, against 3.70 / 1.55 for |t|.
- 2.9 and 2.2 points less underfit detection than |t| (SE 0.4), no
  measurable overfit gain, fewer false underfit flags when items overfit
  (4.7 against 6.2 percent with two overfitting items).
- As seed-stable as |t| at B = 400 (10.3 against 10.0 percent of seed pairs
  disagreeing), more stable at B = 1000 (4.6 against 6.1).

Two findings bear on the shipped design. **Self-influence:** studentising by
the moments of the simulated values only is slightly liberal, because an
extreme simulated value enlarges its own SD while the observed value never
does. This affects the shipped |t| too, by about 0.35 points at B = 400 on
synthetic nulls; on the real null |t| was still 5.0. Pooling the moments
makes the step-down exact under exchangeability and costs 0.4 to 0.8 points
of power. **Seed dependence:** at B = 400 about 10 percent of WY analyses
change with the seed, about 6 percent at B = 1000, and BH is less stable
than WY at both.

Whether to ship the split test, whether to pool moments in
`.bootstrap_pvalues()` for every cutoff function, and the default number of
iterations are open decisions. Full results and implementation notes are in
`HANDOFF-restscore-equal-tails.md`, the data in `restscore_eti_null.qmd`,
`restscore_power_conditional.qmd` and `restscore_stability.qmd`.

The sections below describe the original design.

Adds a parametric bootstrap null for the item-restscore test, and teaches
`RMitemRestscore()` to use it the same way `RMitemInfit()` uses
`RMitemInfitCutoff()`.

## Why

`dev/restscore_asymptotic_null.qmd` shows that the asymptotic p-value from
`iarm::item_restscore()` is miscalibrated under a true Rasch model, in a
direction that depends on the data:

- **Too many flags.** With BH over the items, 13.5% of null datasets get at
  least one item flagged for 20 dichotomous items at n = 150, and 30.5% when
  mistargeted by 1.5 logits.
- **Too few flags.** The same rate is about 1% for 9 polytomous items at
  n ≥ 1000.
- **Underfit is under-flagged in every cell.** Per-item rejection in that tail
  is 0.2 to 1.2% against a nominal 2.5%.

There are two causes, and a bootstrap has to handle both:

1. The expected gamma is treated as fixed. It is estimated from the same data
   and correlates .4 to .55 with the observed gamma, so the SD of
   observed − expected is 0.83 to 0.91 times the ASE, at every n.
2. The difference has a small-sample bias of order 1/n. It pushes z up by .15 at
   n = 150, which falls to .04 at n = 2000.

## Proposed API

### `RMitemRestscoreCutoff()`

```r
RMitemRestscoreCutoff(
  data,
  iterations = 400,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.95,
  dgp = c("resample", "conditional")
)
```

The signature and defaults are identical to `RMitemInfitCutoff()`, so the two
item-fit cutoff functions read as a pair.

**Algorithm per iteration**

1. Simulate a dataset from the CML thresholds, using `.wle_theta_pool()` and
   the same two DGPs as infit.
2. Apply the same validity checks as infit: every category observed, and at
   least 8 positive responses per dichotomous item.
3. **Refit by CML** with `psychotools::pcmodel(hessian = FALSE)`, then call
   `iarm::item_restscore(fit, p.adj = "none")`.
4. Store `Observed`, `Expected` and `Difference = Observed - Expected` per item.

Step 3 is the essential one. Holding the thresholds fixed at the generating
values would leave `Expected` constant across replicates and reproduce cause 1.

**Return value.** The same list shape as `RMitemInfitCutoff()`, so that
consumers and the reproducibility machinery need nothing new:

| element | content |
|---|---|
| `results` | `iteration`, `Item`, `Observed`, `Expected`, `Difference` |
| `item_cutoffs` | `Item`, `diff_low`, `diff_high` (HDCI or quantile of simulated `Difference`) |
| `actual_iterations`, `requested_iterations` | as infit |
| `sample_n`, `sample_n_total`, `sample_has_na`, `sample_summary` | as infit |
| `item_names`, `cutoff_method`, `hdci_width`, `dgp` | as infit |

**Missing data.** Use complete cases, as infit does. `iarm::item_restscore()`
already refits on complete cases when the fitted object contains `NA`, so the
observed statistic in `RMitemRestscore()` is complete-case too, and the two
agree without extra work.

### `RMitemRestscore()` gains a cutoff path

```r
RMitemRestscore(
  data,
  cutoff = NULL,
  p_value = NULL,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  output = "kable",
  sort,
  p_adj = "BH"
)
```

⚠️ Inserting `cutoff` second changes the positional order. Nobody calls
`RMitemRestscore(data, "dataframe")` positionally in the package or the jamovi
module, but a user script might. The alternative is to append the new arguments
after `p_adj`, which is safe but out of line with `RMitemInfit()`. Decision 4
below.

**Behaviour, mirroring `RMitemInfit()`:**

- `cutoff = NULL`: unchanged. The asymptotic iarm p-value is adjusted by
  `p_adj`, and the output is byte-identical to now.
- `cutoff` = full object, `p_value = NULL` resolves to `TRUE`. The following
  columns are added or replaced:
  - `Diff_low` and `Diff_high`, the interval, which is **descriptive**.
  - `p_restscore`, the marginal two-sided Monte Carlo p-value.
  - `padj_restscore`, the corrected p-value.
  - `Flagged`, which becomes `padj_restscore < alpha`.

  The iarm `p_adjusted` column is **dropped, not shown alongside**, the rule the
  gamma functions already follow.
- `cutoff` = the bare `$item_cutoffs` data.frame, or `p_value = FALSE`: flag on
  the interval, with the existing `.notify_band_flagging()` message.
- `p_adj` is ignored when p-values come from the bootstrap. If it is set to a
  non-default value in that case, warn once.

**The statistic being tested.** `Difference`, handed to the shared
`.bootstrap_pvalues(tail = "two.sided")`. It studentises by the bootstrap
**mean** and SD, so it removes the cause-2 bias and rescales for cause 1 without
any restscore-specific code. `Difference` is preferred over the ASE-based z for
three reasons:

- The bootstrap does the scaling, so the ASE adds nothing.
- `Difference` is the column users already see.
- It avoids the infinite z that an ASE of 0 produces when observed gamma = ±1,
  which happened twice in 80,000 null item-fits.

**The direction label comes from the studentised residual, not the raw sign.**
The null mean of `Difference` is positive at small n, so an item can sit
slightly above zero and still be significantly *below* its null mean. `Flagged`
is `"overfit"` when `Difference` exceeds the bootstrap mean and `"underfit"`
when it falls below. Document this, because it is the one place where the
restscore table can look inconsistent (a flagged "underfit" item with a
positive `Difference` is possible at small n, though rare).

**Caption.** The same elements as infit:

- iterations and `dgp`
- correction label via `.correction_label()`
- the `< 1000` reproducibility note
- the `n = X of Y` sample line

**Consistency checks.** Item names in `cutoff$results` must match `data`, the
same check infit uses. Also warn when `cutoff$sample_n` differs from the
complete-case n of `data`, which catches a cutoff object computed on another
dataset. Infit does not check this yet, so either add it to both or to neither.
Decision 5.

### Plot (phase 2)

`RMitemRestscorePlot(simfit, data)` follows the `<base>Plot` convention of its
siblings:

- per-item simulated `Difference` distribution
- the observed value as a point
- the descriptive band

This is optional for the first release. The `Boot` function has no plot either.

## What does not change

- `RMitemRestscoreBoot()`. It answers a different question (how stable is a
  flag under resampling), and uses the asymptotic test inside.
- **Proposed change to the `RMitemRestscoreBoot()` docs:** the sentence saying
  the asymptotic test "can flag items that are not practically misfitting" at
  large n is about power, not calibration, and the null study shows the test
  is *conservative* at large n. Worth rewording in the same release. Decision 6.
- jamovi. Out of scope until the package side is settled.

## Validation before release

The bootstrap should fix both causes by construction, but that is a claim to
test, not to assume. A nested null study reuses `restscore_asymptotic_null.qmd`:

- **Cells.** The worst and the most conservative cells from the null study:
  - dichotomous, μ = 1.5, n = 150
  - dichotomous, μ = 0, n = 300
  - polytomous, μ = 0, n = 150
  - polytomous, μ = 0, n = 1000
- **Arms.** `dgp = "resample"` and `"conditional"`, which also settles
  decision 1.
- **Per cell.** 500 null datasets, each with `RMitemRestscoreCutoff(iterations
  = 400)`. Record FWER from WY step-down and BH, and per-item rates in each
  tail.
- **Target.** Family-wise error within Monte Carlo error of .05 (±1.9 points at
  500 datasets), and symmetric tails.
- **Cost.** About 0.1 s per iteration gives 40 s per dataset and 5.6 CPU-hours
  per cell-arm. The full study is about 45 CPU-hours, or 4.5 hours on 10 cores.
  Halve it by running a single arm if decision 1 is taken on theory.

Power against infit with cutoffs is a separate, later study. (Done 2026-10-01
in `restscore_power_conditional.qmd`. The symmetric-tails target was not met
by |t|; see the 2026-10-03 revision at the top.)

## Performance

`iarm::item_restscore()` computes the expected gamma with R loops over
`pscore_poly()`, and that dominates the cost of one iteration. 400 iterations
take about 40 s on one core, or under 10 s in parallel. That is acceptable. A
vectorised expected-gamma, like the partial gamma one, is a later optimisation
and would need its own validation against iarm.

**Alternative considered: extending `RMitemInfitCutoff()`.** The same
simulated dataset and CML refit could also return the restscore statistics,
saving a refit per iteration. This is rejected for now for two reasons:

- It changes an existing function's return shape.
- It couples two statistics whose defaults (`dgp`, width) may need to diverge
  after validation.

It can be revisited if users routinely run both.

## Decisions for the user

1. **`dgp` default.** Keep `"resample"` for sibling consistency, and let the
   validation study decide whether `"conditional"` is better. Recommended.
   The case for `"conditional"`: the expected gamma is computed from the
   observed score distribution, and that DGP holds the score distribution
   fixed, so it is the matched null in the same sense as for conditional infit.
2. **Validate first, or ship with a caption note.** Recommended: implement,
   then run the validation, then release. Nothing ships unvalidated, following
   the infit precedent.
3. **The asymptotic default path.** Leave the numbers alone, but add one
   sentence to the docs, and possibly the caption, saying the asymptotic
   p-value is miscalibrated in both directions, pointing to
   `RMitemRestscoreCutoff()`. Recommended. Changing the default flagging is a
   separate decision, like Release B for local dependence.
4. **Argument order in `RMitemRestscore()`.** `cutoff` second (matches
   `RMitemInfit()`, small breaking risk) or new arguments last (safe,
   inconsistent). Leaning towards matching infit, with a NEWS line.
5. **The sample-size mismatch check.** Add it to both item-fit consumers or to
   neither.
6. **Reword the `RMitemRestscoreBoot()` large-sample sentence** in the same
   release.
7. **Plot now or in phase 2.**

## Files touched

- New: `R/item_restscore_cutoff.R` with runners modelled on
  `run_single_infit_sim()`, and
  `tests/testthat/test-item_restscore_cutoff.R`.
- `R/item_restscore.R`: the cutoff path, docs and caption.
- `R/zzz_reproducibility.R`: register the new function if the reproducibility
  page lists the bootstrap functions.
- `NEWS.md`, `_pkgdown.yml` (add to the "Item Restscore" section) and
  `README.md` (function list).
- Tests: the new file plus `test-item_restscore.R`, including one check that
  the `cutoff = NULL` output is unchanged.
