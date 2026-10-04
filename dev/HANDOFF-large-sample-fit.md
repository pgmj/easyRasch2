# Handoff: item fit in large samples, and `RMitemRestscoreBoot()`

Written 2026-10-04 at the end of the equal-tailed restscore thread
(`HANDOFF-restscore-equal-tails.md`). Nothing here is built yet. This is a
starting point for a new chat.

## Why this came up

`RMitemRestscoreBoot()` is documented and recommended as the tool for large
samples (n > 1000), where "the asymptotic test underlying `RMitemRestscore()`
can flag items that are not practically misfitting". Each bootstrap
iteration runs that same asymptotic test. Since 2026-09-30 we know the
asymptotic test is miscalibrated (`restscore_asymptotic_null.qmd`), so the
method inherits the problem it is meant to work around.

## What `RMitemRestscoreBoot()` does (`R/bootstrap_restscore.R`)

- Draws `iterations` (default 200) samples of size `samplesize` (default
  600) **with replacement** from the data, refits the PCM by CML, and runs
  `iarm::item_restscore()` with its default BH adjustment across items.
- Classifies an item as overfit or underfit when the BH-adjusted p < .05,
  by the sign of expected − observed.
- Reports the percentage of iterations in which each item was flagged, with
  full-sample infit and relative location alongside.

## What is wrong with it

1. **The inner test is miscalibrated** (`restscore_asymptotic_null.qmd`,
   under a true Rasch model):
   - Too many flags for dichotomous items in small or mistargeted samples:
     13.5 percent of datasets with a BH flag for 20 items at n = 150, 30.5
     percent when 1.5 logits off target, almost all overfit.
   - Too few flags for polytomous items at large n: about 1 percent for 9
     items at n ≥ 1000.
   - Underfit under-flagged in every cell: 0.2 to 1.2 percent per item
     against a nominal 2.5.
   - Causes: the SE ignores that expected gamma is estimated from the same
     data (it correlates .4 to .55 with observed), and the difference has an
     O(1/n) bias.

   What matters is the calibration **at `samplesize`, not at the full n**,
   because each iteration tests a sample of that size. At the default 600,
   the test still under-flags underfit, so a low percentage for underfit is
   weak evidence of fit.
2. **The percentage is not a probability of misfit.** The docs say
   bootstrapping gives "a more nuanced view of the probability of an item
   actually being misfit". What it estimates is how often a miscalibrated
   test, at n = `samplesize`, flags the item when samples resemble the
   observed data. It is closer to an estimate of that test's power at that
   sample size than to a probability that the item misfits.
3. **The rationale is half wrong.** At large n the asymptotic test is
   conservative for polytomous items, not liberal, so "can flag items that
   are not practically misfitting" describes power, not calibration. The
   design record (`restscore-cutoff-design.md`, decision 6) kept that
   sentence on the grounds that it is about practical versus statistical
   significance. That still holds for a *calibrated* test, whose power goes
   to 1 at large n for any misfit, however small.
4. **Resampling with replacement at `samplesize` < n** is an m-out-of-n
   bootstrap. Duplicated respondents make each sample less informative than
   a real sample of that size. Subsampling without replacement would give
   samples that behave like real ones of size `samplesize`.
5. **`samplesize` drives the result** and has no principled default. 600 is
   arbitrary, and the percentages move with it.

## The underlying problem to solve

With large n, a calibrated test (infit or restscore with
`RMitemInfitCutoff()` / `RMitemRestscoreCutoff()`) detects misfit too small
to matter. Users need a way to judge whether misfit is **large enough to
matter**, not whether it exists. The same problem applies to conditional
infit, so a solution should cover both item fit statistics.

## Candidate directions (none tested)

A. **Calibrated subsample flag rate.** Keep the idea of `RMitemRestscoreBoot()`
   but replace the inner test with the calibrated one: subsample without
   replacement at size n_s, and test each subsample against a parametric
   null *at n_s*. Running `RMitemRestscoreCutoff()` per subsample is too
   slow (about 40 s each at 400 iterations). A cheaper version builds one
   null at n_s from the full-sample fit and reuses it for every subsample.
   Its validity needs checking, since the null depends on n and on the
   estimated parameters. The output would answer "how often would this item
   be flagged in a well-calibrated study of size n_s?". This is still a
   power statement, but an honest one.

B. **Effect size with an interval, against a practical threshold.** Report
   the difference (restscore) or the MSQ (infit) with a bootstrap CI, and
   test against a minimum relevant effect (H0: |effect| ≤ δ) instead of
   against zero. The hard part is justifying δ. The descriptive `Diff_low` /
   `Diff_high` band of the cutoff functions is related, but it shrinks with n,
   so it is not a practical threshold.

C. **Impact-based judgement.** Judge misfit by its consequences for
   measurement instead of its size: the change in person locations,
   reliability or targeting when the item is removed or modelled
   differently. This links to the planned DIF impact function
   (memory: `easyrasch2-dif-impact-function-idea`) and the LD magnitude
   versus impact distinction (memory: `ld-magnitude-vs-impact`). It probably
   gives users the most useful answer, and is the most work.

D. **Calibrate δ by simulation.** A hybrid of B and C: find the size of
   gamma difference or infit deviation at which person estimates or
   reliability change by a stated amount, and use that as δ.

Before choosing, check the literature for established large-sample
approaches in Rasch item fit: sample-size-adjusted fit statistics,
subsampling recommendations, and effect sizes for item-restscore gamma.
Verify each source before citing (memory:
`verify-before-reporting-simulation-findings`, `reviews-psychometric-articles`).

## Documentation fix (applied 2026-10-04)

Applied to the help page, NEWS.md and README.md on 2026-10-04 as described below, and to the jamovi bootrestscore analysis (footnote, description, NEWS). Only the optional kable caption line is not done. It changes no results.

- **`RMitemRestscoreBoot()` help:** replace the "Useful with large samples"
  paragraph. Say that each iteration uses the asymptotic test with BH, and
  that this test is miscalibrated at the subsample size, pointing to the
  `RMitemRestscore()` help. Say what the percentage is (how often that test
  flags the item in samples of size `samplesize`) and is not (a probability
  of misfit). Say that a low underfit percentage is weak evidence of fit, and
  that the result depends on `samplesize`. Recommend `RMitemRestscoreCutoff()`
  for testing, with this function as a descriptive check of how stable the
  asymptotic flags are.
- **README line 99:** leave as is, or add "(descriptive; see its help page)".
- **Vignette:** does not mention `RMitemRestscoreBoot()`, so nothing to change.
- **jamovi:** has a `bootrestscore` analysis that calls this function. Its
  notes would need the same caveat, asking first (memory:
  `ask-before-crossing-package-boundary`).
- **Optional:** a caption line in the kable output with the same caveat. That
  is a code change, not only docs.

## Infrastructure that can be reused

- `RMitemRestscoreCutoff()` with `dgp = "conditional"`: a calibrated null at
  the observed n.
- `restscore_asymptotic_null.qmd`: the asymptotic test's null rates by cell.
- `restscore_power_conditional.qmd` and its generator
  `restscore_power_generator.R`: planted underfit (second dimension) and
  overfit (slope) at n = 150 to 600. Large-n cells would need adding.
- Cached simulated matrices in `restscore_split_null/`,
  `restscore_power_split/` and `restscore_power_split_ext/`.

## Working constraints

- Lead with theory and ask before running simulations (memory:
  `ask-before-running-simulations`). Cap at 10 cores.
- Do not change existing package behaviour without confirming first.
- Run the doc review checklist for any exported-function change.
