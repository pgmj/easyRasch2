# Slope ROPE: equivalence decisions on item discrimination

Sketch written 2026-10-04 in the large-sample item fit thread. Nothing has
been run. Background in `conditional-rmsd-design.md` (section 2, the infit
identity) and `easyRaschBayes-fitstat-notes.md`. Nothing here relies on the
superseded `sx2_itemfit_sim.qmd`. The dichotomous generator follows
`sx2_comparison_v2.R`.

## Idea

For dichotomous items, conditional infit is a linear slope contrast of the
conditional ICC residuals. Its parametric counterpart is the item
discrimination. Estimate each item's discrimination relative to the common
level, and decide with a region of practical equivalence (ROPE) instead of a
test against exact Rasch fit.

- **Parameter.** `delta_i = log a_i - mu`, the item's log-discrimination
  relative to the common level. Zero is Rasch-consistent in slope. Negative
  means underfit (flatter) and positive means overfit (steeper).
- **ROPE.** `|delta_i| <= r` on the log scale, symmetric in ratio terms (for
  example `r = log(1.25)`, so `a_i / exp(mu)` lies in [0.80, 1.25]).
- **Decision (Kruschke, HDI+ROPE).** 95% HDI fully inside the ROPE: the item
  is practically Rasch in slope. Fully outside: relevant slope misfit.
  Otherwise: undecided.

Why this could address the large-sample problem: the posterior concentrates
on the population value as n grows. Trivial misfit therefore ends up
decisively *inside* the ROPE instead of being flagged. A test does the
opposite.

What it does not address: the ROPE width is the same `delta` question as
before, now explicit. Misfit that is not a slope (shared residual structure,
threshold problems) only shows up through its slope projection.

## Model

- **Bayesian (easyRaschBayes route).** One family for both formats: a
  dichotomous item is a two-category `acat`.

      bf(resp | thres(gr = item) ~ 1 + (1 | id),
         disc ~ 1 + (1 | item))

  The formula parses and gives the expected prior table with brms 2.23
  (checked 2026-10-04). It has not been fitted yet. `disc` uses a log link,
  so the item effects on `disc` are `delta_i` directly. **Identification:**
  the person SD and the `disc` intercept cannot both be free. Either fix the
  `disc` intercept at 0 (`prior(constant(0), class = Intercept, dpar = disc)`)
  and estimate `sd(id)`, or fix `sd(id) = 1` and estimate the intercept as
  `mu`. The first makes `delta_i` the item random effect.
- **Pooling.** Hierarchical item effects shrink `delta_i` toward 0, which is
  toward the ROPE. That is a risk at small n and with few items (k = 9).
  Compare this with fixed item effects under a sum-to-zero constraint.
- **Frequentist twin.** An MML GPCM/2PL (`mirt`) with Wald or profile CIs on
  `log a_i - mean(log a)`, using the same inside/outside/undecided rule. This
  is TOST equivalence with a 90% CI, or the 95% version to match the HDI. At
  large n with weak priors the two should agree (Bernstein-von Mises). The
  MML version is fast enough for thousands of replicates. It would also make
  an easyRasch2 implementation possible without Stan.

## Steps before any study

1. **Identification and recovery.** Fit the brms model to one simulated
   dataset per format with known slopes and check recovery of `delta_i`.
   Time the fit at n = 500, 2,000 and 5,000. The guess is tens of minutes at
   n = 5,000 with 20 items (100,000 observations).
2. **Slope-to-infit mapping.** At large n, the expected conditional infit
   for each relative slope on the grid below, for both formats. This
   translates a ROPE into infit units and back. For reference, one
   verification dataset gave infit 1.27 for slope 0.5 with 12 dichotomous
   items. The phq9 calibration gave infit 0.79 to 0.80 for slope 1.6.
3. **Pseudo-true `delta` for out-of-model misfit.** Fit the MML 2PL/GPCM at
   n = 200,000 to second-dimension data (rho .50, .75, .85) to find the
   slope value the misfit projects to. Decisions under that DGP are judged
   against this pseudo-true value.

## Study

### Arm 1: slope misfit (the model is correct)

- One target item per dataset with true relative slope `exp(delta)` in
  {0.5, 0.67, 0.80, 0.9, 1, 1.1, 1.25, 1.5, 2}, other items at 1.
  Includes the boundary values 0.80 and 1.25 for `r = log(1.25)`.
- A cascade condition: three items at 0.5 (or 2.0), with decisions recorded
  for all the fitting items.
- Formats: 20 dichotomous items (evenly spaced on [-2, 2], theta SD 1.5),
  and 9 four-category items (phq9-like thresholds from
  `restscore_power_generator.R`). Target items on and off target.

### Arm 2: second-dimension misfit (the model is wrong)

- The `sx2_comparison_v2.R` and phq9 generators, rho .50, .75 and .85, one or
  three items. Decisions are judged against the pseudo-true `delta` from
  step 3.

### Sample sizes and engines

- n = 500, 2,000, 5,000, and 10,000 for the MML engine only.
- MML: 500 replicates per cell for the operating characteristics.
- brms: 30 to 50 replicates in a subset of cells (both formats, n = 500 and
  5,000, slopes 0.5, 0.80, 1, 1.25), to check that HDI+ROPE matches the MML
  rule and to measure shrinkage. Runtime decides the final size.

### Outcomes

1. **Operating characteristics.** P(inside), P(outside) and P(undecided) as
   functions of the true `delta` and n. At a boundary, P(outside) should be
   at most about 2.5%. Well inside the ROPE, P(inside) should go to 1 as n
   grows. That is the large-sample property this approach is meant to
   deliver.
2. **Cascade.** How often a fitting item is declared outside the ROPE when
   one or three items misfit. On a subset, compare with conditional infit
   under Westfall-Young at the same n.
3. **Shrinkage.** Bias of `delta_i` for hierarchical against fixed effects,
   by n and k.
4. **Undecided share.** How much data an inside decision needs at each
   `delta`. This is the practical sample-size guidance.
5. **Runtime** of the brms route at survey sizes.

### Predictions

- P1. HDI+ROPE and MML-TOST give the same decisions in at least 95% of
  replicates at n of 2,000 or more.
- P2. The cascade appears only as a small shift in `mu`, about 1/k of the
  misfit. Fitting items stay inside a ROPE of `log(1.25)`, even at n =
  10,000.
- P3. Hierarchical shrinkage makes inside decisions too easy at n = 500
  with k = 9.

## Open decisions

1. **ROPE width.** `log(1.25)` is a placeholder. Candidates are a fixed
   ratio band, the infit mapping from step 2, or an impact anchor. The
   linear-composite result (keeping an item beats removing it when its
   relative slope exceeds about 0.5) suggests two zones: "practically
   Rasch" inside about [0.8, 1.25], and "worth removing" below about 0.5.
2. **Engine.** Whether MML (easyRasch2) or brms (easyRaschBayes) is the
   primary target. Both are possible, and the study tells whether they
   agree.
3. **Hierarchical or fixed item effects.**
4. **Polytomous scope.** A GPCM slope captures only discrimination misfit.
   Disordered thresholds and category problems need separate parameters.
5. Ask before running anything (memory: `ask-before-running-simulations`).
   Steps 1-3 above are short and would come first.
