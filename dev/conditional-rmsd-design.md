# Conditional RMSD for Rasch item fit: literature and study design

Written 2026-10-04. Follows `HANDOFF-large-sample-fit.md`. Nothing here has
been run.

**Revised 2026-10-05.** An earlier version leaned on a matched-power
reanalysis of the cache of `sx2_itemfit_sim.qmd`. That study is superseded
by `sx2_null_diagnostic.qmd` and `sx2_comparison_v2.qmd` (shared master
dataset, mirt's S-X2 df one too small, infit without multiplicity control),
so the reanalysis is withdrawn and nothing here rests on it. The Part A
design in section 4 is superseded by `conditional_rmsd_partA.qmd`, which
uses the v2 generator, sample sizes up to 20,000, two more conditions, and a
Holm reanalysis of the v2 cache.

## 1. What the RMSD literature says

Nine papers read in full (`articles/RMSD/`).

1. **Definition.** Every RMSD in this literature is marginal. The observed
   item response function comes from MML posterior pseudo-counts on a
   quadrature grid, weighted by the estimated normal theta density (Köhler
   et al., 2020; Robitzsch, 2022; Kim et al., 2026). Almost all work is on the
   2PL, dichotomous items, and misfit that is either a wrong functional form
   (guessing, non-monotone, plateaus) or uniform DIF. No paper defines an
   RMSD conditional on the sum score under CML, and Robitzsch (2025, 2026)
   states that how to extend the RMSD to polytomous items is unresolved.
2. **Positive bias.** The sample RMSD overestimates the population value
   because sampling error enters squared. For fitting items with 50 items,
   the mean is .038 at n = 500 and .012 at n = 5,000 (Köhler et al., 2020,
   Table 1). The bias grows with test length. Nonparametric bootstrap
   correction removes only part of it (Köhler et al., 2020). An analytic
   correction on the squared scale (`CDM`) and Taylor-based corrections on
   the root scale (Robitzsch, 2022, 2025, who recommends his R4 and R5)
   remove most of it at the cost of a slightly larger SD.
3. **Testing.** A parametric bootstrap null gives Type I error of .06 to .07,
   rising to .10 to .24 when 30 percent of items misfit at n = 5,000 (Köhler
   et al., 2020, Table 4). Poor person's PPMC is conservative (Kim et al.,
   2026). Köhler et al. report that the bootstrap RMSD beats infit, outfit
   and S-X2, but their infit was TAM's posterior-based infit under a 2PL with
   a fixed 1.15 cutoff, which detected nothing. That comparison says nothing
   about conditional infit with simulated cutoffs.
4. **Minimum-effect testing already exists.** Robitzsch (2025, PISA example)
   flags an item when the lower CI bound exceeds a cutoff (0.05 to 0.15).
   Robitzsch (2026, preprint) builds CIs from a parametric bootstrap of the
   residuals and recommends percentile intervals because the RMSD is skewed.
   Two limits: the CI covers the n-dependent pseudo-true value, not the
   population value, and item parameters are held fixed. He says that
   single-group item fit with estimated parameters, which is our case, needs
   the method adapted.
5. **The cascade appears here too.** Misfit leaks into fitting items through
   the posterior. Population RMSD of fitting items rises, and that of
   misfitting items falls, as the share of misfitting items grows. Robitzsch
   (2022) concludes that the RMSD "must always be interpreted as a relative
   fit statistic".
6. **Slope misfit gives small values.** Under a 1PL, an item with
   discrimination 0.6 against 1.0 has population RMSD .027 (Robitzsch, 2022,
   Table 2), well under the common .05 cutoff. In a Rasch analysis, low slope
   is the main signature of underfit, including misfit from a second
   dimension.
7. **Targeting.** Distribution weighting shrinks the RMSD for items far from
   the bulk of the respondents, and asymmetrically by the direction of the
   shift (Tijmstra et al., 2020). Weighting centred on the item location
   removes this but adds bias and variance (Robitzsch, 2025).
8. **Cutoffs.** Fixed cutoffs (.12 PISA, .15 PIAAC) miss moderate misfit.
   Köhler et al. (2020) label population RMSD below .02 negligible, .02 to
   .05 small, .05 to .08 medium, and .08 or more large. Kim et al. (2026)
   derive ROC-optimal cutoffs that depend on n and test length (for example
   .013 at n = 5,000 with 20 items) and give a response-surface formula.
   Those cutoffs depend on the misfit types simulated, which assume a
   population RMSD of at least .05.
9. **Relative flagging.** von Davier and Bezirhan (2023) flag items whose
   RMSD lies more than 2.5 MADs above the item median. This answers "which
   items stand out", not "do items fit". False-positive rates were .06 to
   .25 per item with 30 items, and the MAD is unstable for short scales.
10. **Practical significance.** Sinharay and Haberman (2014) judge misfit by
    changes in equating and pass/fail decisions after removing misfitting
    items. They compare against removing the same number of *random* items,
    which isolates the effect of misfit from the effect of a shorter test.
    Köhler and Hartig (2017) bound the change in a correlation with any
    covariate, using a two-dimensional Rasch model with fitting and
    misfitting items on separate dimensions.

Correction to the chat on 2026-10-04: the paper that rescaled PISA 2018
field-trial data with and without misfitting items is the 2022 *Large-scale
Assessments in Education* article, not Köhler and Hartig (2017). Its PDF is
not in the folder.

## 2. Conditional RMSD under the Rasch model

### Definition

For item `i` with maximum score `m_i`, let `s` be the total score, `n_s` the
number of respondents with score `s`, `N` the number with a non-extreme total
score, `O_is` the observed mean item score in group `s`, and `E_is`, `V_is`
the conditional mean and variance from the ESF (`.cond_moments()`). Weight
by the observed score distribution, `w_s = n_s / N`, over non-extreme scores:

    MSD_i  = sum_s w_s (O_is - E_is)^2 / m_i^2
    cRMSD_i = sqrt(MSD_i)

Dividing by `m_i` puts dichotomous and polytomous items on the scale of
proportion of the maximum item score. This is the numeric summary of what
`RMitemICCPlot()` draws. It needs no theta distribution, no quadrature grid
and no normality assumption.

### Three properties that follow from CML

1. **The signed mean deviation is zero.** The CML score equation forces
   `sum_s n_s (O_is - E_is) = 0`, so the conditional MD is 0 by construction
   (verified numerically, section 3). A uniform shift is absorbed into the
   item location. What remains is the shape of the conditional ICC.
2. **The bias has a closed form.** Under the model, responses within a score
   group are independent given `s`, so `E[(O_is - E_is)^2] = V_is / n_s` and

       B_i = sum_s w_s V_is / (n_s m_i^2) = sum_{s in G} V_is / (N m_i^2)

   where `G` is the set of observed non-extreme scores. The bias grows with
   the number of score groups, which explains why the RMSD bias grows with
   test length. Estimating the item location removes about one degree of
   freedom, so the true bias is slightly below `B_i`. This needs checking
   by simulation. The bias-corrected value is `sqrt(max(MSD_i - B_i, 0))`.
   For 20 dichotomous items the raw null value is roughly
   `sqrt(19 * 0.15 / N)`: about .024 at N = 5,000 and .076 at N = 500. That
   is larger than the marginal RMSD bias, because sum-score groups are
   finer than smoothed posteriors. So the statistic is noisy at small n,
   and pooling adjacent scores is a possible fix. The large-n case is the
   target, so this is acceptable.
3. **Conditional infit is a linear slope contrast of the same residuals.**
   For dichotomous items, the within-group variance is fixed by the group
   mean, and algebra gives

       InfitMSQ_i - 1 = sum_s n_s (O_is - E_is)(1 - 2 E_is) / sum_s n_s V_is

   This is exact. It was verified against `iarm::out_infit()` on simulated
   data to four decimals for all 12 items. Infit is a *signed, linear*
   contrast: above 1 when the observed curve is above expectation at low
   scores and below it at high scores, which is a flatter curve. The
   cRMSD is a *quadratic, omnibus* summary of the same residuals. For
   polytomous items, infit also contains within-group dispersion, so the
   identity does not hold.

### What property 3 implies

- Infit MSQ is linear in the residuals, so it has no squared-error bias. Its
  expected value for a given misfit settles at a population value as n
  grows. **Infit MSQ is already an n-stable effect size**, at least for
  dichotomous items. In large samples, the problem is that a test against
  MSQ = 1 rejects trivial values. A new statistic is not obviously needed.
- Infit is a test against the alternative a second dimension produces (a
  flatter curve). So it should be more powerful than cRMSD against that
  misfit, as it was against S-X2 in `sx2_comparison_v2.qmd`. cRMSD should
  win only against shapes that
  are orthogonal to the `(1 - 2E)` contrast (non-monotone, plateau, or a
  curve steep in one half and flat in the other).
- What cRMSD offers is units: "the observed conditional curve deviates from
  the model by about .04 on average" is readable. MSQ = 1.15 is not.
- The cascade shows up as overfit in the other items (0.94 to 0.97 with one
  flat item in the verification data). The cascade is a population effect,
  not an artefact of n. If the leakage stays well below the misfitting
  item's value on an effect-size scale, a floor `delta` between the two
  stops the cascade at any n, for either statistic.

## 3. Verification already done

`devtools::load_all()`, 12 dichotomous items, n = 2,000, item 3 generated with
slope 0.5. Infit from the identity, from direct summation and from
`iarm::out_infit()` agreed to four decimals for all items. Conditional MD
was 0.0000 for all items. Item 3 infit 1.27, the others 0.94 to 1.03.

## 4. The study

### Questions

- Q1. Null: does bias-corrected cRMSD centre near 0 at every n? Does `B_i`
  match the bootstrap mean of MSD?
- Q2. Stability: under misfit, does bias-corrected cRMSD of the misfitting
  items settle at its population value as n grows, as infit MSQ should?
- Q3. Leakage: how large are the values on fitting items when 1 or 3 items
  misfit, relative to the misfitting items, on each scale (cRMSD, MSQ,
  gamma difference, RMSEA.S_X2)? The ratio at population level decides
  whether an effect-size floor can stop the cascade.
- Q4. Discrimination: at matched detection of the misfitting items, which
  statistic flags the fewest fitting items?
- Q5. Targeting: how do values for the misfitting items at 0, -1 and -2
  logits compare across statistics?

### Predictions, stated before running

- P1. Infit at least matches cRMSD on Q4 under second-dimension misfit, and
  beats it at rho = .75.
- P2. Bias-corrected cRMSD and infit MSQ are both stable in n (Q2). Raw cRMSD
  and the gamma difference are not.
- P3. Leakage is small relative to misfit at rho = .50 and closer at
  rho = .75 (Q3).
- P4. cRMSD falls off faster than infit for the -2 logit item (Q5).

### Part A: statistics only, no bootstrap

Cheap. It answers Q1 to Q5 descriptively and through ROC.

- **DGP:** the dichotomous design of `sx2_comparison_v2.R` (20 items evenly
  spaced on [-2, 2], theta SD 1.5, second dimension generated directly,
  misfitting items nearest 0, -1 and -2 logits). Fresh datasets each
  replicate, not subsamples of a master dataset.
- **Conditions:** k = 0 (null), and k = 1 and 3 at rho = .50 and .75.
  Five conditions.
- **n:** 500, 1,000, 2,000, 5,000. 200 replicates per cell (4,000 datasets).
- **Population values:** one dataset of n = 200,000 per condition, with the
  same statistics, for the Q2 and Q3 targets.
- **Recorded per item and dataset:** raw and bias-corrected cRMSD
  (distribution-weighted, unpooled; primary), plus secondary variants
  (pooled to at least 50 per score group, and information-weighted
  `w_s ∝ n_s V_is`). Also conditional infit MSQ, the signed direction of the
  cRMSD residual slope, the item-restscore gamma difference (observed minus
  expected), and S-X2 with RMSEA from `mirt`. All are continuous, so ROC
  curves can be drawn for every statistic.
- **ROC scoring:** misfitting against fitting items within the same
  datasets. Use |MSQ - 1| and |gamma difference| so all statistics are
  direction-free like cRMSD. Report the false-positive rate at matched
  detection and partial AUC at FPR of .05 or less.
- **Estimated cost:** about 20 to 30 minutes on 10 cores. `mirt` at
  n = 5,000 dominates.

### Part B: bootstrap tests, smaller

- Conditional parametric bootstrap (the `dgp = "conditional"` machinery),
  with a CML refit per replicate. In each replicate compute cRMSD and infit
  from the same refit, so the two tests are paired. One-sided upper test for
  cRMSD, two-sided for infit, Westfall-Young step-down for both.
- Conditions: null and rho = .75 (k = 3), n of 500, 2,000 and 5,000,
  200 replicates, B = 400.
- Outcomes: familywise error under the null, detection and false flags
  under misfit.
- **Estimated cost:** about 1.5 to 2 hours on 10 cores.

### Checks before the main run

1. `B_i` against the mean null MSD from 1,000 conditional simulations at
   n = 500 and 5,000.
2. Conditional MD equals 0 on every simulated dataset (sanity check on the
   moments).
3. cRMSD from score-group means against an independent computation from
   `RMitemICCPlot()`'s internal curve data for one dataset.
4. A small pilot (n = 1,000 and 5,000, 50 replicates, Part A only) to time
   the run and look at the distributions before committing.

## 5. Open decisions

1. **Weighting.** Distribution weighting (primary here) measures the
   deviation where respondents are. That fits a population-impact reading
   but loses sensitivity off target (Tijmstra et al., 2020).
   Information-weighted or location-centred weighting is the alternative.
   Both are recorded in Part A, but the package would need one.
2. **Polytomous definition.** Expected-score deviation over `m_i` (proposed),
   or category-wise deviations. The identity in section 2 does not carry
   over, so the polytomous arm (phq9 generator from
   `restscore_power_generator.R`, underfit rho = .85 and overfit slope 1.6)
   is a second phase.
3. **Misfit types.** Part A uses only second-dimension misfit, where infit
   should win. A fair test of cRMSD needs a shape infit cannot see (for
   example a non-monotone or plateau item). Add it to Part A, or leave it
   for phase 2?
4. **Where this leads.** If P1 to P3 hold, the large-sample answer may be
   "conditional infit plus an effect-size floor", with cRMSD as the readable
   companion and `delta` set from impact (option D in the handoff).
   Robitzsch's lower-CI-bound rule would then need adapting to estimated
   item parameters and the population (not pseudo-true) value.

## References (all in `articles/RMSD/`)

- Kim, Y.-K., Cai, L., & Kim, Y. (2026). *Educational and Psychological
  Measurement, 86*(1), 156-190.
- Köhler, C., & Hartig, J. (2017). *Applied Psychological Measurement,
  41*(5), 388-400.
- Köhler, C., Robitzsch, A., & Hartig, J. (2020). *Journal of Educational and
  Behavioral Statistics, 45*(3), 251-273.
- Robitzsch, A. (2022). *Foundations, 2*, 488-503.
- Robitzsch, A. (2025). *Foundations, 5*, 36.
- Robitzsch, A. (2026). Confidence interval estimation for RMSD and MD item
  fit statistics. Preprints.org, doi:10.20944/preprints202603.1174.v1
  (not peer reviewed).
- Sinharay, S., & Haberman, S. J. (2014). *Educational Measurement: Issues
  and Practice, 33*(1), 23-35.
- Tijmstra, J., Bolsinova, M., Liaw, Y.-L., Rutkowski, L., & Rutkowski, D.
  (2020). *Journal of Educational Measurement, 57*(4), 566-583.
- von Davier, M., & Bezirhan, U. (2023). *Educational and Psychological
  Measurement, 83*(4), 740-765.
