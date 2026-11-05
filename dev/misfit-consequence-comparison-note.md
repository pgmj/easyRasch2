# Consequences of misfit: comparing candidate scales

Discussion note written 2026-10-14. No code. It follows the decision to
pause the item impact validation run (`item-impact-validation-design.md`,
section 14) and reframes the large-sample item fit problem around
decisions between candidate scales instead of per-item impact.

## 1. Why reframe

The impact interval work validated the easy case: one strongly misfitting
item, which is visible at any sample size and leads to an obvious decision.
The difficult cases are different.

- **The residual phase.** After the gross misfit is removed, many items
  deviate by .03 to .10 in MSQ. At n of 2,000 and above all of them are
  statistically significant. The question is whether the remaining misfit,
  taken together, changes anything that matters.
- **Attribution.** Misfit is relative. Conditional infit compares an item
  with the score formed by all items, so the data show that the items do
  not form a Rasch scale but not which items are the construct. A cascade
  and an LD pair that makes other items look underfit are the same
  phenomenon.

The cached dependence arm shows the second case. One overfitting pair in
phq9 (uniform response dependence, λ = 1) gives the dependent item MSQ .83,
its partner .96, and every other item 1.01 to 1.05, although all of them
fit in truth. At n = 5,000 the FWB interval is roughly ±.04, so the most
affected of these would be flagged as underfitting. In d20 the effect is
diluted (others .99 to 1.03 at λ = 2).

Per-item impact rests on the counterfactual "this item fitting, everything
else unchanged". With interdependent misfit that is not an option anyone
can choose. The options are scale-level: keep all items, remove one of
several candidates, merge an LD pair into a testlet, collapse categories.

## 2. Sinharay and Haberman (2014)

*How often is the misfit of item response theory models practically
significant?* Educational Measurement: Issues and Practice, 33(1), 23-35.

**Definition.** Practical significance of misfit is the extent to which
decisions made from the test scores are robust against the misfit. The
ideal comparison would be with decisions under a perfectly fitting model,
which does not exist. They replace it with a "model-data combination with
better fit", obtained by a more general IRT model, exclusion of a few
misfitting items, combining unpopular score categories of polytomous items,
or exclusion of a few examinees. If decisions under the operational and the
better-fitting combination agree, the misfit is not practically
significant. Their Figure 2 lays this out as a flow chart.

**What decides the comparison.** Two pieces of information: how the model
is used to compute scores (equating, pattern scoring), and how the scores
are used (pass or fail, proficiency classification). The consequence
measure follows from the use.

**Examples.**

1. English proficiency test (2PL and GPCM, n about 1,500, 28 anchor items).
   Five misfitting anchor items removed. The difference in equating
   functions stayed below the difference that matters (0.5 raw-score
   points, Dorans and Feigenbaum) except at very low scores. Pass or fail
   disagreement over 1,000 random cut scores was 0.0003 percent.
2. State test (1PL operational against a restricted 3PL with a common
   guessing parameter, n 31,000 to 75,000). Severe misfit of the 1PL by
   their generalized residuals. Equating differences exceeded the 0.5
   criterion at some raw scores. Proficiency classification disagreed for
   0 percent in two subjects and 2.4 percent in the third.
3. Two basic skills tests (40 and 38 items, n 2,000 to 3,000). Equated
   scores with and without the five most misfitting items correlated .99.
   Pass or fail disagreement was 3.3 and 2.4 percent, against 2.8 and 1.8
   percent when five items were removed at random, so misfit added about
   0.5 points of disagreement.

**Conclusions.** The model did not fit any of the datasets, but the misfit
was not always practically significant. IRT true-score equating appeared
robust. Pattern scoring and adaptive testing were not examined, and the
authors expect misfit to matter more there. They also note that "how large
is large" remains a judgment, and call the assessment a never-ending task
akin to validation.

## 3. What to take and what to change

**Take.**

- Practical significance is defined by the use of the scores, and measured
  as agreement between the operational scores and those from a
  better-fitting alternative fitted to the same data. No simulated
  mechanism is needed, which removes the curve bias found for the
  dependence side (`item-impact-validation-design.md`, sections 12 and 14).
- The alternatives include item removal and category collapsing, both of
  which stay inside the Rasch family.
- A baseline for what the alternative costs regardless of misfit. Their
  random-removal comparison separated the misfit's contribution (0.5
  points) from the loss of five items.

**Change.**

- **Attribution.** They remove "the most misfitting items". With a
  cascade, the most misfitting items are not necessarily the cause. Our
  version compares several candidate actions and shows which one resolves
  the misfit of the others.
- **More general models.** A 3PL or GPCM changes the scoring rule and gives
  up the sufficiency of the raw score. For a Rasch analysis it can serve as
  a sensitivity reference ("would scores change if slopes were free?") but
  not as the alternative the analyst adopts.
- **Uncertainty.** Their comparisons are point estimates. FWB over persons
  gives intervals for every consequence measure, consistent with the infit
  intervals.
- **Baseline.** Random removal assumes the other items fit, which was the
  objection to it in the impact design note. A parametric bootstrap under
  the fitted Rasch model gives the same baseline without that assumption:
  simulate from the fitted model, apply the same action, and compute the
  same disagreement. The excess over that is the part due to misfit.
- **Uses.** Health surveys are rarely equated. Their uses are cut-offs
  (PHQ-9 at 10), group comparisons, change scores, and person estimates for
  further modelling. The consequence measures should follow these.

## 4. A candidate comparison

For a given dataset the analyst, or the function from the infit and Q3
patterns, lists candidate actions. For each action, the same set of
measures, with FWB intervals, side by side.

**Candidate actions.**

- Keep all items (reference).
- Remove item A, remove item B, and so on, one at a time.
- Merge an LD pair into one polytomous super-item (testlet). This keeps the
  content, which often matters more than a marginal gain in fit.
- Collapse categories of an item with disordered thresholds.
- Sensitivity reference only: a model with free slopes, or two dimensions
  for a content cluster.

**Measures.**

1. **Misfit after the action.** The infit MSQ intervals of the remaining
   items. Whether the deviations of other items resolve is the attribution
   evidence. In the LD example, removing one item of the pair should bring
   the apparently underfitting items back to about 1, while removing one of
   those items leaves the pair unchanged. This is the logic of purification
   in `RMdifTree`, applied to fit. A compact summary is the number of items
   whose interval lies outside a floor band, or the sum of departures beyond
   it.
2. **Person scores.** Change in WLE estimates relative to the full scale,
   on a common metric through the shared items: mean shift in SEM units,
   share of persons moving more than one SEM, KL, and rank correlation.
3. **Uses.** If a cut-off is given, classification agreement with the full
   scale. If a group variable is given, the change in the standardized group
   difference. These follow Sinharay and Haberman's principle that the use
   decides the measure.
4. **Precision and targeting.** Reliability (PSI or RMU) and the
   conditional SEM over the range of the cut-off or the population.
5. **Baseline.** For measures 2 and 3, the value expected from the action
   alone if the Rasch model held, from a parametric bootstrap under the
   fitted model. The report shows observed, expected and the excess.

**What it leaves to the analyst.** Content validity, the intended
construct, and the final choice. The function shows the evidence. It does
not rank the actions, because the data cannot say which items are the
construct.

## 5. The difficult cases

- **LD pair with apparently underfitting items.** Measure 1 separates the
  cause from the cascade. Measure 2 shows whether either removal, or the
  testlet, changes scores more than the action alone would.
- **Residual phase with many small deviations.** Candidate actions may not
  be needed. The question is answered by the sensitivity reference: if
  person scores under a model with free slopes, or with the suspected LD
  modelled, differ from the Rasch scores by less than a small fraction of
  an SEM, the misfit is not practically significant for scoring.
- **Content cluster of two or three items.** Measure 1 shows that removing
  any one member leaves the others misfitting. The two-dimensional
  sensitivity reference shows what the cluster does to scores. The
  decision between removing the cluster and accepting it is a content
  decision.

## 6. Relation to the earlier work

- Kept: infit MSQ as an n-stable effect size, FWB intervals for it, and the
  finding that per-item impact fails for coherent groups (which this
  approach handles by comparing actions, not items).
- Kept as an aid: the calibrated impact curves, as a guide to what MSQ .9 or
  1.1 means in score units under a given mechanism. Useful for teaching and
  for choosing a floor, not as a decision rule.
- Paused: the 31-hour validation run of the per-item function. The
  prototype and pilot remain in `dev/` for reference.

## 7. A simulation, if this direction holds

The data-generating structures would be the difficult cases, with the
truth known:

1. An LD pair plus apparently underfitting items, at several λ and scale
   lengths.
2. Several small misfits of mixed mechanisms (slopes, dependence, a weak
   second dimension), MSQ within about .9 to 1.15.
3. A content cluster of two or three items on a second dimension.
4. The null, for the baseline and the false alarm rate of measure 1.

Criteria: whether measure 1 points to the true cause, whether the excess
over baseline is near 0 when the model holds, and whether the action with
the best evidence is the one whose scores are closest to the true θ1. The
single gross misfit is kept only as a reference.

## 8. Questions to settle

1. **Scope of actions.** Removal, testlets and category collapsing inside
   the Rasch family, with more general models only as a sensitivity
   reference. Is that the right line?
2. **Who proposes the candidates.** The analyst, the function from the
   infit and Q3 patterns, or both?
3. **Default measures.** Person-score change and misfit after the action
   always, uses only when a cut-off or group variable is supplied?
4. **Sensitivity models.** Which more general models, and through which
   package? Free slopes need mirt (already used for the MML option) or a
   native GPCM. LD terms need an LLRM fit.
5. **Relation to existing functions.** Much of measure 1 exists
   (`RMitemInfit`, `RMlocdepQ3`, `RMdifTree`'s purification). The new
   function would mostly orchestrate refits and comparisons. Is that a
   function, a vignette workflow, or both?
