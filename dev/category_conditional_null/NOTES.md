# Category-by-score conditional test (prototype)

Started 2026-11-16. Side quest to understand threshold disorder, not for the
psychometrics skill.

## Idea

Under the PCM, P(X_i = j | total score r) depends only on the item
parameters, so the item-category by score table has model-expected counts
without any theta. This is the table behind Glas's CML test R1c for
polytomous items (Glas 1991, Cito report 91-2, eq. 29-30, scaled deviates
d_rij per item, category and score level). Earlier versions are Glas (1988,
Psychometrika 53, 525-546, dichotomous) and Glas & Verhelst (1989,
Psychometrika 54, 635-659, MML).

What is added here:
- focused statistics per category: rest-score mean and SD within the
  category (observed minus model-expected), the category's count in the
  score tails, and Q_j, the sum of squared Glas deviates for that category;
- a parametric bootstrap that samples response patterns exactly from
  P(x | r) given the observed totals (no theta drawn), and refits the CML
  thresholds in every replicate.

## Files

- `R/functions.R`: conditional probabilities via elementary symmetric
  functions (convolution), exact sampler given totals, statistics, bootstrap,
  and data generators for conditions A-C.
- `R/smoke.R`: one dataset per condition (B = 200) plus a 20-dataset null
  check under A (B = 100). Results in `results/smoke*.rds`.
- `R/study_functions.R`: generators for the 14 study conditions (step-specific
  discriminations as in Glas & Verhelst 1989 sec. 7, swapped categories,
  "don't know" use on a wide or narrow middle, extreme response style),
  Glas item index for every item, category statistics for item 1, Null 1
  adjacent-category AUC, Wald z for threshold order.
- `R/run_study.R`: host run, resumable per cell in `results/cells/`.
  Smoke mode: `SMOKE=1 Rscript R/run_study.R` (writes `results/smoke_cells/`).

Moved from analyses/ to easyRasch2/dev/ on 2026-11-16 (user request).

## Host run

From this folder:

    LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 nohup Rscript R/run_study.R > results/run_study.log 2>&1 &

Design: n = 250, 500, 1000; k = 5, 10 items, 5 categories each; 14
conditions (2 null, 12 alternatives). Null cells 500 reps with B = 199,
alternatives 200 reps with B = 99. 10 cores. Null cells run first.
Timing in the container: one refit at n = 1000, k = 10 takes about 0.07 s,
so a full run should take roughly 3 hours on 10 cores (rough estimate).

Bootstrap replicates with an empty category in some item are dropped
(recorded as B_used). Smoke: 14 of 19 kept in one n = 250, k = 5 dataset,
48 of 49 at n = 250, k = 10. Check B_used in the results; at n = 250 this
conditions the bootstrap on all categories being observed, as the observed
data were.

## Conditions (focus item 1 of 6, 5 categories, n = 500)

Categories are scored 0-4 (0 = lowest) everywhere in this folder. So
"categories 1 and 2" are the second and third response options (on a 5-point
agreement scale: disagree and neutral), the pair Wetzel & Carstensen (2014)
swapped. The middle category is 2. Thresholds t1-t4 sit between categories
0/1, 1/2, 2/3 and 3/4.

- A: PCM with disordered thresholds (-1.5, 0.5, -0.5, 1.5). The null is true.
- B: categories 1 and 2 interchangeable (generated merged, split at random).
- C: "don't know" use. With probability pi the middle category is chosen
  whatever theta.

## Smoke results 2026-11-16 (single datasets, not a study)

- Sampler checked against brute-force simulation: matches.
- A: no flags. Null check over 20 datasets: rejection shares 0-.05.
- B: thresholds ordered but close (-0.97, -0.58). Mean, SD and Q for
  categories 1 and 2 flagged (p about .01). Rest-score means 8.84 vs 8.95.
- C (pi = .25): SD and Q of the middle category flagged (p .04, .025). The
  tail count was not.
- Threshold check without bootstrap (20 reps, n = 1000): "don't know" use
  WIDENS the middle category and squeezes or disorders the outer pairs
  (t1/t2, t3/t4). It never disordered t2/t3. With generating thresholds
  (-1, -.2, .2, 1) and pi = .25, P(t2 < t1) = .85 and P(t4 < t3) = .75.
  So "don't know" use does not by itself CAUSE middle-category disorder. That
  does not rule out the two co-occurring (user's point, 2026-11-16): if few
  respondents are truly neutral, the trait-based middle is narrow and rare,
  and "don't know" use adds theta-unrelated responses on top of it.
- User's scenario, one dataset each (B = 200): middle thresholds (0.3, -0.3),
  pi = 0, .05, .10. Estimated middle pair: (0.27, -0.05), (0.02, 0.19),
  (-0.17, 0.39). "Don't know" use REMOVED the disorder. Middle SD/Q flagged
  only at n = 1000, pi = .10 (p .01, .015), not at n = 500.

## Purpose (user, 2026-11-16)

Main purpose: evidence that the response scale is not working as intended,
for any category. How the middle category is used is a secondary purpose.
Glas & Verhelst (1989, sec. 7) illustrated sensitivity with category-specific
discrimination parameters a_ij in exp(a_ij(j*theta - eta_ij)), a22 = .5,
a23 = 2. That is a natural general alternative: each step should discriminate
like the integer scoring says.

## Open

- Which mechanism produces an underused middle category (midpoint avoidance,
  extreme response style, "don't know" coded elsewhere)?
- Proper simulation: more conditions, pi grid, n and number of items, Null 1
  statistics (adjacent-category AUC/slope) alongside, enough reps. Host run.
- Literature: published version of Glas's polytomous CML test, IRTree models
  for midpoint use.

## Host run 1 (2026-11-16): discarded, RNG bug

`other_items()` called `set.seed(20261116 + k)` inside `gen_data()`. That
reset the global RNG on every call, so after each worker's first replicate
every theta draw and every bootstrap was the same. A null cell of 500
replicates held 16 distinct datasets, and rejection rates came out as 0,
about .1, .25, .33 or 1. The run has no information about the test.

Fixed 2026-11-17: `other_items()` uses `withr::with_seed()`, which restores
the caller's RNG state. Checked that the state is unchanged after the call and
that 40 of 40 replicates under `mclapply` differ. The smoke mode runs from
start to end. `R/summarise.R` now warns if a cell repeats datasets. The
discarded cells, log and summary are in `results/cells_rngbug/`. The smoke
run before the study used `R/functions.R` only, not `other_items()`, so the
smoke results above are not affected.

Summary step: `Rscript R/summarise.R` collects the cells into
`results/summary.rds` (one row per replicate, plus a dated copy).

Why the smoke mode missed it: it ran 4 null replicates on 4 cores, one per
worker, and the bug only repeats data from a worker's second replicate on.
Smoke mode now runs 12 null and 8 alternative replicates on 4 cores, so every
worker runs at least two. `SMOKE=1 Rscript R/summarise.R` checks
`results/smoke_cells` for repeated datasets. Checked 2026-11-17 with the old
`other_items()` put back in: 4 replicates on 4 cores gave 4 distinct
datasets, 12 gave 5. The guard flags the archived host cell and passes the
new smoke cells.

## Host run 2 (2026-11-17): results

Run 11:37-13:53, 84 cells, 20,400 replicates, 16 skipped (empty category in
the data, all n = 250). `R/summarise.R` now uses `dplyr::bind_rows` (skipped
rows have fewer columns). No repeated datasets.

Null (500 reps per cell, B = 199): every category statistic and the item
index reject at .02-.066 per cell, mean .035-.048 (MC SE .01). The PCM with
disordered thresholds (A1) is not flagged. The Wald test of threshold order
flags A1 in .69-1.00 of datasets and A0 in none.

Bootstrap retention: B_used below 90% of B in 7% (k = 5) and 11% (k = 10) of
datasets at n = 250, almost never at n >= 500. Eight datasets at n = 250,
k = 5 kept no bootstrap replicates, which gives p = 1. These should be
reported as not evaluable, not as non-rejections.

Multiplicity: rejecting when any of the 20 category statistics has p < .05
gives .44-.51 under the null. Holm over the 20 can never reject because the
smallest possible p is 1/(B + 1) (.005 or .01) and 20 times that is above
.05. The single calibrated test per item is the Glas item index. A joint test
would need the bootstrap matrices (not saved) for a max-T or min-p
adjustment, or a larger B.

Power (unadjusted, focus item, rejection at .05):
- SWAP_1_2: mean_1 and mean_2 1.00 at every n. Item index .58-1.00. Wald
  disorder 0: thresholds stay ordered (mean -1.71, -0.07, 0.20, 1.39 at
  n = 1000, k = 10).
- S_inner_a0 (categories 1 and 2 interchangeable): mean_1 .83-1.00, Q_1
  .44-1.00, index .25-1.00. Disorder 0. Adjacent AUC for 1/2 rejects at .03-.04,
  so Null 1 (interchangeable) is retained, as it should be.
- S_outer_a0: Q_0 .44-.97 and mean_1 .46-.93. Index weak (.04-.32).
- S_inner_a05 and S_inner_a2: mean statistics .23-.95, index and Q weak.
- DK: sd_2 is the main signal, then Q_2 and tail_2. pi = .25: sd_2 .25-.92
  (wide), .34-.98 (narrow). pi = .05: .05-.28. Wald disorder flags
  DK_narrow_05 in .10-.25 and DK_narrow_10/25 in none: "don't know" use
  removes the disorder that a rare true-neutral middle produces.
- ERS_30 (style on every item): nothing detects it (.03-.10). The focus-item
  test cannot see a person-level style shared by all items.
- Non-focus items: Holm-adjusted index flags no other item in any
  alternative cell.
