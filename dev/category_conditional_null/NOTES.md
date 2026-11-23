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

Report: `dev/category_conditional_null.qmd` (Typst PDF, 21 pages), written
2026-11-18 from `results/summary.rds`. Rates leave out the 8 datasets with
no bootstrap replicate kept. Numbers quoted in the text were rechecked on
the evaluable datasets.

## Diagnostic plots (prototype, 2026-11-18)

`R/plot_functions.R` + `R/plot_examples.R` (one dataset per condition,
n = 500, k = 10, B = 199, 1 min in the container). PNGs in
`results/plots/`, data behind them in `results/plot_examples.rds`.

A bar plot of observed vs expected category counts over the whole sample
shows nothing: the category counts are sufficient for the item's CML
parameters, so marginal expected = observed. Every plot below splits by total
score.

1. `1_cond_prob`: observed share per category at each total score vs
   P(X = j | r), pointwise binomial 95% band given n_r. Swap is visible in
   categories 1 and 2, but per-score points are noisy.
2. `2a_group_bars` / `2b_group_resid`: five score groups of about equal size,
   observed count (bar) vs expected with a 95% interval, or the same as
   standardised residuals. 2b is the clearest single picture: swap gives a
   sign reversal across groups (2 over-used low, 1 over-used high).
3. `3_deviate_map`: Glas deviates per score level and category, diverging
   blue-red. Shows the swap band pattern; noisy at single score levels.
4. `4_adjacent`: share in the upper of two adjacent categories vs its CML
   curve. Interchangeable or swapped 1/2 gives a flat or falling curve.
5. `5_focused`: observed minus expected rest-score mean and SD per category
   against the 2.5-97.5% bootstrap null range. Maps one-to-one to the tests.

None shows step discrimination 2 clearly, matching the low power there.
All bands are pointwise and (1-4) ignore item-parameter uncertainty.

### PAV calibration of the cumulative events (2026-11-18)

User asked whether the Säilynoja et al. (2025) PAV-adjusted calibration plots
apply here. `pav_table()` / `plot_pav()` in `R/plot_functions.R`, demo
`R/plot_pav_examples.R` (same five datasets and seeds, 7 s) ->
`results/plots/6_pav.png`. The forecast for event X_1 >= k is the CML
conditional P(X_1 >= k | r), so no theta enters. The band comes from responses
sampled exactly given totals at the fitted parameters (no refit), PAV in the
same ordering, pointwise 95%.

- Swapped and interchangeable 1/2: the X >= 2 curve is flat around .5 (62%
  and 40% of points outside the band). Other events stay on the diagonal, so
  the plot points to the step.
- Step 2 discrimination 2: S-shape in X >= 2 (13% outside). This is the first
  plot here that shows this condition at all.
- "Don't know" 25%: mild, X >= 3 11%, X >= 2 8%.
- Null (disordered PCM): 0-5% per event.
- Curves cannot be off the diagonal on average: CML fits the category counts,
  so the mean deviation per event is 0. Only shape (slope) shows (same as the
  easyRaschBayes check of 2026-11-06, dev/pav-calibration-plan.md there).

Monotonicity check (exact ESF, no simulation): P(X_i >= k | r) is
nondecreasing in r when all OTHER items have ordered thresholds (their ESF is
then log-concave, so the conditional family is TP2). One disordered other item
broke it in 30 of 73 random cases. The focus item's own disorder does not.
PAV orders by p, not by r, so the plot stays valid either way.

Side finding: `.gamma()` uses `stats::convolve` (FFT). Max error in P(X = j |
r) about 5e-6, but differences of cumulative probabilities at extreme scores
were off by up to 1e-3. Direct polynomial multiplication is exact and cheap
at these sizes. Not changed.

### 2026-11-18 later: report, exact ESF, calibration test

- `.gamma()` in `R/functions.R` now multiplies polynomials directly (FFT
  version kept nowhere; difference to brute-force enumeration 4e-16 vs 2e-10,
  up to 0.0011 in P(X = j | r) at extreme scores for the k = 10 design). In 60
  datasets the observed Q_j moved by at most 0.007, other statistics by
  < 1e-4, so run 2 was not repeated. The report says so in "Run".
  Example plots regenerated with it (old rds in `results/old_fft_gamma/`),
  numbers unchanged.
- Report: new section "Diagnostic plots for single datasets" (2b, 4, 5,
  PAV) before the Discussion, plus Dimitriadis 2021 and Säilynoja 2025 refs.
  Rendered, colour + greyscale checked (residual bars now red vs dark navy so
  the sign legend survives greyscale; dataviz validator not run, no node in
  the container).
- Calibration test: `R/calibration_functions.R` (`.pav()` moved here,
  `mcb_matrix()`, `analyse_calibration()`). Statistic per item and event:
  CORP Brier miscalibration MCB = mean (p - y)^2 - mean (cep - y)^2 with p the
  CML P(X_i >= k | r). Reference: exact conditional bootstrap with refit.
  Raw p per event, single-step min-p over the 4 events of an item, item p =
  smallest adjusted p, Holm over items. PAV checked against isoreg (5e-17).
- Host run: `R/run_calibration_study.R` (same 84 cells, reps and B as run 2,
  new datasets, cells in `results/cal_cells/`, about 2.3 h on 10 cores), then
  `Rscript R/summarise_calibration.R` -> `results/calibration_summary.rds`.
  Smoke run done (120 reps, 1 skipped, no repeated datasets; SWAP at the
  minimum p with B = 19, null mean p .5-.6).

## Package draft: RMitemCatFitPlot() (2026-11-18)

Draft in `package_draft/item_cat_fit_plot.R`, not yet in `R/`. One function,
`type = "residuals" | "adjacent" | "focused" | "calibration"` (plots 2b, 4, 5,
6 of the prototype), items in facet rows, `output = "ggplot" | "dataframe"`.

- Reuses package internals: `.fit_cml_thresholds()`, `.esf_convolve()`,
  `.maybe_parallel_lapply()`, `.n_caption()`, `.attrition_clause()`,
  `er2_caption()`, `er2_plot_caption()`, `er2_axis_margins()`.
- New internals: `.cond_cat_probs()` (P(X_i = j | r), all items),
  `.sim_given_totals()` (vectorised over respondents, unlike the
  one-pattern-at-a-time `.sim_conditional()`), `.pav()` (stack version),
  table and ggplot helpers per type.
- Complete cases only. Extreme total scores excluded. Caption via `.n_caption()`.
- Focused: refit in every simulated dataset, datasets with an empty category
  skipped and reported via `.attrition_clause()`. Calibration: band at the
  fitted parameters without refit, from the same simulation scheme.
- Checked in a scratch copy of the package with load_all: conditional
  probabilities match `cond_cat_prob()` to 4e-16, `.pav()` matches the dev
  `.pav()`, sampler keeps totals and reproduces P(X | r), focused observed
  values equal `item_stats()`, same seed gives identical output, Rd builds.
- Timing at n = 500, 10 items, 5 categories, 100 iterations: residuals and
  adjacent 0.2 s, calibration 1.6 s, focused 7.7 s (sequential).
- Not done: tests, NAMESPACE/pkgdown/NEWS, palette validator (no node).
- Renamed from RMitemCategoryPlot() to RMitemCatFitPlot() (user, 2026-11-18).
  If the calibration test is adopted, it becomes RMitemCatFit() (table of
  item p-values), and RMitemCatFitPlot() gets a `fit` argument: reuse the
  test's simulated datasets for the bands (refitted) and put each item's
  adjusted p in its row strip. Whether RMitemCatFit() also reports the
  focused statistics is decided after the calibration run.

## Calibration test run, results (2026-11-18 evening)

Host run finished. `summarise_calibration.R` now drops and warns about NULL
results instead of failing.

- **SWAP_1_2 n1000 k10 is empty.** All 200 results NULL, saved 19:44. The log
  was restarted, so the run was most likely interrupted during that cell and
  mclapply returned NULL for killed workers. The cell needs a rerun.
- Other losses: 16 reps skipped (empty category, n = 250). 5 reps at n250 k5
  had B_used = 0 (every refit failed).
- **Raw per-event p calibrated:** P(p <= .05) = .049-.050 for all four events
  under the null (6000 item-reps).
- **Item-level min-p p is too small in the lower tail:** P(p <= .01) = .016
  (k5) / .019 (k10), P(p <= .05) = .057 / .056. Run 2's focused index: .008,
  .049. Cause: bootstrap p's on denominator Bu, observed on Bu + 1, so with
  `minp <= p_raw` replicates tied on the same count are not counted. Toy
  check with 4 exchangeable correlated normals (scratchpad minp_check.R)
  reproduces it: .0168 at .01, .0564 at .05. Fix: rank observed and
  replicates in one pool of Bu + 1 and take p = share of the pool with min-p
  <= the observed min-p. Valid: at B = 199 it gives .0013 / .045 at .01 / .05
  because the grid is coarse, at B = 999 .0088 / .0497.
- **Holm columns are uninformative.** At B = 199 the smallest p is .005, and
  times 10 items that is .05, so nothing can be rejected at k = 10. At B = 99
  (alternatives) not even at k = 5. any_other_holm at k5 null = .084 is the
  min-p bug (4 x ~.02).
- **Power, item 1 at .05 per item, against run 2's focused index** (independent
  datasets, 200 reps per cell, calibration test's ~.057 null rate noted):
  - outer categories interchangeable: .42-1.00 vs .04-.32
  - step discrimination 2: n1000 k10 .555 vs .025, n1000 k5 .345 vs .085
  - inner categories interchangeable k10: n250 .72 vs .25, n500 .98 vs .55
  - swap n250 k10: .99 vs .58
  - don't-know 25% k10 n1000: narrow .505 vs .285, wide .35 vs .145
  - focused ahead only in don't-know 5-10% narrow, by .03-.05, close to null
    rates for both. ERS and don't-know 5-10% wide: neither detects.
- The event that rejects points at the category: outer -> X >= 1, inner and
  swap -> X >= 2, don't-know narrow -> X >= 2 and X >= 3.

### Fix and rerun prepared (2026-11-18, late)

- `analyse_calibration()` now uses `.minp_pooled()`: observed and replicates
  ranked in one pool, item p = share of the pool with min-p <= the observed.
  Matches the toy `pooled_p()` exactly (200 random matrices with ties, and
  one event). `p_adj` is gone. `boot_mcb` (k x 4 x B_used, 7 significant
  digits) is stored per replication. Old version in scratchpad
  calibration_functions_v1.R.
- First-run results moved to `results/old_minp_bug/` (cells, summaries, log).
- `summarise_calibration.R`: p_adj columns dropped, p_item_all_1..k added for
  null calibration over all items.
- Smoke run in the container: 120 reps, 1 skipped, no repeated datasets.
  No item rejects at B = 19 (smallest item p about 4/20), as expected.
- Same seeds and cores as the first run, so the datasets are the same
  (data_id can be compared). Expected time about 2.3 h. Disk: about 0.6 GB
  for the cells, twice that with the dated copies.
- Run script reads CORES=<n> (default 10, smoke 4); a different count gives different datasets than the first run.

### TODO: clean up results/ (user, 2026-11-18)
After the calibration rerun is analysed: decide what to keep of
results/old_minp_bug/, results/old_fft_gamma/, smoke summaries and dated
duplicates, and which rds files the qmd and a rasch_book chapter need.

## Calibration test rerun, results (2026-11-19)

Rerun on 4 cores (so new datasets), 22:25-04:28, all 84 cells complete,
19 reps skipped (empty category), 5 with B_used = 0 (n250 k5). cal_cells 1.2 GB.

- Null, item p over all items (k5: 14985 item-reps, k10: 29960):
  P(p <= .05) = .039 / .041, P(p <= .10) = .095 / .094, P(p <= .025) = .019,
  P(p <= .01) = 0 (smallest item p seen .015). Same at every n. Valid, and
  below nominal at small alpha because of the grid at B = 199, as the toy
  check predicted (.045 at .05, .001 at .01).
- Raw event p: .045-.050 at .05.
- Holm over items cannot reject at B = 199 (needs p <= .01 at k5), so
  FWER = 0 is uninformative. The toy check at B = 999 gave .0088 at .01 and
  .0497 at .05, so the package needs iterations of about 1000.
- Power at .05 per item barely moved from the first run (differences within
  +-.07, 200 reps), and the comparison with run 2's focused index stands:
  outer interchangeable .40-.99 vs .04-.32, step discrimination 2 n1000 k10
  .53 vs .03, inner interchangeable k10 n250 .73 vs .25, swap n250 k10 .97 vs
  .58, don't-know 25% narrow n1000 k10 .48 vs .29. Focused ahead only in
  don't-know 5-10% narrow, by .03-.07, where both are near their null rates.
  SWAP n1000 k10 (missing in the first run): 1.00.

### Write-up (2026-11-19)
Section "A calibration test from the PAV curves" added to
category_conditional_null.qmd (after the diagnostic plots): statistic,
p-values, run, null figure (share of item p <= alpha), power figure against
the item index, per-event figure at n1000 k10, appendix table per cell.
All numbers in the text are inline R from calibration_summary.rds. Discussion
gets one paragraph and an extra sentence in the multiplicity bullet. Toy
check saved as R/minp_toy_check.R. Reference added: Ge, Dudoit & Speed
(2003), checked against Crossref.

## Package draft: RMitemCatFit() (2026-11-19)

`package_draft/item_cat_fit.R`: RMitemCatFit(data, iterations = 1000,
correction = fwer/fdr_bh/fdr_by/none, item_labels, parallel, n_cores, seed,
output = kable/dataframe/list) plus the internals shared with the plot
(.catfit_prepare, .catfit_simulate, .mcb_vector, .catfit_pvalues).
RMitemCatFitPlot() now uses the same preparation and takes `fit` (the list
output): no new simulations, the corrected item p under each item label,
check that fit and data match.

- Every p computed on the pool of observed + B simulated datasets (event p,
  item p by single-step min-p, fwer = .wy_stepdown() on -item p), so ties
  count alike for observed and simulated. MCB equals the dev mcb_matrix()
  (1e-17), item p equals .minp_pooled() exactly, fwer >= p, fwer = p for one
  item, same seed gives identical output, mismatched data stops.
- Timing (n = 500, 10 items, sequential): 200 iterations 11 s, 1000
  iterations 52 s. Plots with `fit` take 0.1-0.2 s.
- **Floor of the corrected p.** Item p cannot go below m/(B + 1) (m events,
  one simulated dataset holds rank 1 for each event), and the step-down
  floor is about m k/(B + 1). Swap example, 1000 iterations: item p .004
  (= 4/1001), fwer .040. So reaching padj < .05 needs B > 20 m k: 800 for 10
  items with 5 categories, 2400 for 30 items. Studentised max-T within the
  item (like .bootstrap_pvalues) would give a floor of 1/(B + 1). It can be
  evaluated from the stored boot_mcb of the rerun without new simulations.

## Min-p vs studentised max-T (2026-11-19)

`R/compare_combination.R` recomputes both combinations from the stored
boot_mcb of the rerun (no new simulation, 20 s on 4 cores) and writes
results/compare_combination.rds. Recomputed min-p equals the stored item p
(2 of 152,805 items differ by up to .01 from the 7-digit rounding, no
rejection changes).

- Null, item p over all items: min-p .039 / .040 at .05 (k5 / k10), 0 at
  .01. Max-T .044 / .046 at .05, .0049 / .0044 at .01, per cell .039-.051 at
  .05 and .004-.007 at .01. Smallest p: min-p .015, max-T .005 (= 1/200).
- Null family-wise rate (any item, WY step-down): min-p .003 / .007 (the
  floor), max-T .045 / .038, per cell .030-.074. The .074 cell is the ordered
  PCM at n1000 k5 (500 datasets, about 2.5 Monte Carlo SEs above .05). All
  other cells <= .048.
- Power per item at .05: the same for both (differences within noise).
- Power after the correction, B = 99: max-T far ahead. Swap .85-1.00 vs
  .04-.19, outer interchangeable .23-.97 vs .01-.09, inner interchangeable
  .36-1.00 vs .03-.09, step discrimination 2 n1000 k10 .24 vs .02.
- Other items under alternatives, any flagged after the correction: max-T up
  to .075 per cell (mean .038).
- Decision: the draft RMitemCatFit() now uses max-T (old version in
  scratchpad item_cat_fit_minp.R). Draft p and fwer equal the comparison
  script exactly on 20 stored datasets. Swap example, 200 iterations:
  item p = padj = .005.
- Open: confirm the family-wise rate at a larger B, especially the ordered
  PCM at n1000 k5. The qmd still describes the min-p item p.

## Large-B null check prepared + qmd switched to max-T (2026-11-19)

- `R/run_calibration_B999.R` (host): [revised same day: n = 1000 only by default, both conditions, about 3.1 h, NS env var adds 250/500] six null cells at k = 5 (A0 ordered,
  A1 disordered, n = 250/500/1000), B = 999, 1000 new datasets per cell,
  about 7 h on 10 cores. Env vars CORES, REPS, B, SMOKE. Stores per dataset
  the item statistics (N x 5), min-p and max-T item p and the max-T
  step-down p (not the full MCB array, about 40 MB per cell). Output in
  results/cal_b999_cells/. Summary: `R/summarise_calibration_B999.R` ->
  results/calibration_b999_summary.rds. Smoke-tested in the container.
- New helpers in R/calibration_functions.R: `.stepdown_maxT()`,
  `.maxT_items()`, `analyse_calibration_maxT()`. `.maxT_items()` gives the
  same p and padj as compare_combination.R on a stored dataset.
- Convention: rejection is p < alpha (strict), as in easyRasch2. With
  p <= alpha the max-T item p of A0 n1000 k5 is .016 at .01 and .057 at
  .05 (from ties at the 1/(B+1) grid), which the B = 999 run should settle
  too. FWER of that cell: .074 strict, .078 with <=.
- qmd: calibration section now reports the max-T item p and step-down
  (from results/compare_combination.rds), min-p kept as the first version
  with its numbers. New appendix tables: null rates per cell (both
  combinations, FWER) and focus-item power after the step-down. Discussion
  multiplicity bullet updated. Backup of the previous qmd in the session
  scratchpad (category_conditional_null_pre_maxT.qmd). Renders, 38 pages.
- Corrected claim: after the step-down max-T keeps most of the power only
  for swapped and interchangeable categories, about a third for the weaker
  departures (step 2 discrimination 2 n1000 k10: .57 -> .24).

## B = 999 null check, results (2026-11-19)

Host run, 10 cores, 40 min per cell (estimate was 94, so about 3.1 h was
pessimistic: 81 min total). 1000 datasets per cell, none skipped, all data_id
unique, min B_used 996. p recomputed from stored Tm equals the stored p.

| Cell (n1000 k5) | item p < .05 | item p < .01 | FWER (padj < .05) | MC SE |
|---|---|---|---|---|
| A0 ordered | .052 | .010 | .052 | .007 |
| A1 disordered | .048 | .008 | .045 | .007 |

With <= instead of <: A0 .053 / .011 / .053, A1 .050 / .009 / .045.
Min-p item p at B = 999: .050 / .048 at .05 (calibrated once the grid is fine).

Same datasets, first 199 replicates only (studentised with the full-pool
moments, so an approximation): FWER .043 (A0) and .045 (A1), item p < .01
.005. So the .074 at B = 199 in the rerun was chance in those 500 datasets,
not an effect of B. At B = 199 the item p is slightly below nominal at .01
(.005) from the coarse grid.

Conclusion: max-T with the WY step-down holds the family-wise rate at the
n1000 k5 null cells. No change to RMitemCatFit needed.
- qmd updated with these results (inline from results/calibration_b999_summary.rds,
  whose `strict` table uses p < .05 and holds the B = 199 subsample), and the
  Discussion multiplicity bullet. summarise_calibration_B999.R now writes `strict`.

## Moved into the package (2026-11-19)

R/item_cat_fit.R and R/item_cat_fit_plot.R are now the source of truth;
package_draft/ is the pre-move copy. Changes on the move:
- `.catfit_prepare()` runs `.check_complete_case_categories()` (classed
  category error, highest category only in incomplete rows) with the calling
  function named in the message.
- Score-group labels built without tapply names (en dash as –, no
  re-encoding warning in a C locale).
- Roxygen: simulation numbers now say 199 simulated datasets and add the
  B = 999 check.
- NAMESPACE, man/, _pkgdown.yml (Item Fit, after RMitemCatProb), NEWS.md (two
  bullets), README.md (category section).
- tests/testthat/test-item_cat_fit.R: exactness of .cond_cat_probs (brute
  force), .sim_given_totals (totals kept, frequencies), .pav vs isoreg,
  .catfit_pvalues floors and step-down vs an independent implementation,
  argument checks, category error, output structure, seed reproducibility,
  swapped-category detection, complete cases, low-iteration notice,
  parallel = sequential, plots for every type, residual counts summing to
  zero (CML sufficiency), fit reuse equals the plot's own simulation, strips
  with p, fit from other data. All pass (126 expectations).
- The parallel test needs the dev build installed (mirai workers load the
  installed package). In the container 1.4.0 is installed, so it was run with
  the dev build in a scratch library via R_LIBS.
- Found while testing: the adjacent-category binomial band need not contain
  the expected share when it is near 0 or 1 (discrete quantiles). Test
  changed, function unchanged.

## Book chapter drafted (2026-11-19)

rasch_book/chapters/category-fit.qmd (Item fit part), scripts in
rasch_book/scripts/category-fit/ (paths rewritten, dev-history comments
trimmed), results in rasch_book/results/category-fit/ (summary,
calibration_summary, compare_combination, calibration_b999_summary). Book
.gitignore now covers results/**/*cells*/, *smoke*, *.log. Five refs added
to references.bib (Crossref/arXiv checked). Audit of the dev qmd before
promotion: every numeric claim in the first study matched the summary.
Open: study 1 results were computed with the FFT ESF, the scripts now use
direct multiplication (one sentence in the chapter's Reproduction section).
The chapter example (n = 1000, 8 items, Q1 swapped, Q2 categories 0/1
interchangeable, seed 2026) runs RMitemCatFit() at render, about 1 min.

## Code review fixes (2026-11-19)

`/code-review high` on RMitemCatFit/RMitemCatFitPlot, log in dev/code-review/item_cat_fit.md. All seven findings fixed (CF1-CF7). The main one: a category used only by extreme scorers now stops in `.catfit_prepare()` with a classed error naming the item, before any simulation. Before the fix, RMitemCatFit gave a generic refit failure and the adjacent/residual plots were drawn on a diverged threshold with no warning. The chapter's results are unchanged.

## Vignette (2026-11-19)

New subsection "Category fit given the total score" in section 3 of vignettes/easyRasch2.Rmd.orig: kable table, calibration plot with `fit`, adjacent plot. PHQ-9, iterations 200, seed 9. Items 3 and 9 adjusted p < .01 (also at 1000 iterations). Andrich 2013 and Adams et al. 2012 added to vignettes/references.bib. The .Rmd and figures are not regenerated yet (precompile on the host).

## Rename (2026-11-19)

`RMitemCatFitPlot(type = "focused")` renamed to `type = "restscore"`, since the plot shows the rest-score mean and SD per category and "focused" was a dev label. Internal `.catplot_focused*` renamed to `.catplot_restscore*`. Tests, Rd, vignette and the book chapter (label fig-cfit-example-restscore) updated. Earlier sections of this file keep the old name.

## Plot layout and p marking (2026-11-20)

Calibration plot uses coord_cartesian instead of coord_equal, so panels fill the width. With `fit`, strips show the corrected p in bold when < .05 (ggtext, labels escaped as HTML entities), or with " *" without ggtext. The caption says which. Test added. Vignette calibration fig.height 13 -> 11.
- Rest-score plot: outside points now filled vermillion #D55E00 (as RMdifAnchor flags), inside open black; bars grey60 -> grey85 for contrast. Shape still carries the distinction (CVD-safe, checked with a deutan simulation).
- 2026-11-20: help pages of both functions reread and revised for clarity (terms defined, Data paragraph, p formula, FDR note, book link, Glas 1988 now cited, refs alphabetical).

## 2026-11-20: null check of the residual plot's ±1.96 lines

Question (user): are the standardised residuals in `RMitemCatFitPlot(type = "residuals")` N(0, 1) under the PCM? Not checked in any earlier study: the main study and the calibration study judged their statistics against the conditional bootstrap only.

Script `R/residual_null_check.R` (300 reps per cell, container, 0.9 min). PCM data, theta ~ N(0, 1), n = 200/500/1000, k = 5/10, m = 2/4 thresholds, score groups 3/5/10. z computed as the plot does (CML estimates, binomial variance) and with the true thresholds. Datasets with an empty category skipped (7 of 3600). Summary `results/residual_null_summary.rds`; cells in `results/residual_null_cells/` (gitignored).

- With the true thresholds, z is fine: SD 0.98-1.01, 4.6-5.2% beyond ±1.96. So the binomial variance is right and discreteness does not distort the overall rate.
- With the estimates (what the plot shows), z is under-dispersed. Default 5 groups: SD 0.86-0.90, 2.4-2.9% beyond. 3 groups: SD 0.77-0.80, 1.4-2.0%. 10 groups: SD 0.91-0.95, 3.2-3.8%. Close to sqrt((G-1)/G) (0.82, 0.89, 0.95), the loss from each category's residuals summing to zero over the G groups at the CML estimates. n and k barely matter.
- By expected count (5 groups): E > 20 SD ~0.82, 1.7-1.9% beyond; E < 5 4.1% beyond but nearly all positive (upper ~4.0%, lower ~0.1%). Sparse cells can only exceed upward.
- So the ±1.96 lines are conservative, by an amount that depends on score_groups and cell size, and asymmetric for sparse cells. Options put to the user: bootstrap bands from the simulations RMitemCatFit already stores, or keep the lines and say this in the help text. Not decided.
- Implemented same day (user's choice): `type = "residuals"` draws no ±1.96 lines. With `fit`, grey boxes give each cell's pointwise range from RMitemCatFit's simulations. `.catfit_simulate()` now also returns `seeds` and `thr` (the re-estimated thresholds) of each refitted dataset, and `.catplot_residual_band()` regenerates the datasets for any `score_groups` (sampling is given the observed totals, so the groups are fixed) without refitting. It restores the caller's .Random.seed. `.catplot_residuals()` rewritten on `.residual_sums()` (tabulate + rowsum), output checked identical to the old aggregate() version on pcmdat2, phq9 and random 3-category data at 2, 5, 10 groups. Dataframe output gains lower, upper, iterations_used (NA without fit). Cost: 0.18 s for 200 iterations on pcmdat2.
- Band null check `R/residual_band_check.R` (100 datasets, n = 400, k = 8, 4 categories, 200 iterations): 5.25% outside (2.50% below, 2.75% above) vs 2.61% beyond ±1.96. E < 5 cells: 1.3% below / 2.5% above vs 0.04% / 4.2%. Results `results/residual_band_check.rds`.
- Dev qmd section "Null distribution of the residuals in the residual plot" (before the calibration test) with fig-resid-null, tables tbl-resid-null and tbl-resid-band in the appendix. Rendered.

## 2026-11-20: rest-score plot and whole-item discrimination

User question: can the SD panel be read like discrimination (higher SD than expected = lower discrimination, underfit)? Checked with `R/restscore_discrimination_check.R` (container, 38 min): 10 items, 5 categories, focus item GPCM slope a = 0.5/1/2 (thresholds centred), n = 300/1000, 100 datasets per cell, RMitemCatFit 200 iterations, none skipped. Summary `results/restscore_discr_summary.rds`, cells gitignored.

- Direction confirmed in every category. a = 0.5: SD +0.80 to +0.99 in all categories, means compressed (cat 0/1 above, 3/4 below, middle ~0). a = 2: SD -0.9 to -1.4, means spread. a = 1: all ~0, per-category outside rates 1-9%.
- Correction of what I told the user beforehand: the SD is NOT the weaker signal. Mean shifts are largest in the extreme categories and zero in the middle one, SD shifts in every category. n = 1000, a = 2, middle category: SD below range 100%, mean outside 9%. n = 300, a = 2, inner categories: SD 56-94%, mean 1-48%.
- RMitemCatFit (calibration test) flags the focus item 36% (a = 0.5) / 8% (a = 2) at n = 300, 96% / 83% at n = 1000. Rest-score plot ranges are pointwise and uncorrected: under a = 1, 18-28% of datasets have some category of the focus item outside on a panel, so the two rates are not comparable (dev qmd wording fixed 2026-11-20 after the user asked which plot was meant).
- Dev qmd section "Rest-score plot and item discrimination" with fig-discr, tbl-discr and tbl-discr-test in the appendix. Rendered. Help page and book text proposed to the user, not yet written.

## 2026-11-20: docs, vignette and book

- Both residual checks rerun after the grouping change (`.catplot_residuals()` now takes a group factor, equal-count groups from the shared `.cicc_class_bins()`): 5 groups 2.5-2.8% beyond ±1.96 (SD 0.86-0.89), ranges 5.2% outside (2.4% below, 2.9% above). Dev qmd re-rendered.
- Help: "What a flag means" in RMitemCatFit (flag not specific to categories, check item fit first, slope numbers); rest-score SD reading in RMitemCatFitPlot.
- Vignette (.orig and knitted .Rmd, prose only): flag paragraph using items 3 (first step only, category) and 9 (all steps, item), rest-score reading sentence.
- rasch_book chapters/category-fit.qmd: summary + purpose sentences, residual example now counts residuals outside the simulated range (Q1 9 of 25, PCM items 0-2), rest-score reading, new sections "Reference ranges for the residual plot" and "Item-level misfit", discussion paragraph on what a flag means, two limits bullets, reproduction bullet, four appendix tables. Scripts copied to scripts/category-fit/ (installed easyRasch2), summaries to results/category-fit/. Rendered in the container with the dev package installed.
