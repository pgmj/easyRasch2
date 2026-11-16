# Threshold-level DIF tests (RMdifAnchor development)

Status 2026-11-14: scripts written and smoke-tested in the container (up to 10 replications per cell, 4 workers, no failures), including the unequal-size and tie-break additions. The full run is still to be done on the host.

Status 2026-11-15: full host run done (1000 reps x 144 cells, 42 min on 10 cores). Results are written up in `../dif_anchor_tests.qmd` (Study 2). Main points: the spacing test holds the FWER with single-step or Holm (BH does not under DSF); forward is mostly calibrated but liberal with three non-symmetric DSF items (up to .097) and with three DIF items at 150 per group; Gini is conservative; forward anchor sets and flags depend on the random tie-break (some seed changes a threshold flag in up to 61% of datasets), Gini does not. Up to 47% of replications were skipped for an empty category in the smallest unequal cells.

This study is written to become a chapter in the planned "Differential item functioning" part of `rasch_book`, together with the item-level anchor study in `../dif_lr_wald/`.

## Why

`RMdifLR()` will be reduced to Andersen's LR test (decided 2026-11-14). Item- and threshold-level analysis moves to a new `RMdifAnchor()`, with forward (MTT) and Gini anchors. The old threshold table compared per-group thresholds on the centred scale, which has the artificial-DIF problem shown in `../dif_lr_wald/`. This study checks two replacements before they are implemented.

## Methods (`dif_thresh.R`)

Partial credit model, item i, threshold k: delta_ik = beta_i + tau_ik, where beta_i is the mean of the item's thresholds. Fixing a group's scale adds one constant to every delta in that group. Thresholds and their covariance come from `psychotools::threshpar(type = "mode", vcov = TRUE)`, and everything below is a contrast that does not depend on how each group is normalised.

- **Spacing test (DSF, no anchor).** The group differences in the adjacent distances delta_i,k+1 - delta_ik, tested jointly per item with a Wald chi-square on m_i - 1 df. The scale shift cancels, so there is no anchor and no artificial DIF. Uniform DIF leaves the spacing unchanged and must not trigger it. Corrections across items: Holm, BH, and a Monte-Carlo single-step min-p from the joint normal of all contrasts.
- **Anchored threshold differences.** Each group is shifted so that the mean location of the anchor items is zero. delta_ik(2) - delta_ik(1) is then tested per threshold. This is the location difference plus the spacing difference. Corrections across all thresholds: Holm and single-step max-|z|.
- **Anchored item locations.** The mean of an item's anchored threshold differences, with the first anchor item untested. This reproduces `psychotools::anchortest()` (z identical to 1e-14 for forward and Gini, 10 x 4 and 20 x 3, anchors identical).
- **Anchors** are selected at item level by `psychotools::anchor()`: centred (all items, the old RMdifLR scale, included as the reference that fails), forward (MTT) and Gini (constant, length 1, which is Strobl et al. 2021's anchor point).

The single-step adjustment is computed by Monte Carlo (B draws from the joint normal, |z| compared with the simulated max-|z|). multcomp's `pmvnorm` integration gives the same adjustment (largest difference 0.003 with B = 20,000 on a 30-threshold example) but took 1-2 s per call with 30-40 correlated tests, which would have made the run about 11 hours. The covariance of the anchored differences is singular (the anchor items' location differences sum to zero). Both methods handle that.

## Design (`dif_thresh_sim.R`, 144 cells)

- mean group size n = 150, 300, 600
- group sizes equal (n, n) or unequal (1.5 n, 0.5 n). Group 2 is the smaller focal group with the impact and the DIF. This extends the item-level study in `../dif_lr_wald/` (equal sizes only) through the location family here, for polytomous items.
- 10 items x 4 categories (3 thresholds), or 20 items x 3 categories (2 thresholds)
- impact 0 or 0.5 (group 2 latent mean, SD 1)
- scenarios, as shifts of group 2's thresholds per DIF item:
  - `null`
  - `unif1`, `unif3`: every threshold +0.6 in one or three items (location DIF only)
  - `sym1`: lowest threshold -0.6 and highest +0.6 in one item (symmetric DSF, location unchanged)
  - `nonsym1`, `nonsym3`: highest threshold +0.6 in one or three items (non-symmetric DSF, location +0.2 with 3 thresholds, +0.3 with 2)

Null sets: for spacing, items whose thresholds all shift by the same amount; for thresholds, thresholds with shift 0; for locations, items with mean shift 0. Replications with an unused category in either group are skipped, as RMdifLR() requires all categories. Each replication is seeded by cell and replication.

**Tie-break check.** psychotools breaks ties in the anchor ranking at random, and only the ranking step draws random numbers. In the first 250 replications of each cell, the forward and Gini rankings are redone under 20 seeds (auxiliary tests computed once). For each seed the anchors and single-step flags are compared with the main run. The Monte-Carlo adjustments use a fixed seed per replication, so flags change only through the anchor. Reported: unique anchor sets, the share of seeds with another set or another untested first item, the share of replications where some seed changes a flag, and the share where the number of true DIF units found depends on the seed. A first look on five 10 x 4 datasets: forward gave 8 to 20 distinct anchor orderings in 20 seeds, Gini always one.

Output per replication, family (spacing, thresholds, locations), anchor and correction: counts of false and true flags, null-set z moments, anchor length, and whether the anchor holds a location-DIF item (`anchor_locdif`) or any shifted item (`anchor_anyshift`, which includes symmetric DSF items that are valid anchors).

## Host run

From the package root (smoke-test timing suggests about 70 minutes on 10 cores, rds about 30 MB):

```
nohup Rscript dev/dif_anchor_threshold/dif_thresh_sim.R 1000 10 > dev/dif_anchor_threshold/run.log 2>&1 &
```

Arguments: replications per cell, cores, output folder, B (Monte-Carlo draws, default 5000), tie_reps (default 250). The run writes `dif_thresh_sim_<date>.rds`, `dif_thresh_sim.rds`, `dif_thresh_sim_summary.csv` and `dif_thresh_sim_ties.csv` here, and prints the single-step and tie-break summaries at the end of `run.log`.

## Still open

- Kopf et al. (2015a/b) labels: not an error. Strobl et al. (2021) and Henninger et al. (2025) label the EPM anchor-selection paper 2015a (APA ignores a leading "A"), while the book's apa.csl orders by the full title and labels the APM framework paper 2015a, as our docs do. Each document is consistent within itself. A real error: `R/dif_lr.R` cites Kopf et al. (2015a) for random tie-breaking, but neither paper nor the psychotools help says so. The source is psychotools' `anchorselect()` (it calls `sample()` on tied ranks). Fix when the code moves to `RMdifAnchor()`.
- Not covered: more than two groups, DIF sizes other than 0.6, unequal group sizes for dichotomous items, non-uniform DIF in the discrimination sense (outside the PCM).
