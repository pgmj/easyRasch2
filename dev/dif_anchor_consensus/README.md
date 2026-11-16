# Consensus forward anchor over random tie-breaks

Study 3 for `RMdifAnchor()` development. Development only, not package code.

Study 2 (`../dif_anchor_threshold`) showed that the forward (MTT) anchor depends on psychotools' random tie-break in the item ranking. The same data give different anchors and, with DIF present, different flags under different seeds. This study checks whether a consensus over many seeds removes that dependence without losing calibration or power.

## Method (`dif_consensus.R`)

The ranking is redone under S seeds (auxiliary tests computed once, as in Study 2). Per item:

- `freq`: share of seeds whose anchor contains the item
- `first`: share of seeds where the item is the first (untested) anchor item

The consensus anchor is the items with `freq >= .5` (if none, the most frequent). The reference item is the anchor item with the highest `first`, ties broken by `freq`, then item number. A second variant, `cons_all`, tests every item when the anchor has two or more items, so no reference item is needed. Testing an anchor item against the anchor mean is well defined (the covariance is singular, which the single-step adjustment handles).

## Design (`dif_consensus_sim.R`)

72 cells: Study 2 with impact 0 only (impact 0.5 changed the Study 2 error rates by at most .027). The data seeds are Study 2's, so the datasets are the same as the impact-0 datasets there.

Methods per replication: `single` (one seed), `cons50`, `cons100`, `cons200` (first 50, 100, 200 seeds of set A), `cons_all` (200 seeds, every item tested), `gini`. Anchored location and threshold tests with the single-step correction.

Reproducibility: two disjoint seed sets A and B of 200 seeds each. For S = 1, 50, 100, 200 (and `cons_all`, coded S = -200) it records whether A and B give the same anchor set, the same reference item, the same location flags and the same threshold flags.

## Checks done

- Container smoke test 2026-11-15: 3 replications per cell, 4 workers, no errors, 2 skipped for an empty category.
- The Gini rows reproduce Study 2 exactly (428 of 428 rows, same flags), so data generation and the Monte-Carlo seeds match.
- Smoke-test hint only (3 replications per cell, not a result): with 200 seeds the anchor set agreed between seed sets in about 78% of datasets and the reference item in about half. The reference rule may need replacing, which is why `cons_all` is included.

## Result (host run 2026-11-15, 1000 reps x 72 cells, 134 min)

Written up as Study 3 in `../dif_anchor_tests.qmd`. The consensus keeps the forward anchor's FWER and power (location FWER mean .051 for both). Flag disagreement between two seed sets: thresholds .037 (max .267) with one seed, .011 (max .056) with 200 seeds; most of the gain by 50 seeds. The rest comes from items with an inclusion frequency near .5 (44% of datasets have one in .4-.6; disagreement .041 there, .002 elsewhere), so more seeds do not remove it. The most-often-first reference item is unstable (same in .475 of datasets); testing all items (`cons_all`) avoids it at no cost. Contamination with three non-symmetric DSF items is unchanged (location FWER up to .088).

## Host run

From the easyRasch2 package root:

```
nohup Rscript dev/dif_anchor_consensus/dif_consensus_sim.R 1000 10 > dev/dif_anchor_consensus/run.log 2>&1 &
```

Arguments: replications per cell, cores, output folder, B (Monte-Carlo draws, default 5000), seeds per set (default 200). Estimated run time about 4 hours on 10 cores (about 1.9 core-seconds per replication in the container). With 500 replications per cell it is about 2 hours, with a Monte-Carlo SE of about .01 for a rate near .05. Output: `dif_consensus_sim_<date>.rds`, `dif_consensus_sim.rds`, `dif_consensus_sim_summary.csv`, `dif_consensus_sim_agree.csv`.
