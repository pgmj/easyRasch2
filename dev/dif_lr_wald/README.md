# Anchored Wald DIF tests (RMdifLR development)

Status 2026-11-13: full run done on the host (1000 replications per cell, 60.8 min). Results below.

- `dif_wald.R`: the prototype. `dif_wald_fit(data, group)` fits a CML partial credit model per group (two groups). `dif_wald(fit, anchors)` returns, per item and anchor, the group 2 minus group 1 location difference with SE, z, raw p, Holm and single-step (multcomp max-|z|) adjusted p.
  - `"centred"` is RMdifLR()'s current scale (all items) and tests every item.
  - `"constant"` (Gini) and `"forward"` (MTT) use `psychotools::anchortest()` and leave the first anchor item untested.
  - `"all-other"` is supported but gives the same z and p as `"centred"`.
- `dif_wald_sim.R`: null calibration and power. 48 cells:
  - n per group 150, 300, 600
  - 10 items × 4 categories, or 20 dichotomous items
  - impact 0 or 0.5
  - DIF: none, one item, three items in one direction, or two items in opposite directions, at 0.6 logits

  Each replication is seeded by cell and replication, so results do not depend on the number of cores.

Host run, from the package root (estimated about 70 minutes on 10 cores, about 250 MB in memory):

```
nohup Rscript dev/dif_lr_wald/dif_wald_sim.R 1000 10 > dev/dif_lr_wald/run.log 2>&1 &
```

The run writes `dif_wald_sim_<date>.rds`, `dif_wald_sim.rds` and `dif_wald_sim_summary.csv` here. The summary is also printed at the end of `run.log`.

The difference-plot sketch is in `../dif_lr_diffplot/`.

## Results (1000 replications per cell, run 2026-11-13)

The Monte-Carlo SE of a rate near .05 is about .007, so values from .036 to .064 are consistent with .05. "False flags" means any DIF-free item flagged after correction. Single-step and Holm gave nearly identical results, so single-step is shown.

**Null (no DIF).**
- The Wald z is calibrated for centred and forward (SD of z 0.98-1.03).
- After correction, the error rate is .033-.058 for centred and .039-.068 for forward. The highest values, .057 and .068, are both in the 10 x 4, n = 600, impact 0.5 cell.
- Constant (Gini, Strobl et al. 2021, a single anchor item by construction) is conservative: SD of z .82-.86 and an error rate of .008-.027.
- Without correction, at least one item is flagged in 38-43% of datasets with 10 items and 61-68% with 20.

**Artificial DIF, centred anchor.** It flags DIF-free items as soon as DIF is unbalanced, and more so with larger n.
- One DIF item, 10 x 4: .06, .09, .135 at n = 150, 300, 600.
- Three DIF items in the same direction, 10 x 4: .24, .47, .87. With 20 dichotomous items: .08-.155.
- With balanced DIF it holds (.02-.07).
- So it should not be used for item-level tests, which is the Hagquist and Andrich point in numbers.

**Forward (MTT) anchor.**
- It holds the error rate in every scenario. The range is .03-.07, apart from .083 for three DIF items at n = 150 with 10 x 4.
- It has the highest power for unbalanced DIF (e.g. .87 vs .73 centred at n = 300, 10 x 4).
- It matches centred for one DIF item, and is slightly below centred for balanced DIF (.70 vs .84 at n = 150).

**Constant (Gini) anchor.** It holds the error rate but has much lower power throughout (e.g. one DIF item, 20 dichotomous items, n = 300: .26 vs .56).

**Run notes.**
- Replications with an empty category in either group were skipped, up to 16% of them in 10 x 4 at n = 150 with impact 0.5. Results there are conditional on all categories being used, which RMdifLR() requires anyway.
- One replication in 48,000 failed inside multcomp ("'lower' not specified or contains NA"). A package version must catch that.
- A DIF item was the untested first anchor in up to 2.8% of replications for forward (n = 150) and 8.4% for constant.

**Not covered:** more than two groups, threshold-level DIF, non-uniform DIF, DIF sizes other than 0.6, and unequal group sizes.

## Literature check (2026-11-13, PDFs in articles/DIF/)

- **Kopf, Zeileis & Strobl (2015a), APM 39(2).** Introduces the iterative forward anchor class.
  - Design: dichotomous Rasch only. 40 items (20-80 checked), two groups, focal group mean -1, 15-45% DIF of ±0.6, item-wise Wald tests at .05 with no multiplicity correction.
  - Recommends forward with single-anchor (SA) selection when nothing is known about the direction of DIF.
  - All-other and forward held the level under the null. Constant anchors stayed below it.
- **Kopf, Zeileis & Strobl (2015b), EPM 75(1).** Introduces the MPT, MTT and Gini-type threshold selections.
  - Design: dichotomous Rasch only, 40 items.
  - Recommends MTT with the forward class, and MPT for a constant anchor of four.
  - Advises against all-other, which is inflated under unbalanced DIF.
  - Forward was robust to the selection strategy except with 40% unbalanced DIF. In that setting MTT was best.
  - Calls for future work on other IRT models.
- **Wang (2004), J Exp Educ 72(3).** Equal-mean-difficulty (= our "centred"), all-other and constant anchors.
  - Includes a partial credit case: 20 four-category items, ConQuest MML, step-level tests.
  - EMD and all-other are biased with inflated Type I error whenever the groups' mean item difficulties differ (unbalanced DIF), in the PCM too. The constant method with a correctly chosen DIF-free anchor works.
  - The anchors were known to be DIF-free (no selection strategy), so this supports anchoring in the PCM but not forward-MTT selection.
- **Hothorn, Bretz & Westfall (2008), Biom J 50(3)** (journal version, doi:10.1002/bimj.200810425). Single-step max-|t| adjusted p-values for asymptotically normal estimates. This is what `"fwer"` would use.
- **Glas & Verhelst (1995)** (book chapter, not downloaded). The source of the Wald test, cited by Kopf et al. and by psychotools.

What the literature does not cover, and the simulation here does: forward-MTT for partial credit item locations, and the single-step correction applied to the tests after anchor selection. Kopf et al. tested each item at .05 without correction. Our run is the only evidence that selection plus correction holds the family-wise rate.
- **Strobl, Kopf, Kohler, von Oertzen & Zeileis (2021), APM 45(3).** The Gini index as an anchor *point*: one shift that best overlaps the two groups' item parameters, chosen from a criterion plot. It does not choose a set of anchor items.
  - Design: dichotomous Rasch, two groups, 40 items, Kopf (2015a)-like setup, 10,000 replications. Per-item false-alarm and hit rates at .05, uncorrected.
  - With 20% unbalanced DIF, all anchors except all-other hold .05, and forward (MTT) has the highest hit rate.
  - With 40% unbalanced DIF, forward's false-alarm rate is "notably above" .05 at small and medium n. This matches our .083 for 3 of 10 items at n = 150.
  - With balanced DIF, Gini is outperformed by the traditional selection methods.
  - The authors stress that the "most items are DIF-free" assumption cannot be tested, and that a majority cluster with equal DIF is an equally valid solution. They recommend reporting the whole criterion plot.
  - Correction 2026-11-14: our "constant" arm (psychotools `class = "constant", select = "Gini"`, length 1) *is* Strobl's anchor point. Online Appendix F shows that every optimum of the Gini criterion falls where one item interlocks, so the anchor point is always a single anchor item. psychotools computes the Gini index of the absolute location differences for each candidate item and takes the maximum, which is the sparse grid of Appendix F. An earlier note here said the opposite, which was wrong. So our results are a direct evaluation of the Gini method for partial credit items with a correction: conservative (FWER .004-.027, .07 at n = 150 with 3 of 10 DIF items) and much lower power than forward. Strobl et al. find the same direction (Gini null false alarms about 2%, lower hit rates than forward, outperformed under balanced DIF) and attribute it to single-item anchors being more affected by sampling error.
  - Supplement (Online Appendices, PDF in articles/DIF/): the final tests leave untested the item with the smallest difference at the optimum. Sample sizes 250-1500 in total, 20/40/60 items, 0-60% DIF of ±0.6, focal mean -1. With 60% balanced DIF the Gini index has a strongly inflated false alarm rate, because the items form three clusters.
