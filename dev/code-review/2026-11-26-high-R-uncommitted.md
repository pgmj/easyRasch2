# Code review: uncommitted changes after f44db3e (v9009)

Date: 2026-11-26. Level: high. Scope: the uncommitted working tree. The main
change is `dif_var` in RMitemCatFitPlot() (new `R/item_cat_fit_dif.R`). The
tree also has the RMdifGamma() asymptotic removal (reviewed in
`pg-dif-asymptotic-removal.md`) and the specificity study scripts.

Checked and found correct: the partial-gamma column changes in RMdifGamma()
and RMitemICCPlot() across all in-package callers, the within-score
permutation, the critical value against the permutation p, residual flat
indexing and the weighted PAV.

## Findings

| ID | Severity | File | Finding | Status |
|----|----------|------|---------|--------|
| CD1 | important | R/item_cat_fit_dif.R:444 | Literal `≥` in `paste0("X ≥ ", b$k)`. R CMD check warns about non-ASCII characters in R code, and the label is garbled in a non-UTF-8 locale. The non-DIF path (`item_cat_fit_plot.R:846`) uses `"\u2265"`. | fixed |
| CD2 | nit | R/item_cat_fit_plot.R:159 | The help calls the DIF band "a simultaneous 95% consistency band", but its coverage follows `conf_level`, as the caption does. | fixed |
| CD3 | nit | R/dif_tree.R:94 | The `p_adj` help points to "the asymptotic p-values of [RMdifGamma()]", which no longer exist. | fixed |
| CD4 | nit | R/item_cat_fit_dif.R:369 | In the DIF residual plot, groups differ by fill colour only. That fails the greyscale check, and the calibration and adjacent layouts also use line type and shape. | fixed: groups 2-4 get diagonal, cross-hatch and horizontal line patterns over their colour (grid patterns, ggplot2 >= 3.5.0, colours alone on older versions); `grid` added to Suggests; greyscale checked |
| CD5 | nit | R/item_cat_fit_dif.R:252 | `.catdif_bands()` holds several dense nrow(lay) x B matrices at once. 20 items x 4 thresholds x 4 groups at B = 1000 needs about 1 GB. The comment calls permutations "cheap". | fixed with option A: moments and maxima one curve at a time; peak 1197 MB to 332 MB at B = 1000, time 11.0 s to 10.3 s, output identical on all three layouts |
| CD6 | nit | R/item_cat_fit_plot.R:302 | `.prepare_dif_data()` runs before `validate_response_data()`, so bad `data` gives a misleading length error. The caption counts NA groups but not dropped all-empty rows. | fixed: validate_response_data() first; caption adds "complete cases" when empty rows were dropped |
| CD7 | nit | dev/category_conditional_null/R/summarise_specificity.R:261 | Formula `aggregate()` defaults to na.omit, so a dataset with any NA non-focus flag drops out of `any_oth`, and the merge then drops it from the rates too. The "datasets used" message does not show this. | checked, no effect: none of the 16767 used datasets has an NA infit or catfit flag, and the rates rest on 16767 datasets, the same number as the log. The script line is 47, not 261. Script left unchanged |

## CD5 measurement

Calibration layout, 20 items with 5 categories, 4 groups, n = 2000:
22,380 rows. Peak R memory 181 MB at B = 100 and 310 MB at B = 200, about
1.3 MB per permutation, so about 1.3 GB at the default B = 1000 (about 7
copies of the 179 MB permutation matrix). Time 1.7 s at B = 100, 2.6 s at
B = 200.

- Option A, one matrix: keep only the permutation matrix, compute the mean
  with rowMeans(na.rm = TRUE), the SD and the per-curve maxima column by
  column. Peak about 0.2 GB, same results, no extra time.
- Option B, two passes: keep no permutation matrix, store the permuted
  labels (n x B integers, 8 MB) and run the layout twice. Peak a few tens of
  MB, about twice the time.

Tests after the fixes: test-item_cat_fit_dif.R 62 pass, test-item_cat_fit.R
179 pass, 1 skip.

Option A implemented 2026-11-26. Same data at B = 1000: old 1197 MB peak,
11.0 s; new 332 MB, 10.3 s. `identical()` output from the old and new code
on calibration, adjacent and residuals (8 items, 3 groups, missing data,
B = 199). test-item_cat_fit_dif.R 62 pass.
