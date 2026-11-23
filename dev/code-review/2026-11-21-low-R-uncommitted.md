# /code-review low R/ (2026-11-21)

Target: uncommitted diff of `R/` on top of 831dba8 (23 files, +656 / -249).
Level: low (single pass, no verification step).

## Result

No findings.

## Coverage

The diff passed to the reviewer was capped at 60 KB and stopped partway
through `R/item_cat_prob.R`. These files were **not reviewed**:

- R/item_restscore.R
- R/item_restscore_cutoff.R
- R/ld_partgam.R
- R/local_dependence.R
- R/martin_lof.R
- R/person_change.R
- R/person_fit.R
- R/reliability_curve.R
- R/residual_pca.R
- R/tile_plot.R
- R/utils-theta.R
- R/utils-validation.R

## Reviewer caveats, checked afterwards

- `data_mat` before `.extreme_only_categories()` in
  `R/bootstrap_restscore.R`: defined at line 187, used at 201. Fine.
- `n_sims` in `RMitemCatFitPlot()` for `type = "residuals"` with
  `fit = NULL`: initialised to `NA_integer_` at line 246 of
  `R/item_cat_fit_plot.R`, and every later use guards with `!is.na()` or
  is only reached once it has been set. Fine.

## Follow-up: /code-review low R/utils-theta.R R/utils-validation.R

No findings. The two helper diffs are small (+39 / -3 lines together), so
they were reviewed in full. Still not reviewed: item_restscore.R,
item_restscore_cutoff.R, ld_partgam.R, local_dependence.R, martin_lof.R,
person_change.R, person_fit.R, reliability_curve.R, residual_pca.R,
tile_plot.R.

## Follow-up: batch 1 (item_restscore, item_restscore_cutoff, ld_partgam, local_dependence, martin_lof)

No findings; all five diffs reviewed in full. Reviewer notes:

- item_restscore.R: `Flagged` and the caption now both use `alpha`
  instead of a hard-coded .05.
- item_restscore_cutoff.R, ld_partgam.R, local_dependence.R: only add
  `.check_iterations(iterations)` at the top.
- martin_lof.R: `n_failed` set in both branches; with sequential stopping
  it counts failures up to the stop, so early stopping is not counted as
  failure.

Still not reviewed: person_change.R, person_fit.R, reliability_curve.R,
residual_pca.R, tile_plot.R.

## Follow-up: batch 2 (person_change, person_fit, reliability_curve, residual_pca, tile_plot)

No findings. With this batch, every changed file in R/ has had a low-level
review.
