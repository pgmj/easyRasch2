# Code review: R/item_parameters.R

Started 2026-11-10, main session. `/code-review` pass on this file not yet run. Checks in R in the container on eRm `raschdat1` (dichotomous, 30 items), `pcmdat2` (4 items, 0-2) and a mixed set (6 dichotomous + 4 polytomous items, n = 100). Single datasets.

Status key as in `utils-theta.md`. No fixes applied. The user commits the utils-theta.R work first.

Callers: `.rasch_fit_cml()` is used by RMitemParameters, RMpersonParameters, RMpersonFit, RMpersonChange, RMretestSD and RMscoreSE. `.rasch_fit_mml()` is used by the same functions except RMscoreSE. `.sparsity_warning()` is also called from cicc_plot.R, local_dependence.R and residual_pca.R.

## Summary

| # | Line | Finding | Status | Reaches users? |
|---|------|---------|--------|----------------|
| IP1 | 284 | `.rasch_fit_cml()` has no category guard (= R1 in utils-theta.md) | confirmed | yes, silent |
| IP2 | 339 | MML path reshapes items with empty categories the same way | confirmed | yes, silent |
| IP3 | 465 | Sparsity warning advises MML, which does not help with empty categories | confirmed | yes, wrong advice |
| IP4 | 373, 143 | MML SEs are not propagated through the centring | confirmed | yes, SEs 2-19% too large |
| IP5 | 297 | CML SEs and centring | no issue | |
| IP6 | 261, 264 | Stale comments: "(eRm)" header, "delta-method SEs" for CML | confirmed | no, doc only |
| IP7 | 445 | `.sparsity_warning()` duplicates `.check_cml_categories()` (= R6) | confirmed | no, maintenance |

## IP1. No category guard in `.rasch_fit_cml()` (confirmed)

The same psychotools behaviour as findings 2, 2b and 11 of utils-theta.md. `pcmodel()` drops the threshold of an empty interior category, rescales an item with an empty lowest category, and omits a constant item. `.rasch_fit_cml()` only calls `.sparsity_warning()`, so:

- **Empty interior category** (pcmdat2 item 1 recoded to 0/2): `RMitemParameters()` reports I1 with one row, "threshold 1", location 0.01, SE 0.099, for an item whose data run 0-2. `RMpersonParameters()` (CML) puts 56 respondents at the +10 boundary.
- **Empty lowest category** (item 1 recoded to 1/2): `RMpersonParameters()` puts 37 respondents at +10 (also reproduced by the diff review as R1, via `RMitemRestscore()`).
- **Constant item**: `setNames()` fails with "'names' attribute [4] must be the same length as the vector [3]".

## IP2. MML reshapes the same items (confirmed)

`mirt::mirt()` recodes item categories to be consecutive from 0, so `.rasch_fit_mml()` returns one threshold for the 0/2 item and for the 1/2 item, while the data keep the original codes. `RMpersonParameters(estimator = "MML")` puts 56 (interior) and 37 (lowest) respondents at +10, the same counts as CML. A constant item stops with mirt's clear "only one response category" message.

So the problem is how the data are represented downstream, not the estimator. A guard has to run for both estimators.

## IP3. The warning recommends MML for empty categories (confirmed)

`.sparsity_warning()` flags constant items and empty categories and says "CML estimates may be unstable; consider estimator = "MML"". By IP2, MML gives the same wrong person estimates for empty categories. The advice fits sparse but non-empty categories only.

## IP4. MML SEs ignore the centring (confirmed)

With `center = TRUE` (default), MML thresholds are shifted to a grand mean of zero, but `.mml_threshold_se()` propagates only each item's own block of the mirt covariance, on mirt's scale (latent mean fixed at 0). The SE of a centred threshold needs the full covariance, `C J V J' C'` with `C = I - 11'/n`. Computed that way:

| Data | Reported SE | SE of centred thresholds | Reported / centred |
|---|---|---|---|
| raschdat1 (dich) | 0.240-0.339 | 0.216-0.317 | 1.07-1.11 |
| pcmdat2 (poly) | 0.157-0.339 | 0.139-0.291 | 1.04-1.19 |
| mixed | 0.213-0.486 | 0.202-0.455 | 1.02-1.07 |

The centred dichotomous MML SEs (0.216-0.317) are close to the CML SEs (0.216-0.316), as expected for the same contrast. The current Wald CIs for MML are therefore too wide. The per-item mapping itself is correct: `sqrt(diag(J V J'))` equals the reported SEs in all three datasets, so the row alignment of `d` parameters is right, including mixed formats. With `center = FALSE` the current SEs are the correct ones for mirt's scale.

Proposed fix: build `J` over all items, form `J V J'`, and apply `C` when `center = TRUE`. That means `.rasch_fit_mml()` returns the full threshold covariance (or the centring happens inside it), since centring now happens in `RMitemParameters()`.

## IP5. CML SEs and centring (no issue)

`psychotools::threshpar()` already returns thresholds with grand mean 0 (both datasets), and the SEs of centred thresholds equal the reported ones exactly. Dichotomous SEs equal `eRm::RM()` and `psychotools::itempar()` sum-zero SEs (0.2161-0.3159). PCM SEs differ from `eRm::thresholds()` (0.142-0.356 against 0.148-0.337), as the code comment says. eRm's threshold SEs depend on its own constraint, so this is not evidence of an error here. The roxygen claim "CML estimates are already centred" holds.

## IP6. Stale comments (confirmed, doc only)

- Line 261: section header "CML fit + threshold extraction (eRm)". The code uses psychotools.
- Line 264: "(+ delta-method SEs)". The CML SEs come straight from `threshpar(vcov = TRUE)`, with no delta step.

## IP7. Duplicate category scan (= R6, maintenance)

`.sparsity_warning()` and the new `.check_cml_categories()` scan the same shapes, one to warn and one to stop.

## Proposed resolution (needs a decision)

IP1, IP2, IP3 and IP7 together with R1, R6 and R8 from utils-theta.md point to one change. A shared check (the current `.check_cml_categories()`) runs before every CML and MML fit, including the seven direct `pcmodel()` call sites R1 lists. `.sparsity_warning()` then keeps only the warn-level shapes: sparse but non-empty categories and extreme-only categories. Its MML advice applies only to those. Optionally, one shared CML fit helper (drop all-NA rows, check, fit, `threshpar()`, return centred thresholds plus optional vcov) replaces the direct calls, as R1 suggests.

IP4 is independent and local to this file. IP6 is a doc edit.

## Not yet done for this unit

- `/code-review high R/item_parameters.R` in a fresh session.
- Test gaps: no direct test of `.rasch_fit_cml()`; no test of MML SEs. Suggested: SE of centred MML thresholds against `C J V J' C'`; the IP1/IP2 shapes once guarded.
