# B4 Local dependence: help-page review

2026-11-20. RMlocdepQ3, RMlocdepQ3Cutoff, RMlocdepQ3Plot, RMlocdepGamma,
RMlocdepGammaCutoff, RMlocdepGammaPlot, plus decision 3 of the plan
(internal help pages). Roxygen only. Rd regenerated, checkRd clean, Rd2txt
scanned, test-documentation passes, all six examples run
(`run.donttest = TRUE`). NAMESPACE unchanged (the S3 registrations stay).

## Decision 3 (done for the whole package, not only B4)

- All 16 `run_*` simulation helpers and the 2 `.format_*_label` helpers
  switched from `@keywords internal` to `@noRd`. Their 18 Rd files were
  deleted by roxygen. Nothing in man/, `_pkgdown.yml` or the vignette
  linked to them. (Other run_* helpers, for CFA, PCA, Martin-Löf,
  reliability and person change, were already `@noRd`.)
- `print.RMlocdepGamma` and `knit_print.RMlocdepGamma` folded into the
  RMlocdepGamma page with `@rdname`. Their two Rd files were deleted. The
  page now documents `x` and `...` and shows both methods in Usage, and
  `?print.RMlocdepGamma` opens it. The Value section says how the object
  prints.

## Changes

| Page | Changes |
|---|---|
| all | Sentence-case titles (`\eqn{Q_3}` removed from titles, "Yen's Q3"), `[fn()]` links, dashes, prose semicolons and "->" in example comments removed, "credible range/interval" (a Bayesian term) replaced by "simulated interval", example comments say 1000 or more iterations (they said "500+" and "100+", while the pages recommend 1000 to 2000). |
| RMlocdepQ3 | Description says it tests pairs with p-values given the cutoff object. `data`: Q3 is pairwise complete, all-NA respondents dropped (the page said only "NA allowed"). `cutoff` as a three-item list. `n_pairs`: the correction covers all pairs. Typo "\eqn{Q_3}matrix". "true to the Rasch tradition" removed. Details given Statistic, Estimation, Two views and Bootstrap p-values headings. Ferreira (2024) cited in the text. References alphabetical. **No See also before**, now RMlocdepQ3Cutoff, RMlocdepQ3Plot, RMlocdepGamma. |
| RMlocdepQ3Cutoff | Description now says what is returned and what it is for. Johansson (2026) and Christensen et al. (2017) were cited or used but the page had **no references**, both added. "see the package's comparison study" (not findable by users) removed. "Silently discarded" replaced (RMlocdepQ3 reports attrition). `requested_iterations` added to Value (returned, not documented). Redundant "Supports both dichotomous and polytomous" paragraph removed. |
| RMlocdepQ3Plot | Description said the observed Q3 comes from `mirt::residuals()`. It uses the estimator of `simfit` (CML by default). The page said the segments come from `simfit$pair_cutoffs`; they are fixed 66% and 99% quantile intervals, and without data fixed 66/95/99% HDCI (F20). Says the plot shows no p-values. |
| RMlocdepGamma | Description says it tests pairs with p-values given the cutoff object. `data`: complete cases. `gamma_pair` added to the data.frame columns (always returned, not documented). Details split into Statistic, Asymptotic p-values and Bootstrap p-values. Johansson (2026) was cited but missing from the references, added. Westfall & Young cited. See also adds RMlocdepQ3. Print methods folded in (above). |
| RMlocdepGammaCutoff | Description no longer says the cutoffs are "via iarm::partgam_LD()" (the simulations use an internal routine). The two overlapping Details paragraphs merged. `requested_iterations` added to Value. Johansson (2026) added to the references (cited). "Silently discarded" replaced. |
| RMlocdepGammaPlot | Description said the observed values come from `iarm::partgam_LD()`; they come from the same internal routine as the simulations (iarm is only checked for). Segment description corrected as for Q3 (F20). iarm needed only with `data`. |

## Findings in the code (not changed)

- **F20** RMlocdepQ3Plot and RMlocdepGammaPlot draw fixed interval widths:
  66/95/99% HDCI without `data`, 66% and 99% quantile intervals with
  `data`. Neither follows `hdci_width` or `pair_cutoffs`, which the tables
  use. Same family as F3 (infit, restscore plots) and F13 (CFA plot).
- **F21** `Flagged` has three types across families: character
  `"overfit"`/`"underfit"` (infit, restscore), `"above"`/`"below"` (Q3,
  CFA loadings) and logical (RMlocdepGamma, RMdimResidualPCA). For B10.
