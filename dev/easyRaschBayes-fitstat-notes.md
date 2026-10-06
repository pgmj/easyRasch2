# easyRaschBayes: two issues in the item fit functions

Found 2026-10-04 while reading easyRaschBayes 0.3.0 for the large-sample item
fit thread. Not changed. These belong to the easyRaschBayes package and need
the user's decision.

## 1. `infit_statistic()` conditions on theta, not on the total score

`R/brms_fitstats.R`, roxygen around lines 594-640 and the function at 708.

- The documentation describes the procedure as adapting "the conditional
  infit/outfit statistics (Christensen et al., 2013; Kreiner & Christensen,
  2011; Müller, 2020)". The `infit_post()` plot labels the x-axis "expected
  conditional item infit" (`R/post_helpers.R`).
- The computation takes `E` and `Var` for each observation from
  `brms::posterior_epred()`. Those are expectations given each person's
  theta draw. In the cited references, "conditional" means conditional on
  the total score: CML item parameters and expectations from the elementary
  symmetric functions, with no person parameters.
- So the statistic is a theta-based (unconditional) infit, evaluated draw by
  draw, with a posterior predictive reference. As a posterior predictive
  check it is still valid, but the label may mislead, and values will not
  match `easyRasch2::RMitemInfit()` or `iarm::out_infit()`.
- Options: relabel it as posterior predictive infit and drop "conditional",
  or implement total-score conditioning per draw. Under the Rasch model the
  conditional expectations depend only on the item parameters, so they could
  be computed from each draw's thresholds with an ESF.

## 2. The observed-gamma interval in `item_restscore_statistic()` has zero width

`R/brms_fitstats.R`, function at line 1047. The relevant lines are about 73-87
and 125-137 within the function.

- `gamma_obs` is computed once from the observed data and then repeated
  across all draws (`matrix(gamma_obs, nrow = n_draws, byrow = TRUE)`). So
  `gamma_obs_q025 == gamma_obs_q975` and `gamma_obs_q005 == gamma_obs_q995`.
  The documentation calls these "credible interval for the observed gamma".
- `gamma_diff_q025` and `gamma_diff_q975` are therefore `gamma_obs` minus
  quantiles of the *replicated* gamma. That is a posterior predictive
  interval for the gamma of a new dataset of the same size. It is not a
  credible interval for a misfit parameter, and it narrows at about
  1/sqrt(n). Read as "misfit size with uncertainty", it would flag trivial
  misfit at large n.
- Options: drop the `gamma_obs_q*` columns, or document them as the
  observed value. Rename the `gamma_diff` interval as a predictive interval,
  and say in the docs that ppp-based and interval-based rules behave like
  tests as n grows.

## Related context

Both functions are posterior predictive checks. Neither gives a posterior
for the size of the misfit, so a ROPE cannot be applied to them directly.
The design in `slope-rope-design.md` puts the misfit into the model as a
parameter.
