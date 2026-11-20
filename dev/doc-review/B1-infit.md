# B1 Infit: help-page review

2026-11-20. RMitemInfit, RMitemInfitCutoff, RMitemInfitMI,
RMitemInfitCutoffMI, RMitemInfitPlot (+ deprecated RMitemInfitCutoffPlot).
Roxygen only. Rd regenerated, checkRd clean (only the usual non-ASCII note
for "Müller", the package declares UTF-8), Rd2txt scanned, test-documentation
passes, all five examples run (`run.donttest = TRUE`).

## Changes

| Page | Changes |
|---|---|
| all | Titles in sentence case. `\code{\link{}}` replaced by `[fn()]` links. `--` dashes and prose semicolons removed. "Suggests, not Imports" shortened. |
| RMitemInfit | Description now says what the function is for (flagging, p-values by default). `cutoff`: stale "identical to the current version" removed, the two forms explained by what they enable. `p_value`: tightened, version note kept in one sentence. `correction` names BH and BY. `alpha`: "used only with p-values". `sort`: says other values are ignored. `statistic`: 2k-family sentence clarified. Value: stale "(not the previous logical)" removed, kable p columns added, `p_infit` is "uncorrected" not "marginal". Details split into Statistic, Data and model (complete cases for infit, category check stops naming the item, relative location from all respondents' WLE), Bootstrap p-values, Iterations. Westfall & Young (1993) now cited in the text. References alphabetical. See also adds RMitemInfitPlot, RMitemInfitMI, RMitemRestscore. |
| RMitemInfit, inherited | "Multiple comparisons" section rewritten without dashes. It is inherited by RMitemRestscore, RMlocdepQ3, RMlocdepGamma, RMdifGamma and RMdimCFA, so it stays general ("comparisons", "item or item pair"). Those Rd files changed with it. |
| RMitemInfitCutoff | Description says what the result is used for. `iterations`: 400 is the recommended minimum, 1000 to 2000 for a final analysis. `cutoff_method`: needs ggdist. `dgp` rephrased. Details: "silently discarded" was wrong (RMitemInfit reports attrition in its caption); now says what happens when none can be refitted. Johansson (2025) was in the references but not cited, now cited. **Example bug fixed**: "Use the cutoffs in RMitemInfit()" called `RMitemInfit(sim_data)` without the cutoff. |
| RMitemInfitMI | Description says it flags against the interval. `mids_object`: factor items and the level requirement (`.mi_completed_to_numeric()`). Stale "not the previous logical" removed. SE caveat split into sentences. Johansson (2026) was cited in a section but missing from the references, added. Example comment "250+ iterations" contradicted the 0.999 width, now points to RMitemInfitCutoffMI. |
| RMitemInfitCutoffMI | Value no longer claims "the same structure" as RMitemInfitCutoff (it lacks `requested_iterations`, `sample_n_total`, `sample_has_na`, `dgp`) and no longer says it can be passed to RMitemInfit (see F7). Em-dash removed. Interval-width section now says the 0.999 interval needs roughly 5000 iterations and gives its family-wise error rate. Johansson (2026) added to the references (cited, was missing). Example comment fixed as above. |
| RMitemInfitPlot | `statistic`: without `data` it is ignored and both are shown (F4). Details rewritten to match the code: without `data` the intervals are HDCI (66% and the simfit width), with `data` the segments are central quantile intervals, not "the 66% HDCI from item_cutoffs" (F3). Says the plot shows no p-values. |

## Findings in the code (not changed)

- **F1** RMitemInfitCutoff does not validate `iterations` (0, negative or
  non-numeric give base R errors from `sample.int()`). RMitemCatFit checks
  `>= 20`.
- **F2** RMitemInfitCutoff has no category check on the observed data before
  simulating. An empty or extreme-only category is reported only after all
  iterations have failed (`.all_sims_failed_message()`), which wastes the
  run. RMitemInfit checks first (`.check_complete_case_categories()`).
- **F3** RMitemInfitPlot with `data` draws central quantile intervals, while
  RMitemInfit tabulates the HDCI. The outer segment could use
  `simfit$item_cutoffs`, so plot and table show the same bounds.
- **F4** RMitemInfitPlot ignores `statistic` when `data` is missing (both
  panels always). Could honour it.
- **F5** `n_cores = NULL` with no `mc.cores`: RMitemInfitCutoff warns and runs
  sequentially, RMitemCatFit uses two workers. Cross-package consistency,
  for B10.
- **F6** RMitemInfitCutoffMI's result lacks `requested_iterations`,
  `sample_n_total`, `sample_has_na` and `dgp`, so callers cannot report
  attrition or the complete-case sample for it.
- **F7** RMitemInfit accepts a RMitemInfitCutoffMI result (it has
  `$results`), and with `p_value = NULL` then computes p-values of the
  complete-case infit against simulations stacked over imputations. That
  combination is not calibrated (Johansson 2026 covers complete data). Should
  RMitemInfit refuse or warn for an MI cutoff?

## For B10

- "marginal p-value" is now "uncorrected p-value" in RMitemInfit and the
  inherited section. Other pages and column descriptions may still say
  "marginal".
