# Help-page review: findings for decision

2026-11-20. Collected from B1 to B10. Each finding names the batch log with
the details. Nothing here has been changed in code. "Doc" findings are
reference or wording questions only the user can settle.

## A. Plots that disagree with their tables (one fix pattern)

DONE 2026-11-20. All six plots draw the outer segment from the cutoff
object's stored interval (`item_cutoffs`, `pair_cutoffs`,
`loading_cutoffs`), the inner 66% segment and the dot shading use the same
method, and the no-data Q3 and gamma plots show 66% plus the stored width
instead of 66/95/99%. Captions name the interval and say the table flags by
adjusted p-values by default (decision 2a, RMdimCFA flags against the
interval so its caption omits that). Shared helpers in
`R/utils-plot-intervals.R`, tests in `tests/testthat/test-plot_intervals.R`,
NEWS. jamovi effects logged in `easyRasch2jmv/dev/PLAN-easyRasch2-next-release.md`
section G2 (two module analyses draw their own copies of the old plot).

| # | Function | Was |
|---|---|---|
| F3 | RMitemInfitPlot, RMitemRestscorePlot (with `data`) | DONE: central quantile intervals, table uses HDCI |
| F13 | RMdimCFAPlot | DONE: loadings segments fixed at 0.1/99.9%, ignored `percentile` |
| F20 | RMlocdepQ3Plot, RMlocdepGammaPlot, RMdifGammaPlot | DONE: fixed 66/95/99% HDCI without data, 66/99% quantiles with data |
| F4 | RMitemInfitPlot | DONE: without `data` an explicit `statistic` is followed, the default stays both panels |

## B. Output only in attributes (the 2026-11-15 rule)

DONE 2026-11-20. Decision: `output = "list"` is the house pattern, the
attributes stay for compatibility (the jamovi module reads `prior` and
`variance_partition`). Tests added to each function's test file. NEWS.

| # | Function | Outcome |
|---|---|---|
| F11 | RMitemICCPlot | DONE: no code change, @return says `dif_gamma` is `RMdifGamma(output = "dataframe")` on the same respondents and points there |
| F15 | RMdimResidualPCA | DONE: `output = "list"` with `components` and `variance` |
| F24 | RMdifTree | DONE: list gains `model`, `effect_size`, `stability`, `rescaled_terminal_items`, example uses the list |
| F29 | RMpersonParameters | DONE: left as attribute (user), @return now says what `prior` holds and that the CSV and histogram lack it |
| F30 | RMreliabilityCurve | DONE: `output = "list"` with `curve`, `summary` (one row) and `benchmark_range` |
| F34 | RMpersonChange | DONE: `output = "list"` with `persons` and `settings` (one row, adds `estimator` and `n_persons`) |

## C. Colour as the only cue

DONE 2026-11-21. F14: observed loadings all `sienna2`, filled diamond
(shape 18) when flagged, open diamond (shape 5) otherwise; caption and
@return updated. F28: new `vary_linetype = FALSE` argument, line type per
category from `scales::linetype_pal()` (13, recycled), merged legend, also
in path mode. F36: low counts bold as well as red. Tests added in all three
test files; NEWS entry added. Also, at the user's request: the CFA fit
plot uses the same filled/open diamonds in `sienna2` (was red vs grey30),
and the loadings plot uses `theme_bw()` (was `theme_minimal()`).

| # | Function | Problem |
|---|---|---|
| F14 | RMdimCFAPlot | flagged loading red vs orange (B3) |
| F28 | RMitemCatProb | categories by viridis colour only, no line type (B6) |
| F36 | RMplotTile | low counts red vs orange text (B9) |

## D. Error handling and input checks

DONE 2026-11-21 (all, user chose the wide scope for F1, F2, F16).
- F1: `.check_iterations()` (utils-validation.R) in all 13 functions with
  `iterations`; `min = 20` in RMitemCatFit(+Plot), `min = 0` in RMpersonFit
  (0 = statistics only). RMitemInfitCutoffMI also needs `iterations >= m`.
- F2: the extreme-only check sits in `.wle_theta_pool()`, the CML-DGP entry
  point of every `*Cutoff()`; the DIF gamma permutation DGP catches it and
  runs. RMitemRestscoreBoot checks before resampling. Found while checking:
  an extreme-only top category let 1 of 20 RMitemInfitCutoff iterations
  succeed and the function returned cutoffs from that one, silently.
  RMlocdepGammaCutoff (no refit) used to run on such data and now stops.
  RMdimMartinLof not included (CML on subscales, different failure mode).
- F6: requested_iterations, sample_n_total, sample_has_na, dgp added.
- F8: asymptotic path flags at `alpha`; caption shows it.
- F10: `.all_sims_failed_message(what = "bootstrap")`.
- F12: message with the number of dropped rows.
- F16: requested_iterations stored by PCA and CFA cutoffs, attrition
  clause in RMdimResidualPCA and both RMdimCFA tables. Q3 already did it.
  RMdimMartinLof: `requested_iterations`, `failed_iterations` and a message
  (sequential stopping makes actual < requested by design).

| # | Function | Problem |
|---|---|---|
| F1 | RMitemInfitCutoff | no validation of `iterations` |
| F2 | RMitemInfitCutoff | no category check before simulating (fails only after all iterations) |
| F10 | RMitemRestscoreBoot | "All bootstrap iterations failed" does not name extreme-only categories |
| F12 | RMitemICCPlot | drops incomplete rows without a message |
| F16 | RMdimResidualPCACutoff, RMdimMartinLof | drop failed simulations without a message |
| F6 | RMitemInfitCutoffMI | result lacks `requested_iterations`, `sample_n_total`, `sample_has_na`, `dgp` |
| F7 | RMitemInfit | DONE 2026-11-20: refuses an RMitemInfitCutoffMI result (error), test added, NEWS |
| F8 | RMitemRestscore | asymptotic test flags at a fixed .05, ignores `alpha` |

## E. Consistency across functions

| # | Topic | Detail |
|---|---|---|
| F5 | `n_cores = NULL` without `mc.cores` | cutoff functions warn and run sequentially, RMitemCatFit and RMpersonFit use two workers, the reliability and change functions run sequentially without a warning |
| F9 | sign of the restscore difference | RMitemRestscoreBoot `diff` is expected minus observed, the others observed minus expected (documented, changing it breaks user code) |
| F19 | `p_value` default | RMdimCFA `FALSE`, infit/restscore/Q3/gamma `NULL` (use p-values when the full cutoff object is given) |
| F21 | `Flagged` type | `"overfit"`/`"underfit"` (infit, restscore), `"above"`/`"below"` (Q3, CFA loadings), logical (gamma LD, gamma DIF, residual PCA, DIF tree) |
| F27 | "reversal" | DONE 2026-11-20: RMtargeting caption says "threshold disorder", NEWS |
| F18 | fixed cutoff | RMdimMartinLofResiduals flags cells at fixed abs(residual) > 2 |

## F. References (doc decisions)

| # | Page | Question |
|---|---|---|
| F17 | RMdimCFACutoff | DONE 2026-11-20: Yuan & Bentler (2000) removed |
| F22 | RMdifGammaCutoff | DONE 2026-11-20: Henninger et al. (2025) cited next to Bjorner et al. in the description |
| F23 | RMdifTree | DONE 2026-11-20: new Performance paragraph citing Asamoah et al. (2025), from the abstract the user supplied |
| F25 | RMdifGamma | DONE 2026-11-20: class list removed, one sentence points to RMdifTree and says the classes are conventions the function does not use |
| F26 | RMitemParameters | DONE 2026-11-20: Mair & Hatzinger (2007) removed |
| F31 | RMreliability | DONE 2026-11-20: Adams (2005) removed |
| F33 | RMpersonFit | DONE 2026-11-20: de la Torre & Deng (2008) cited in the lz paragraph |
| -- | RMdimMartinLof, RMdifTree | CONFIRMED by the user 2026-11-20: Tjur (1982) Scand J Stat 9(1):23-30; Strobl et al. (2015) doi 10.1007/s11336-013-9388-3; Komboz et al. (2018) doi 10.1177/0013164416664394 |

## G. Examples and small items

| # | Item |
|---|---|
| F32 | RMUreliability has no example |
| F35 | RMpersonChange and RMretestSD examples use `easyRasch2:::sim_partial_score()` |
| F37 | easyRasch2-reproducibility See also points to the renaming page |
