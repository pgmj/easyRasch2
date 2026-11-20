# B6 Categories, targeting, parameters: help-page review

2026-11-20. RMitemCatProb, RMitemHierarchy, RMtargeting, RMitemParameters,
RMpersonParameters, RMscoreSE. Roxygen only. Rd regenerated, checkRd clean,
Rd2txt scanned, test-documentation passes, all six examples run
(`run.donttest = TRUE`).

## Changes

| Page | Changes |
|---|---|
| all | Sentence-case titles ("Rasch or partial credit model" instead of "Rasch / Partial Credit Model"), dashes ("---", "--", em-dash) and prose semicolons removed, "e.g."/"i.e." mostly replaced, references alphabetical. |
| RMitemCatProb | Details said dichotomous items are fitted by "a Rasch model" with \eqn{\delta = -\beta}. The code fits one PCM (`.fit_cml_thresholds()`) for all items, fixed. "a la classic IRT trace plots" rephrased. See also adds RMtargeting and RMitemCatFit. |
| RMitemHierarchy | `data`: dichotomous data give an error (the code stops). Description and Details without "--" and "i.e.". See also adds RMitemParameters and RMitemCatProb. |
| RMtargeting | Wright & Stone (1979) was listed but not cited, now cited in the description. "reversal" now "disorder" (the user's wording rule, F27). "sit" (style list) replaced. `height_ratios` called the bottom panel "dot-whisker", stale since the category bands became the default. See also listed only `psychotools::pcmodel()` and `mirt::mirt()`; now the related easyRasch2 functions. |
| RMitemParameters | Andrich (1978) was listed but not cited, now cited for Andrich thresholds. Response categories paragraph adds the warning for categories used only by extreme scorers (`.sparsity_warning()`). Mair & Hatzinger (2007), the eRm paper, is still listed but not cited (F26). |
| RMpersonParameters | Bock & Mislevy (1982) and Magis (2015) were listed but not cited, now cited for EAP and for WLE with polytomous items. |
| RMscoreSE | Description now says the table applies to complete responses and points to RMpersonParameters for incomplete ones (the RMpersonParameters page said so, this page did not). Em-dash in Value removed. Warm and Bock & Mislevy cited in `method`. **No See also before**, now RMpersonParameters, RMitemParameters, RMreliabilityCurve. |

## Findings (not changed)

- **F26** RMitemParameters lists Mair & Hatzinger (2007, eRm) without
  citing it, and eRm is no longer used for estimation. Drop it?
- **F27** RMtargeting's caption says "Red arrows give the size in logits of
  each threshold reversal" (targeting_plot.R:527). Your wording rule is
  "disordered", not "reversed". The vignette passage I wrote copies it
  (easyRasch2.Rmd.orig:362). Change both to "threshold disorder"?
- **F28** RMitemCatProb tells categories apart by colour only (an ordered
  viridis ramp, no line type). The figure rule is colour plus line type
  plus shape per series. With 4 to 5 ordered categories the ramp may be
  enough, but greyscale and colour-vision readers lose the order of
  neighbouring curves.
- **F29** RMpersonParameters returns the EAP prior only as
  `attr(result, "prior")` (attributes rule).
