# B5 DIF (remaining): help-page review

2026-11-20. RMdifGamma, RMdifGammaCutoff, RMdifGammaPlot, RMdifTree. Roxygen
only. Rd regenerated, checkRd clean, Rd2txt scanned, test-documentation
passes, all four examples run (`run.donttest = TRUE`).

## Changes

| Page | Changes |
|---|---|
| all | Sentence-case titles, `[fn()]` links instead of `\code{\link{}}`, dashes and prose semicolons removed, "credible range" now "simulated interval", "Monte-Carlo" spelled out, `n_cores` wording aligned with the other batches, example comments recommend 1000 or more iterations. |
| RMdifGamma | Description says it tests items with p-values given the cutoff object. `data` and `dif_var`: complete cases in both, two or more groups. `cutoff` and `p_value` restructured like the infit and LD pages. Bjorner's three DIF classes kept but introduced as "Bjorner et al. (1998) proposed" (was "Recommended interpretive thresholds") and followed by a sentence that they are conventions, not calibrated (F25). Details split into Statistic, Asymptotic p-values, Bootstrap p-values, Real and artificial DIF, which now points to RMdifAnchor. Westfall & Young cited. References alphabetical (Johansson was last). See also adds RMdifGammaPlot, RMdifAnchor, RMdifImpact. |
| RMdifGammaCutoff | Description cites Bjorner et al. (1998), listed but not cited. Percent signs, "null" wording simplified. Henninger et al. (2025) is still listed but not cited (F22). |
| RMdifGammaPlot | Same correction as the other simulation plots: the segments with `data` are fixed 66% and 99% quantile intervals, not `item_cutoffs`, and without `data` the HDCI widths are fixed at 66/95/99% (F20). Says the plot shows no p-values. |
| RMdifTree | Rewritten without the many "--" asides. The classification credit was muddled ("ETS A/B/C categories (Bjorner et al., 1998)"); now ETS for MH and Bjorner for partial gamma. Holland & Thayer was cited as 1988 in Details but listed as 1986, now 1986 throughout. Strobl et al. (2015) and Komboz et al. (2018), the raschtree and pctree methods, were only in the source-file header; now cited and listed (DOIs added from memory, **please check**: 10.1007/s11336-013-9388-3 and 10.1177/0013164416664394). Henninger et al. (2023) now cited (MH effect size in trees). `covariates`: rows with a missing covariate are dropped with a message, at least 30 must remain. `output` options moved to Value. See also adds RMdifAnchor and RMdifImpact. |

## Findings (not changed)

- **F22** RMdifGammaCutoff lists Henninger et al. (2025) but never cites it.
  Cite it where it supports something (partial gamma for DIF?), or drop it.
- **F23** RMdifTree lists Asamoah et al. (2025, item balance and effect
  sizes with the Rasch tree) but never cites it. A sentence on what it found
  would earn its place, or drop it.
- **F24** RMdifTree returns the stability results only as attributes
  (`stability`, `stability_kable`), against the attributes rule. With
  `output = "list"` they could be list elements.
- **F25** RMdifGamma's page lists Bjorner's fixed DIF classes. They are
  framed now as a proposal and as uncalibrated conventions, matching the
  caveat in RMdifTree. Given the package line on fixed cutoffs, keep, or drop
  them from RMdifGamma (RMdifTree uses them in code, so they stay there)?
- RMdifGammaPlot has the fixed-interval issue F20.
