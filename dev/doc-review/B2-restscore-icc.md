# B2 Restscore and ICC: help-page review

2026-11-20. RMitemRestscore, RMitemRestscoreCutoff, RMitemRestscoreBoot,
RMitemRestscorePlot, RMitemICCPlot. Roxygen only. Rd regenerated, checkRd
clean (only the UTF-8 note for the en dash in Kreiner's title), Rd2txt
scanned, test-documentation passes, all five examples run
(`run.donttest = TRUE`).

## Changes

| Page | Changes |
|---|---|
| all | Sentence-case titles, `[fn()]` links, dashes and prose semicolons removed, "percent" as "%" in running numbers, "Suggests, not Imports" shortened, reference lists alphabetical. |
| RMitemRestscore | Title "Item-restscore association". Description said "absolute difference", but `Difference` is signed. It now also says the bootstrap test is recommended. `alpha`: the asymptotic test always flags at .05 (F8). `sort`: other values ignored. Value: "marginal" now "uncorrected". Details given Statistic and Data and model headings (complete cases, category stop, relative location). See also adds RMitemRestscoreBoot and RMitemInfit. |
| RMitemRestscoreCutoff | Description says what the result is used for and cites Kreiner (2011), which was listed but not cited. Johansson (2026) was listed but not cited; now cited where the interval-versus-p-value argument is made, "as shown for infit", since the preprint covers infit. "Silently discarded" was wrong (RMitemRestscore reports attrition), fixed. Details given Why simulate and Simulation headings. |
| RMitemRestscoreBoot | Title now says what it is: "Bootstrap flag rates of the asymptotic item-restscore test". The long second Description paragraph moved to Details ("What the percentage means"). `cutoff`: the kable keeps items with more than `cutoff` percent (`>`), the old text implied `>=`. Value: `diff` is expected minus observed, the opposite sign of RMitemRestscore's `Difference`, now stated (F9). History about the eRm engine removed. **Example** used `theme_minimal()` and colour alone. Now `theme_bw()`, Okabe-Ito colours plus shape, and the axis label says gamma, not "correlation". |
| RMitemRestscorePlot | Same as RMitemInfitPlot in B1: without `data` the intervals are HDCI, with `data` the segments are central quantile intervals (the old text said the shaded interval matches the table in both cases). Says the plot shows no p-values. |
| RMitemICCPlot | New Data paragraph (complete cases, rows with missing `dif_var` dropped, empty category stops, extreme-only category warns). Andersen (1995) was listed but not cited, now cited in Conditioning. "sit" replaced (style list). Dashes removed. See also adds the item fit tests, RMdifAnchor and RMitemCatFitPlot. |

## Findings in the code (not changed)

- **F8** RMitemRestscore without a `cutoff` flags at a fixed .05
  (`p_adjusted < 0.05`) and ignores `alpha`. Documented. Could use `alpha`.
- **F9** RMitemRestscoreBoot's raw output has `diff = expected - observed`,
  the opposite sign of `Difference` in RMitemRestscore and
  RMitemRestscoreCutoff. Documented. Changing it would break user code that
  reads `diff`.
- **F10** RMitemRestscoreBoot stops with "All bootstrap iterations failed.
  Check your data." without naming categories used only by extreme scorers,
  unlike the other simulation functions (`.all_sims_failed_message()`).
- **F11** RMitemICCPlot returns the DIF partial gamma table only as
  `attr(., "dif_gamma")`. That conflicts with the rule that nothing users
  need lives only in attributes (no-output-in-attributes memory,
  2026-11-15). A list output with the table would fit the rule.
- **F12** RMitemICCPlot drops rows with missing responses or missing
  `dif_var` without a message. Other functions say how many rows they drop.
- RMitemRestscorePlot with `data` has the same quantile-versus-HDCI mismatch
  as F3 in B1.
