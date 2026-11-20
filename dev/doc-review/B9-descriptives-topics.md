# B9 Descriptives and topics: help-page review

2026-11-20. RMplotBar, RMplotStackedbar, RMplotTile, phq9, easyRasch2-package,
easyRasch2-reproducibility, easyRasch2-renaming. Roxygen only. Rd
regenerated, checkRd clean, Rd2txt scanned, test-documentation passes, the
three plot examples run.

## Changes

| Page | Changes |
|---|---|
| RMplotBar, RMplotStackedbar | Sentence-case titles, shorter descriptions. The "Input requirements" list became one sentence. Caption paragraph without "--" and semicolon. RMplotStackedbar `show_percent` said percentages are shown "instead of (or in addition to) counts", now says that with `show_n = TRUE` the label is "n (p%)", as the code does. |
| RMplotTile | Sentence-case title, description and Details without "--" asides. **No See also before**, now the other descriptive plots and the DIF tests. |
| phq9 | "Respondents with a sum score of zero were reduced by hand to 8" now "The number of respondents with a sum score of zero was reduced to 8". One semicolon. The em-dashes in the q6 and q8 descriptions are the instrument's own wording and stay. |
| easyRasch2-reproducibility | The list of functions the page applies to **lacked RMpersonChange (with `critical = "simulate"`), RMreliabilityCurve and RMretestSD**, which take a `seed` and pin the generator per iteration in the same way. Added, list alphabetical. One semicolon. |
| easyRasch2-renaming | Said "22 exported functions were renamed"; the table has 25 rows, now 25. "---" dashes replaced. The "No deprecation aliases" section now says the old names were removed without aliases, and that `RMitemInfitCutoffPlot()` (the new name of `RMinfitcutoffPlot()`) has since become a deprecated alias of `RMitemInfitPlot()`. "ahead of the first CRAN submission" now "before the first CRAN release". |
| easyRasch2-package | Generated from DESCRIPTION (updated earlier today, estimation packages), nothing else to change. |

## Findings (not changed)

- **F36** RMplotTile highlights low counts by text colour alone, red against
  the default orange. Red and orange are hard to tell apart with colour
  vision deficiency. A bold face or a marker would carry it.
- **F37** easyRasch2-reproducibility has `@seealso [easyRasch2-renaming]`,
  which is unrelated. Point it to [set.seed()] and
  [withr::with_preserve_seed()] instead?
