# Help-page review: plan

Started 2026-11-20. Goal: every exported help page is accurate against the
code, coherent and clear, in the same way as the pass on RMitemCatFit,
RMitemCatFitPlot, RMdifAnchor, RMdifImpact and RMdifLR (done 2026-11-20).

## Procedure per batch

1. Read the full roxygen block of each function and the code paths it
   describes (argument checks, defaults, stop/warning/message calls, output
   columns).
2. Check each page against the checklist below. Note every finding with the
   evidence from the code.
3. Edit the roxygen only. If the code looks wrong, or a fix would change
   behaviour, log it as a finding for the user and do not change the code.
4. `roxygen2::roxygenise()`, `tools::checkRd()` on each changed Rd, and
   `R CMD Rd2txt` to catch rendering artefacts (for example a wrapped line
   that starts with `+` or `-` becomes a bullet). Run
   `tests/testthat/test-documentation.R`, which guards against links to
   dotted helpers.
5. Log the batch in `dev/doc-review/<batch>.md`: a table of changes per
   function, findings left open, and decisions needed. Report to the user,
   who reviews before the next batch.

The container's C locale is enough for the checks. Run Rscript with
`LANG=C.UTF-8 LC_ALL=C.UTF-8`.

## Checklist

Accuracy

- Every default, allowed range and minimum (for example `iterations >= 20`)
  matches the code.
- Every error, warning and notice the user can meet is mentioned where it
  matters: empty categories, categories used only by extreme scorers,
  complete-case handling, dropped rows, low-iteration notices.
- Value: every `output` option is described, and the column names match the
  code.
- Statements about methods and simulation results match the cited source or
  the book chapter.

Clarity and coherence

- Terms are defined at first use (events, steps, padj, rest score, anchor).
- The same thing has the same name across a family of functions (test,
  cutoff, plot) and across families ("corrected p", "simulated datasets",
  "iterations").
- Interactions between arguments are stated (for example `fit` overriding
  `iterations` and `seed`).
- Sentences carry one idea each. Long parentheses are split out.

Structure

- Details use bold run-in headings where there are several topics, with a
  "Data" paragraph for missing data, extreme scores and category checks.
- Where a rasch_book chapter backs the function, the evidence paragraph
  links it with `\url{}` (never under `@references`). See the
  doc-review-checklist memory.
- References: every entry is cited in the text, every citation has an
  entry, entries are alphabetical and peer-reviewed where possible.
- See also links the family (test, cutoff, plot) and the natural next step
  (for example RMdifAnchor to RMdifImpact).
- Examples use the current API and run.

Style (user's writing rules)

- No em-dashes (`--` or `---` in prose), no semicolons in prose, no
  announcing sentences, "exactly" only when load-bearing.
- "Disordered", not "reversed", thresholds.
- Line length of roxygen comments within about 80 characters.

## Batches

Roxygen line counts in brackets. Done: RMitemCatFit, RMitemCatFitPlot,
RMdifAnchor, RMdifImpact, RMdifLR.

| Batch | Functions | Status |
|---|---|---|
| B1 Infit | RMitemInfit (227), RMitemInfitCutoff (120), RMitemInfitMI (141), RMitemInfitCutoffMI (115), RMitemInfitPlot (82), RMitemInfitCutoffPlot alias (7) | done 2026-11-20, log B1-infit.md |
| B2 Restscore and ICC | RMitemRestscore (198), RMitemRestscoreCutoff (138), RMitemRestscoreBoot (149), RMitemRestscorePlot (66), RMitemICCPlot (169) | done 2026-11-20, log B2-restscore-icc.md |
| B3 Dimensionality | RMdimCFA (89), RMdimCFACutoff (140), RMdimCFAPlot (44), RMdimResidualPCA (126), RMdimResidualPCACutoff (81), RMdimMartinLof (209), RMdimMartinLofResiduals (123) | done 2026-11-20, log B3-dimensionality.md |
| B4 Local dependence | RMlocdepQ3 (169), RMlocdepQ3Cutoff (127), RMlocdepQ3Plot (97), RMlocdepGamma (163), RMlocdepGammaCutoff (122), RMlocdepGammaPlot (87) | done 2026-11-20, log B4-local-dependence.md (incl. decision 3 package-wide) |
| B5 DIF, remaining | RMdifGamma (168), RMdifGammaCutoff (151), RMdifGammaPlot (73), RMdifTree (261) | done 2026-11-20, log B5-dif.md |
| B6 Categories, targeting, parameters | RMitemCatProb (141), RMitemHierarchy (78), RMtargeting (153), RMitemParameters (111), RMpersonParameters (125), RMscoreSE (95) | done 2026-11-20, log B6-categories-targeting-parameters.md |
| B7 Reliability | RMreliability (185), RMUreliability (35), RMreliabilityCurve (185) | done 2026-11-20, log B7-reliability.md |
| B8 Persons | RMpersonFit (193), RMpersonChange (296), RMretestSD (154) | done 2026-11-20, log B8-persons.md |
| B9 Descriptives and topics | RMplotBar (92), RMplotStackedbar (91), RMplotTile (109), phq9, easyRasch2-package, easyRasch2-reproducibility, easyRasch2-renaming | done 2026-11-20, log B9-descriptives-topics.md |
| B10 Cross-cutting | consistency across batches, README "Functions by domain", _pkgdown.yml, vignette passages | done 2026-11-20, log B10-cross-cutting.md, findings in FINDINGS.md |

Batches B1 to B5 are families of test, cutoff and plot, so consistency
within the family is the main check there. B7 and B8 are smaller in function
count but long and method-heavy.

## Decisions (user, 2026-11-20)

1. Titles in sentence case on all pages ("Conditional item infit MSQ").
   Keep acronyms and proper names (MSQ, DIF, CFA, Martin-Löf, Q3).
2. Convert every `--`/`---` dash and every semicolon in prose, also on
   pages not otherwise rewritten. (Semicolons in R code and in
   APA-style abbreviation lists stay.)
3. Switch the 16 run_* simulation helpers and the 2 `.format_*_label`
   helpers to `@noRd` (no help page, like the other internal helpers).
   Fold `print.RMlocdepGamma` and `knit_print.RMlocdepGamma` into the
   RMlocdepGamma page with `@rdname` and one sentence on the printed
   output. Done in B4 for the helpers in its files, the rest in B10.
