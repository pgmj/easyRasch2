# easyRasch2: Psychometric Analysis with Rasch Measurement Theory

Streamlines reproducible Rasch measurement theory analyses for ordinal
item-response data, combining estimation routines from 'eRm',
'psychotools', 'mirt', 'iarm', and 'lavaan' with consistent diagnostic,
plotting, and reporting layers. Covers the four basic psychometric
criteria summarised by Christensen et al. (2021)
[doi:10.1111/sms.13908](https://doi.org/10.1111/sms.13908) –
unidimensionality, local independence, ordered response category
thresholds, and invariance across subgroups – together with item fit,
targeting, reliability, category functioning, and descriptive
item-response plots. A distinguishing feature is the use of
simulation-based critical values to replace rule-of-thumb cutoffs for
conditional infit mean-square, Yen's Q3 local-dependence statistic, the
largest residual-PCA eigenvalue, ordinal CFA fit indices, and
partial-gamma DIF and local-dependence coefficients, optionally
augmented with multiplicity-corrected bootstrap p-values. Outputs are
knitr::kable() tables and 'ggplot2' figures suitable for direct
inclusion in 'Quarto' and 'R Markdown' reports.

## See also

Useful links:

- <https://github.com/pgmj/easyRasch2>

- <https://pgmj.github.io/easyRasch2/>

- Report bugs at <https://github.com/pgmj/easyRasch2/issues>

## Author

**Maintainer**: Magnus Johansson <pgmj@pm.me>
([ORCID](https://orcid.org/0000-0003-1669-592X))

Authors:

- Magnus Johansson <pgmj@pm.me>
  ([ORCID](https://orcid.org/0000-0003-1669-592X))

Other contributors:

- Nicklas Korsell (PCM simulation code) \[contributor\]

- Mirka Henninger ([ORCID](https://orcid.org/0000-0003-4676-2361)) (MH /
  partial-gamma effect-size and ETS-classification algorithms in
  dif_tree.R, adapted under MIT licence from the raschtreeMH and
  effecttree packages) \[contributor\]

- Jan Radek ([ORCID](https://orcid.org/0009-0003-8842-9206))
  (partial-gamma effect-size and ETS-classification algorithms in
  dif_tree.R, adapted under MIT licence from the effecttree package)
  \[contributor\]
