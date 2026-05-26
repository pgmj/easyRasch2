# knitr knit_print method for RMlocdepGamma kable output

Inside a knitr / Quarto / R Markdown chunk, returns the pre-combined
two-table asis string so pandoc renders them as two distinct pipe
tables. Outside knitr, R's normal dispatch falls back to
[`print.RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/reference/print.RMlocdepGamma.md).

## Usage

``` r
# S3 method for class 'RMlocdepGamma'
knit_print(x, ...)
```

## Arguments

- x:

  An object of class `"RMlocdepGamma"`.

- ...:

  Further arguments passed to
  [`knitr::asis_output()`](https://rdrr.io/pkg/knitr/man/asis_output.html).

## Value

A `knit_asis` character object.
