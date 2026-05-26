# Print method for RMlocdepGamma kable output

Prints the two rest-score direction tables stacked vertically with a
blank line between them. Each table renders via `knitr_kable`'s own
print method as a clean pipe-markdown table.

## Usage

``` r
# S3 method for class 'RMlocdepGamma'
print(x, ...)
```

## Arguments

- x:

  An object of class `"RMlocdepGamma"` returned by
  [`RMlocdepGamma`](https://pgmj.github.io/easyRasch2/reference/RMlocdepGamma.md)
  with `output = "kable"`.

- ...:

  Further arguments (currently unused).

## Value

Invisibly returns `x`.
