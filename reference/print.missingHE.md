# Print method for objects in the class `missingHE`

Produces a printout of the (posterior mean) estimate for each regression
coefficient included in the outcome and missingness/structural value
models based on the health economic evaluation probabilistic model run
using the function
[`selection()`](https://angabrio.github.io/missingHE/reference/selection.md),
[`pattern()`](https://angabrio.github.io/missingHE/reference/pattern.md),
[`hurdle()`](https://angabrio.github.io/missingHE/reference/hurdle.md)
or [`lmdm()`](https://angabrio.github.io/missingHE/reference/lmdm.md).

## Usage

``` r
# S3 method for class 'missingHE'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A `missingHE` object containing the results of the Bayesian modelling
  and the economic evaluation.

- digits:

  Integer indicating the number of decimal places to be used for
  rounding (default = 3).

- ...:

  Additional arguments affecting the summary produced.

## Value

Prints all fixed regression coefficients based on the type of model
fitted using the function
[`selection()`](https://angabrio.github.io/missingHE/reference/selection.md),
[`pattern()`](https://angabrio.github.io/missingHE/reference/pattern.md),
[`hurdle()`](https://angabrio.github.io/missingHE/reference/hurdle.md)
or [`lmdm()`](https://angabrio.github.io/missingHE/reference/lmdm.md).
Summary estimates based on posterior means of each fixed regression
coefficient included into the fitted models are provided.

## References

Baio, G.(2012). *Bayesian Methods in Health Economcis*. CRC/Chapman
Hall, London.

## See also

[`selection()`](https://angabrio.github.io/missingHE/reference/selection.md)
[`pattern()`](https://angabrio.github.io/missingHE/reference/pattern.md)
[`hurdle()`](https://angabrio.github.io/missingHE/reference/hurdle.md)
[`lmdm()`](https://angabrio.github.io/missingHE/reference/lmdm.md)
[`diagnostic()`](https://angabrio.github.io/missingHE/reference/diagnostic.md)
`print.missingHE()`

## Author

Andrea Gabrio

## Examples

``` r
# For examples see the function [selection()], [pattern()], [hurdle()] or [lmdm()].
#
#
```
