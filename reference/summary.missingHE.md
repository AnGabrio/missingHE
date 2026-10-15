# Summary method for objects in the class `missingHE`

Produces a table printout with some summary results of the health
economic evaluation probabilistic model run using the function
[`selection`](https://angabrio.github.io/missingHE/reference/selection.md),
[`pattern`](https://angabrio.github.io/missingHE/reference/pattern.md),
[`hurdle`](https://angabrio.github.io/missingHE/reference/hurdle.md) or
[`lmdm`](https://angabrio.github.io/missingHE/reference/lmdm.md).

## Usage

``` r
# S3 method for class 'missingHE'
summary(object, incremental = FALSE, prob = c(0.025, 0.975), digits = 3, ...)
```

## Arguments

- object:

  A `missingHE` object containing the results of the Bayesian modelling
  and the economic evaluation

- incremental:

  Logical. If `incremental` is `TRUE`, incremental CE results are
  printed.

- prob:

  A numeric vector of probabilities within the range (0, 1),
  representing the upper and lower CI quantiles to be calculated and
  returned for the posterior estimates.

- digits:

  Integer indicating the number of decimal places to be used for
  rounding (default = 3).

- ...:

  Additional arguments affecting the summary produced.

## Value

Prints tables with information on the CE results based on a model fitted
using the function `selection`, `pattern` or `hurdle`. Summary
information on the main parameters of interests is provided.

## References

Baio, G.(2012). *Bayesian Methods in Health Economcis*. CRC/Chapman
Hall, London.

## See also

[`selection`](https://angabrio.github.io/missingHE/reference/selection.md)
[`pattern`](https://angabrio.github.io/missingHE/reference/pattern.md)
[`hurdle`](https://angabrio.github.io/missingHE/reference/hurdle.md)
[`lmdm`](https://angabrio.github.io/missingHE/reference/lmdm.md)
[`diagnostic`](https://angabrio.github.io/missingHE/reference/diagnostic.md)
[`plot.missingHE`](https://angabrio.github.io/missingHE/reference/plot.missingHE.md)

## Author

Andrea Gabrio

## Examples

``` r
# For examples see the function selection, pattern or hurdle
#
#
```
