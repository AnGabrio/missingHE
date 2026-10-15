# Postsum method for inspecting summary posterior estimates of a Bayesian model fitted in JAGS

Prints the summary table for the model fitted, with the estimate of the
parameters and/or missing values.

## Usage

``` r
postsum(x, ...)

# S3 method for class 'missingHE'
postsum(x, display = "fixed", digits = 3, ...)
```

## Arguments

- x:

  A `missingHE` object containing the results of the Bayesian model run
  using the function
  [`selection`](https://angabrio.github.io/missingHE/reference/selection.md),
  [`pattern`](https://angabrio.github.io/missingHE/reference/pattern.md),
  [`hurdle`](https://angabrio.github.io/missingHE/reference/hurdle.md)
  or [`lmdm`](https://angabrio.github.io/missingHE/reference/lmdm.md).

- ...:

  additional arguments affecting the output produced. For example:
  `param = ` JAGS name assigned to the specific model parameters. to be
  shown in the printed table.

- display:

  Set of parameters for which posterior summaries are displayed. Choices
  are: fixed effects ('fixed'), random effects ('random'), individual
  log likelihood values ('loglik'), individual conditional parameter
  values ('conditional') and individual imputed values ('mis').

- digits:

  Integer indicating the number of decimal places to be used for
  rounding (default = 3).

## Author

Andrea Gabrio
