# A function to execute a hurdle model and get posterior results

This function fits a `JAGS` using the
[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html) function
and obtain posterior inferences.

## Usage

``` r
run_hurdle(data_model, type, dist_e, dist_c, model_info)
```

## Arguments

- data_model:

  list containing the data for the model to be passed to `JAGS`.

- type:

  Type of structural value mechanism assumed. Choices are Structural
  Completely At Random (`SCAR`), Structural At Random (`SAR`),

- dist_e:

  distribution assumed for the effects. Current available choices are:
  Normal (`norm`), Beta (`beta`), Gamma (`gamma`), Exponential (`exp`),
  Weibull (`weib`), Logistic (`logis`), Poisson (`pois`), Negative
  Binomial (`negbin`) or Bernoulli (`bern`).

- dist_c:

  Distribution assumed for the costs. Current available choices are:
  Normal (`norm`), Gamma (`gamma`) or LogNormal (`lnorm`).

- model_info:

  list containing model and MCMC information to be passed to `JAGS`.

## Examples

``` r
#Internal function only
#No examples
#
#
```
