# A function to execute a selection model and get posterior results

This function fits a `JAGS` using the
[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html) function
and obtain posterior inferences.

## Usage

``` r
run_selection(data_model, type, dist_e, dist_c, model_info)
```

## Arguments

- data_model:

  list containing the data for the model to be passed to `JAGS`.

- type:

  Type of missingness mechanism assumed. Choices are Missing At Random
  (`MAR`), Missing Not At Random for the effects (`MNAR_eff`), Missing
  Not At Random for the costs (`MNAR_cost`), and Missing Not At Random
  for both (`MNAR`).

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
