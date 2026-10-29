# A function to change the hyperprior parameters in a hurdle model

This function modifies default hyper prior parameter values in the type
of hurdle model selected according to the type of structural value
mechanism and distributions for the outcomes assumed.

## Usage

``` r
prior_hurdle(type, dist_e, dist_c, model_txt_info, model_string_jags)
```

## Arguments

- type:

  Type of structural value mechanism assumed. Choices are Structural
  Completely At Random (`SCAR`) and Structural At Random (`SAR`). For a
  complete list of all available hyper parameters and types of models
  see the manual.

- dist_e:

  distribution assumed for the effects. Current available choices are:
  Normal (`norm`), Beta (`beta`), Gamma (`gamma`), Exponential (`exp`),
  Weibull (`weib`), Logistic (`logis`), Poisson (`pois`), Negative
  Binomial (`negbin`) or Bernoulli (`bern`)

- dist_c:

  Distribution assumed for the costs. Current available choices are:
  Normal (`norm`), Gamma (`gamma`) or LogNormal (`lnorm`)

- model_txt_info:

  list containing model specification information used to write the
  `txt` file of the `JAGS` model.

- model_string_jags:

  text file of the model.

## Examples

``` r
#Internal function only
#no examples
#
#
```
