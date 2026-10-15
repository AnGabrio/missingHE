# A function to read and re-arrange the data in different ways for the hurdle model

This internal function imports the data and outputs only those variables
that are needed to run the hurdle model according to the information
provided by the user.

## Usage

``` r
data_read_hurdle(
  data,
  model.eff,
  model.cost,
  model.se,
  model.sc,
  se,
  sc,
  cov_matrix,
  type,
  center,
  fixed_e,
  fixed_c,
  random_e,
  random_c,
  trt_lev,
  trt_pos
)
```

## Arguments

- data:

  A data frame in which to find variables supplied in `model.eff`,
  `model.cost` (model formulas for effects and costs) and `model.se`,
  `model.sc` (model formulas for the structural effect and cost models)
  . Among these, effectiveness, cost and treatment indicator variables
  must always be provided and named 'e', 'c' and 'trt' respectively.

- model.eff:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be a health economics effectiveness outcome ('e') whose
  name must correspond to that used in `data`, and any covariates are
  given on the right-hand side. If there are no covariates, specify `1`
  on the right hand side. Random effects can also be specified for each
  model parameter.

- model.cost:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be a health economics cost outcome ('c') whose name must
  correspond to that used in `data`, and any covariates are given on the
  right-hand side. If there are no covariates, specify `1` on the right
  hand side. By default, covariates are placed on the "location"
  parameter of the distribution through a linear model. Random effects
  can also be specified for each model parameter.

- model.se:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be a health economics effectiveness outcome ('e') whose
  name must correspond to that used in `data`, and any covariates used
  to estimate the probability of structural effects are given on the
  right-hand side. If there are no covariates, specify `1` on the right
  hand side. By default, covariates are placed on the "probability"
  parameter for the structural effects through a logistic-linear model.
  Random effects can also be specified for each model parameter.

- model.sc:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be a health economics cost outcome ('c') whose name must
  correspond to that used in `data`, and any covariates used to estimate
  the probability of structural costs are given on the right-hand side.
  If there are no covariates, specify `1` on the right hand side. By
  default, covariates are placed on the "probability" parameter for the
  structural costs through a logistic-linear model. Random effects can
  also be specified for each model parameter.

- se:

  Structural value to be found in the effect data defined in `data`. If
  set to `NULL`, no structural value is chosen and a standard model for
  the effects is run.

- sc:

  Structural value to be found in the cost data defined in `data`. If
  set to `NULL`, no structural value is chosen and a standard model for
  the costs is run.

- cov_matrix:

  Covariance matrix containing all model covariates.

- type:

  Type of structural value mechanism assumed, either 'SCAR' (Structural
  Completely At Random) or 'SAR' (Strcutural At Random).

- center:

  Logical. If `center` is `TRUE` all the covariates in the model are
  centered.

- fixed_e:

  Formula of the fixed effects specified in the effectiveness model.

- fixed_c:

  Formula of the fixed effects specified in the cost model.

- random_e:

  Formula of the random effects specified in the effectiveness model.

- random_c:

  Formula of the random effects specified in the cost model.

- trt_lev:

  Levels of the treatment indicator factor.

- trt_pos:

  Numeric position of the treatment indicator within the model formula.

## Examples

``` r
#Internal function only
#no examples
#
#
```
