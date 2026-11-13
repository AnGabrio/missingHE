# Bayesian Hurdle Models

Full Bayesian cost-effectiveness models to handle missing data in the
outcomes using Hurdle models under a variety of alternative parametric
distributions for the effect and cost variables. Alternative assumptions
about the mechanisms of the structural values are implemented using a
hurdle approach. The analysis is performed using the `BUGS` language,
which is implemented in the software `JAGS` using the function
[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html). The output
is stored in an object of class `missingHE`.

## Usage

``` r
hurdle(
  data,
  model.eff,
  model.cost,
  model.se = se ~ 1,
  model.sc = sc ~ 1,
  se = 1,
  sc = 0,
  dist_e,
  dist_c,
  type,
  prob = c(0.025, 0.975),
  n.chains = 2,
  n.iter = 10000,
  n.burnin = floor(n.iter/2),
  inits = NULL,
  n.thin = 1,
  save_model = FALSE,
  prior = "default",
  center = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame in which to find the variables supplied in `model.eff`,
  `model.cost` (model formulas for effects and costs) and `model.se`,
  `model.sc` (model formulas for the structural effect and cost models).
  Among these, effectiveness, cost and treatment indicator variables
  must always be provided and named `e`, `c` and `trt`, respectively.

- model.eff:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be a health economic effectiveness outcome (`e`) whose
  name must correspond to that used in `data`. Any covariates in the
  model must be provided on the right-hand side of the formula. If there
  are no covariates, `1` should be specified on the right hand side of
  the formula. By default, covariates are placed on the "location"
  parameter of the distribution through a linear model. Random effects
  can also be specified for each model parameter. See details for how
  these can be specified.

- model.cost:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be a health economic cost outcome (`c`) whose name must
  correspond to that used in `data`. Any covariates in the model must be
  provided on the right-hand side of the formula. If there are no
  covariates, `1` should be specified on the right hand side of the
  formula. By default, covariates are placed on the "location" parameter
  of the distribution through a linear model. A joint bivariate
  distribution for effects and costs can be specified by including `e`
  on the right-hand side of the formula for the costs model. Random
  effects can also be specified for each model parameter. See details
  for how these can be specified.

- model.se:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be indicated with the term `se` (structural effects).
  Any covariates in the model must be provided on the right-hand side of
  the formula. If there are no covariates, `1` should be specified on
  the right hand side of the formula. By default, covariates are placed
  on the "probability" parameter for the structural effects through a
  logistic-linear model. Random effects can also be specified for each
  model parameter. See details for how these can be specified.

- model.sc:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be indicated with the term `sc` (structural costs). Any
  covariates in the model must be provided on the right-hand side of the
  formula. If there are no covariates, `1` should be specified on the
  right hand side of the formula. By default, covariates are placed on
  the "probability" parameter for the structural costs through a
  logistic-linear model. Random effects can also be specified for each
  model parameter. See details for how these can be specified.

- se:

  Structural value to be found in the effect variables defined in
  `data`. If set to `NULL`, no structural value is chosen and a standard
  model for the effects is run.

- sc:

  Structural value to be found in the cost variables defined in `data`.
  If set to `NULL`, no structural value is chosen and a standard model
  for the costs is run.

- dist_e:

  Distribution assumed for the effects. Current available choices are:
  Normal (`norm`), Beta (`beta`), Gamma (`gamma`), Exponential (`exp`),
  Weibull (`weib`), Logistic (`logis`), Poisson (`pois`), Negative
  Binomial (`negbin`) or Bernoulli (`bern`).

- dist_c:

  Distribution assumed for the costs. Current available choices are:
  Normal (`norm`), Gamma (`gamma`) or LogNormal (`lnorm`).

- type:

  Type of structural value mechanism assumed. Choices are Structural
  Completely At Random (`SCAR`), and Structural At Random (`SAR`).

- prob:

  A numeric vector of probabilities within the range `(0, 1)`,
  representing the upper and lower CI sample quantiles to be calculated
  and returned for the imputed values.

- n.chains:

  Number of chains.

- n.iter:

  Number of iterations.

- n.burnin:

  Number of warmup iterations.

- inits:

  A list with elements equal to the number of chains selected; each
  element of the list is itself a list of starting values for the `JAGS`
  model, or a function creating (possibly random) initial values. If
  `inits` is `NULL`, `JAGS` will generate initial values for all the
  model parameters.

- n.thin:

  Thinning interval.

- save_model:

  Logical. If `save_model` is `TRUE` a `txt` file containing the model
  code is printed in the current working directory.

- prior:

  A list containing the hyperprior values provided by the user. Each
  element of this list must be a vector containing the user-provided
  prior distribution and parameter values and must be named with the
  name of the corresponding parameter. For example, the hyperprior
  values for the standard deviation parameter for the effects can be
  provided using the list
  `prior = list("sigma.prior.e" = c("unif", 0, 10))`. For more
  information about how to provide prior hypervalues for different types
  of parameters and models see details. If `prior` is set to `default`,
  the default values will be used.

- center:

  Logical. If `center` is `TRUE`, all covariates in both the effect and
  cost models will be centred.

- ...:

  Additional arguments that can be provided by the user. Examples are
  the additional arguments that can be provided to the function
  [`BCEA::bcea()`](https://n8thangreen.github.io/BCEA/reference/bcea.html)
  to summarise the health economic evaluation results.

## Value

An object of the class `missingHE` containing the following elements:

- **data_set**: A list containing the original data set provided in
  `data` (see Arguments). Additional information is also included about,
  among others, the number of observed and missing individuals, the
  total number of individuals by treatment arm and the indicator vectors
  for the missing values.

- **model_output**: A list containing the output of a `JAGS` model
  generated from the functions
  [`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html), and the
  posterior samples for the main parameters of the model.

- **cea**: A list containing the output of the economic evaluation
  performed using the function
  [`BCEA::bcea()`](https://n8thangreen.github.io/BCEA/reference/bcea.html).

- **type**: A character variable that indicate which type of missing
  value mechanism used to run the model, either `MAR` or `MNAR` (see
  details).

- **data_format**: A character variable that indicates which type of
  analysis was conducted, either using a `wide` or `long` dataset.

- **se_prior**: A numeric value that indicates the value assumed for the
  standard deviation associated with the degenerate structural effects
  prior distribution to implement the model in `JAGS` (deafult =
  0.00001).

- **sc_prior**: A numeric value that indicates the value assumed for the
  standard deviation associated with the degenerate structural costs
  prior distribution to implement the model in `JAGS` (deafult =
  0.00001).

- **call**: The matched call.

- **prior**: The custom prior choices specified by the user. If no
  custom prior specified, default priors are used (see details).

## Details

Depending on the distributions specified for the outcome variables in
the arguments `dist_e` and `dist_c` and the type of structural value
mechanism specified in the argument `type`, different hurdle models are
built and run in the background by the function `hurdle`. These are
mixture models defined by two components: the first one is a mass
distribution at the spike, while the second is a parametric model
applied to the natural range of the relevant variable. Usually, a
logistic regression is used to estimate the probability of incurring a
"structural" value (e.g. `0` for the costs, or `1` for the effects);
this is then used to weight the mean of the "non-structural" values
estimated in the second component. A simple example can be used to show
how hurdle models are specified. Consider a data set comprising a
response variable \\y\\ and a set of centered covariates \\X_j\\, for
\\i = j, ..., J\\.Specifically, for each subject in the trial \\i = 1,
..., n\\ we define an indicator variable \\s_i\\ taking value `1` if the
\\i\\-th individual is associated with a structural value and `0`
otherwise. This is modelled as:

\$\$s_i ~ Bernoulli(\pi_i)\$\$ \$\$logit(\pi_i) = \sum\gamma_j X_j\$\$
where

- \\\pi_i\\ is the individual probability of a structural value in
  \\y\\.

- \\\gamma_j\\ represents the impact on the probability of a structural
  value in \\y\\ of the covariate \\X_j\\.

When \\\gamma_j = 0\\, the model assumes a `SCAR` mechanism, while when
\\\gamma_j != 0\\ the mechanism is `SAR`. For the parameters indexing
the structural value model, the default prior distributions assumed are:

- \\\gamma_j ~ Normal(0, 0.01)\\.

When user-defined hyperprior values are supplied via the argument
`prior` in the function `hurdle`, the elements of this list (see
Arguments) must be vectors containing the user-provided distribution
name and hyperprior values and must take specific names according to the
parameters they are associated with. Specifically, the names accepted by
**missingHE** are the following:

- auxiliary parameters \\\sigma\\: `sigma.prior.e` (effects) and/or
  `sigma.prior.c` (costs).

- covariate parameters \\\alpha_j\\ and \\\beta_j\\: `alpha.prior`
  (effects) and/or `beta.prior` (costs).

- covariate parameters in the structural values model \\\gamma_j\\ (if
  covariate data provided): `gamma.prior.e` (effects) and/or
  `gamma.prior.c` (costs).

For simplicity, here we assumed that the set of covariates \\X_j\\ used
in the models for the effects/costs and in the model of the structural
effect/cost values is the same. However, it is possible to specify
different sets of covariates for each model using the arguments in the
function `hurdle` (see Arguments).

For each model, random effects can also be specified for each parameter
by adding the term `+ (x | z)` to each model formula, where `x` is the
fixed regression coefficient for which also the random effects are
desired and `z` is the clustering variable across which the random
effects are specified (must be the name of a factor variable in the
dataset).

For an overview of the general modeling framework and additional package
utilities, see
[missingHE-package](https://angabrio.github.io/missingHE/reference/missingHE-package.md).

## References

Gabrio, A., Mason, AJ. Baio, G. (2019). *A full Bayesian model to handle
structural ones and missingness in economic evaluations from
individual-level data*. Statistics in medicine, 38(8), 1399-1420
\<10.1002/sim.8045\>

Mason, AJ. Gomes, M. Carpenter, J. Grieve, R. (2021). *Flexible Bayesian
longitudinal models for cost‐effectiveness analyses with informative
missing data*. Health economics, 30(12), 3138-3158. \<10.1002/hec.4408\>

Daniels, MJ. Hogan, JW. *Missing Data in Longitudinal Studies:
strategies for Bayesian modelling and sensitivity analysis*, CRC/Chapman
Hall. <doi:10.1201/9781420011180>

Baio, G.(2012). *Bayesian Methods in Health Economics*. CRC/Chapman
Hall, London. <doi:10.1201/b13099>

Plummer, M. *JAGS: A program for analysis of Bayesian graphical models
using Gibbs sampling*. (2003).
<https://www.r-project.org/conferences/DSC-2003/Proceedings/Plummer.pdf>

## See also

[missingHE-package](https://angabrio.github.io/missingHE/reference/missingHE-package.md)
[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`BCEA::bcea()`](https://n8thangreen.github.io/BCEA/reference/bcea.html).

## Author

Andrea Gabrio

## Examples

``` r
if (isTRUE(tryCatch(requireNamespace("rjags", quietly = TRUE), error = function(e) FALSE))) {
# Quick example to run using the MenSS dataset

# Run the model using the hurdle function assuming a SCAR mechanism
# Use subset of dataset to produce quick results
MenSS_ex <- MenSS[MenSS$age > 25, ]

# Use only 50 iterations to run a quick check
model.hurdle <- hurdle(data = MenSS_ex, model.eff = e ~ trt, model.cost = c ~ trt,
   model.se = se ~ 1, model.sc = sc ~ 1, se = NULL, sc = 0, dist_e = "norm", dist_c = "gamma",
   type = "SCAR", n.chains = 2, n.iter = 50)

# Print the results of the JAGS model
print(model.hurdle)
#

# Use dic information criterion to assess model fit
pic.dic <- pic(model.hurdle, criterion = "dic", cases = "cc")
pic.dic
#
# Extract regression coefficient estimates
coef(model.hurdle)
#

DONTSHOW({
# Use waic information criterion to assess model fit
pic.waic <- pic(model.hurdle, criterion = "waic", cases = "cc")
pic.waic
})

# Assess model convergence using graphical tools
# Produce histograms of the posterior samples for the mean effects
diag.hist <- diagnostic(model.hurdle, type = "histogram", param = "mu.e")
#

# Compare observed effect data with imputations from the model
# using plots (posterior means and credible intervals)
p1 <- plot(model.hurdle, class = "scatter", outcome = "effects")
#

# Summarise the CEA information from the model
summary(model.hurdle)

# \donttest{
# Further examples which take longer to run
model.hurdle <- hurdle(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt + e,
   model.se = se ~ age, model.sc = sc ~ age, se = 1, sc = 0, dist_e = "beta", dist_c = "gamma",
   type = "SAR", n.chains = 2, n.iter = 500)
#
# Print results for all imputed values
print(model.hurdle)

# Use looic to assess model fit
pic.looic <- pic(model.hurdle, criterion = "looic", cases = "cc")
pic.looic

# Show density plots for mean costs parameters
diag.den <- diagnostic(model.hurdle, type = "denplot", param = "mu.c")

# Plots of imputations for all data
p1 <- plot(model.hurdle, class = "scatter", outcome = "both")

# Summarise the CEA results
summary(model.hurdle)

# }
}
#> module glm loaded
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 87
#>    Unobserved stochastic nodes: 220
#>    Total graph size: 1884
#> 
#> Initializing model
#> 
#> Registered S3 method overwritten by 'BCEA':
#>   method     from
#>   plot.evppi voi 
#> Call:
#> hurdle(model.eff = e ~ trt, model.cost = c ~ trt, model.se = se ~ 
#>     1, model.sc = sc ~ 1, data = MenSS_ex, se = NULL, sc = 0, 
#>     dist_e = "norm", dist_c = "gamma", type = "SCAR", n.chains = 2, 
#>     n.iter = 50)
#> 
#> Coefficients:
#> $model.eff
#> (Intercept)        trt2 
#>       0.885      -0.004 
#> 
#> $model.cost
#> (Intercept)        trt2 
#>       6.835      -0.100 
#> 
#> $model.se
#> [1] NA
#> 
#> $model.sc
#> (Intercept) 
#>      -0.945 
#> 
#> $Effects
#>               Mean    SD     QL    QU
#> (Intercept)  0.885 0.036  0.822 0.950
#> trt2        -0.004 0.043 -0.064 0.085
#> 
#> $Costs
#>               Mean    SD     QL    QU
#> (Intercept)  6.835 0.980  5.567 8.065
#> trt2        -0.100 0.237 -0.463 0.326
#> 
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under SCAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.885 0.034 0.816 0.966
#> 2 0.883 0.038 0.833 0.969
#> 
#>  
#>  Mean costs by intervention 
#>   Mean  SD    QL   QU
#> 1  712 582 123.6 1739
#> 2  655 561  92.9 1809
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 43530 1609 40577 46986
#> 2 43510 1972 39934 47496
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 184
#>    Unobserved stochastic nodes: 463
#>    Total graph size: 6241
#> 
#> Initializing model
#> 
#> Call:
#> hurdle(model.eff = e ~ trt, model.cost = c ~ trt + e, model.se = se ~ 
#>     age, model.sc = sc ~ age, data = MenSS, se = 1, sc = 0, dist_e = "beta", 
#>     dist_c = "gamma", type = "SAR", n.chains = 2, n.iter = 500)
#> 
#> Coefficients:
#> $model.eff
#> (Intercept)        trt2 
#>       1.779      -0.166 
#> 
#> $model.cost
#> (Intercept)        trt2           e 
#>       6.788       0.093      -0.583 
#> 
#> $model.se
#> (Intercept)         age 
#>       0.207      -0.026 
#> 
#> $model.sc
#> (Intercept)         age 
#>      -1.200       0.004 
#> 
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under SAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.940 0.014 0.911 0.965
#> 2 0.934 0.018 0.896 0.964
#> 
#>  
#>  Mean costs by intervention 
#>   Mean  SD    QL   QU
#> 1  731 532  87.3 1675
#> 2  761 546 110.1 1765
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 46254  923 44378 47969
#> 2 45926 1111 43440 47739
#
```
