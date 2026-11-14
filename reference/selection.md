# Bayesian Selection Models

Full Bayesian cost-effectiveness models to handle missing data in the
outcomes under different missing data mechanism assumptions, using
alternative parametric distributions for the effect and cost variables
and using a selection model approach to identify the model. The analysis
is performed using the `BUGS` language, which is implemented in the
software `JAGS` using the function
[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html) The output
is stored in an object of class `missingHE`.

## Usage

``` r
selection(
  data,
  model.eff,
  model.cost,
  model.me = me ~ 1,
  model.mc = mc ~ 1,
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
  `model.cost` (model formulas for effects and costs) and `model.me`,
  `model.mc` (model formulas for the missing effect and cost models).
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

- model.me:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be indicated with the term `me` (missing effects) and
  any covariates must be provided on the right-hand side of the formula.
  If there are no covariates, `1` should be specified on the right hand
  side of the formula. By default, covariates are placed on the
  "probability" parameter for the missing effects through a
  logistic-linear model. Random effects can also be specified for each
  model parameter. See details for how these can be specified.

- model.mc:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be indicated with the term `mc` (missing costs) and any
  covariates must be provided on the right-hand side of the formula. If
  there are no covariates, `1` should be specified on the right hand
  side of the formula. By default, covariates are placed on the
  "probability" parameter for the missing costs through a
  logistic-linear model. Random effects can also be specified for each
  model parameter. See details for how these can be specified.

- dist_e:

  Distribution assumed for the effects. Current available choices are:
  Normal (`norm`), Beta (`beta`), Gamma (`gamma`), Exponential (`exp`),
  Weibull (`weib`), Logistic (`logis`), Poisson (`pois`), Negative
  Binomial (`negbin`) or Bernoulli (`bern`).

- dist_c:

  Distribution assumed for the costs. Current available choices are:
  Normal (`norm`), Gamma (`gamma`) or LogNormal (`lnorm`).

- type:

  Type of missingness mechanism assumed. Choices are Missing At Random
  (`MAR`) and Missing Not At Random (`MNAR`).

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

  Logical. If `save_model` is `TRUE`, a `txt` file containing the model
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

- **call**: The matched call.

- **prior**: The custom prior choices specified by the user. If no
  custom prior specified, default priors are used (see details).

## Details

Depending on the distributions specified for the outcome variables in
the arguments `dist_e` and `dist_c` and the type of missingness
mechanism specified in the argument `type`, different selection models
are built and run in the background by the function `selection`. These
models consist in logistic regressions to estimate the probability of
missingness in one or both outcomes. A simple example can be used to
show how selection models are specified. Consider a data set comprising
a response variable \\y\\ and a set of covariates \\X_j\\, for \\i = j,
..., J\\. For each subject in the trial \\i = 1, ..., n\\ we define an
indicator variable \\m_i\\ taking value `1` if the \\i\\-th individual
is associated with a missing value and `0` otherwise. This is modelled
as:

\$\$m_i ~ Bernoulli(\pi_i)\$\$ \$\$logit(\pi_i) = \sum\gamma_j X_j +
\delta y_i\$\$ where

- \\\pi_i\\ is the individual probability of a missing value in \\y\\.

- \\\gamma_j\\ represents the impact on the probability of a missing
  value in \\y\\ of the covariate \\X_j\\.

- \\\delta\\ represents the impact on the probability of a missing value
  in \\y\\ of the outcome values.

When \\\delta = 0\\ the model assumes a `MAR` mechanism, while when
\\\delta != 0\\ the mechanism is `MNAR`. For the parameters indexing the
missingness model, the default prior distributions assumed are:

- \\\gamma_j ~ Normal(0, 0.01)\\.

- \\\delta ~ Normal(0, 1)\\.

When user-defined hyperprior values are supplied via the argument
`prior` in the function `selection`, the elements of this list (see
Arguments) must be vectors containing the user-provided distribution
name and hyperprior values and must take specific names according to the
parameters they are associated with. Specifically, the names for the
parameters indexing the model which are accepted by **missingHE** are
the following:

- auxiliary parameters \\\sigma\\: `sigma.prior.e` (effects) and/or
  `sigma.prior.c` (costs).

- covariate parameters \\\alpha_j\\ and \\\beta_j\\: `alpha.prior`
  (effects) and/or `beta.prior` (costs).

- covariate parameters in the missingness model \\\gamma_j\\ (if
  covariate data provided): `gamma.prior.e` (effects) and/or
  `gamma.prior.c` (costs).

- mnar parameter \\\delta\\: `delta.prior.e` (effects) and/or
  `delta.prior.c` (costs).

For simplicity, here we have assumed that the set of covariates \\X_j\\
used in the models for the effects/costs and in the model of the missing
effect/cost values is the same. However, it is possible to specify
different sets of covariates for each model using the arguments in the
function `selection` (see Arguments).

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
# Quick example to run using subset of MenSS dataset
MenSS.subset <- MenSS[50:100, ]

# Run the model using the selection function assuming a SCAR mechanism
# Use only 100 iterations to run a quick check
model.selection <- selection(data = MenSS.subset, model.eff = e ~ trt, model.cost = c ~ trt,
   model.me = me ~ 1, model.mc = mc ~ 1, dist_e = "norm", dist_c = "norm",
   type = "MAR", n.chains = 2, n.iter = 100)

# Print the results of the JAGS model
print(model.selection)
#

# Use dic information criterion to assess model fit
pic.dic <- pic(model.selection, criterion = "dic", cases = "cc")
pic.dic
#

# Extract regression coefficient estimates
coef(model.selection)
#

DONTSHOW({
# Use waic to assess model fit
pic.waic <- pic(model.selection, criterion = "waic", cases = "cc")
pic.waic
})

# Assess model convergence using graphical tools
# Produce histograms of the posterior samples for the mean effects
diag.hist <- diagnostic(model.selection, type = "histogram", param = "mu.e")
#

# Compare observed effect data with imputations from the model
# using plots (posterior means and credible intervals)
p1 <- plot(model.selection, class = "scatter", outcome = "effects")
#

# Summarise the CEA information from the model
summary(model.selection)

# \donttest{
# Further examples which take longer to run
model.selection <- selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt + e,
   model.se = me ~ age, model.mc = mc ~ age, dist_e = "norm", dist_c = "norm",
   type = "MAR", n.chains = 2, n.iter = 500)
#
# Print results for all imputed values
print(model.selection)

# Use looic to assess model fit
pic.looic <- pic(model.selection, criterion = "looic", cases = "cc")
pic.looic

# Show density plots for mean costs parameters
diag.den <- diagnostic(model.selection, type = "denplot", param = "mu.c")

# Plots of imputations for all data
p1 <- plot(model.selection, class = "scatter", outcome = "both")

# Summarise the CEA results
summary(model.selection)

# }
}
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Fully observed stochastic nodes: 132
#>    Unobserved stochastic nodes: 80
#>    Total graph size: 755
#> 
#> Initializing model
#> 
#> Call:
#> selection(model.eff = e ~ trt, model.cost = c ~ trt, model.me = me ~ 
#>     1, model.mc = mc ~ 1, data = MenSS.subset, dist_e = "norm", 
#>     dist_c = "norm", type = "MAR", n.chains = 2, n.iter = 100)
#> 
#> Coefficients:
#> $model.eff
#> (Intercept)        trt2 
#>       0.908       0.060 
#> 
#> $model.cost
#> (Intercept)        trt2 
#>      48.852      15.051 
#> 
#> $model.me
#> (Intercept) 
#>       0.876 
#> 
#> $model.mc
#> (Intercept) 
#>       0.887 
#> 
#> $Effects
#>              Mean    SD     QL    QU
#> (Intercept) 0.908 0.031  0.856 0.964
#> trt2        0.060 0.061 -0.060 0.166
#> 
#> $Costs
#>               Mean     SD      QL      QU
#> (Intercept) 48.852 35.003 -17.606 107.821
#> trt2        15.051 33.737 -49.552  85.196
#> 
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under MAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.908 0.032 0.854 0.967
#> 2 0.967 0.050 0.867 1.060
#> 
#>  
#>  Mean costs by intervention 
#>   Mean   SD    QL  QU
#> 1 52.3 41.2 -27.4 126
#> 2 63.3 48.4 -34.5 158
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 45343 1612 42685 48318
#> 2 48266 2481 43354 52994
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Fully observed stochastic nodes: 410
#>    Unobserved stochastic nodes: 236
#>    Total graph size: 3109
#> 
#> Initializing model
#> 
#>   |                                                          |                                                  |   0%  |                                                          |**                                                |   4%  |                                                          |****                                              |   8%  |                                                          |******                                            |  12%  |                                                          |********                                          |  16%  |                                                          |**********                                        |  20%  |                                                          |************                                      |  24%  |                                                          |**************                                    |  28%  |                                                          |****************                                  |  32%  |                                                          |******************                                |  36%  |                                                          |********************                              |  40%  |                                                          |**********************                            |  44%  |                                                          |************************                          |  48%  |                                                          |**************************                        |  52%  |                                                          |****************************                      |  56%  |                                                          |******************************                    |  60%  |                                                          |********************************                  |  64%  |                                                          |**********************************                |  68%  |                                                          |************************************              |  72%  |                                                          |**************************************            |  76%  |                                                          |****************************************          |  80%  |                                                          |******************************************        |  84%  |                                                          |********************************************      |  88%  |                                                          |**********************************************    |  92%  |                                                          |************************************************  |  96%  |                                                          |**************************************************| 100%
#>   |                                                          |                                                  |   0%  |                                                          |**                                                |   4%  |                                                          |****                                              |   8%  |                                                          |******                                            |  12%  |                                                          |********                                          |  16%  |                                                          |**********                                        |  20%  |                                                          |************                                      |  24%  |                                                          |**************                                    |  28%  |                                                          |****************                                  |  32%  |                                                          |******************                                |  36%  |                                                          |********************                              |  40%  |                                                          |**********************                            |  44%  |                                                          |************************                          |  48%  |                                                          |**************************                        |  52%  |                                                          |****************************                      |  56%  |                                                          |******************************                    |  60%  |                                                          |********************************                  |  64%  |                                                          |**********************************                |  68%  |                                                          |************************************              |  72%  |                                                          |**************************************            |  76%  |                                                          |****************************************          |  80%  |                                                          |******************************************        |  84%  |                                                          |********************************************      |  88%  |                                                          |**********************************************    |  92%  |                                                          |************************************************  |  96%  |                                                          |**************************************************| 100%
#> Call:
#> selection(model.eff = e ~ trt, model.cost = c ~ trt + e, model.mc = mc ~ 
#>     age, data = MenSS, dist_e = "norm", dist_c = "norm", type = "MAR", 
#>     n.chains = 2, n.iter = 500)
#> 
#> Coefficients:
#> $model.eff
#> (Intercept)        trt2 
#>       0.904      -0.001 
#> 
#> $model.cost
#> (Intercept)        trt2           e 
#>      83.545      24.824      -5.340 
#> 
#> $model.me
#> (Intercept) 
#>       0.893 
#> 
#> $model.mc
#> (Intercept)         age 
#>       1.415      -0.017 
#> 
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under MAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.905 0.024 0.857 0.950
#> 2 0.902 0.027 0.849 0.957
#> 
#>  
#>  Mean costs by intervention 
#>   Mean   SD   QL  QU
#> 1   84 26.6 35.2 137
#> 2  108 36.5 38.1 176
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 45151 1176 42756 47434
#> 2 45013 1371 42337 47760
#
```
