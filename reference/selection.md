# Bayesian Selection Models

Full Bayesian cost-effectiveness models to handle missing data in the
outcomes under different missing data mechanism assumptions, using
alternative parametric distributions for the effect and cost variables
and using a selection model approach to identify the model. The analysis
is performed using the `BUGS` language, which is implemented in the
software `JAGS` using the function
[`jags`](https://rdrr.io/pkg/R2jags/man/jags.html) The output is stored
in an object of class 'missingHE'.

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
  Among these, effectiveness, cost and treatment indicator (only two
  arms) variables must always be provided and named 'e', 'c' and 'trt',
  respectively.

- model.eff:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be a health economic effectiveness outcome ('e') whose
  name must correspond to that used in `data`. Any covariates in the
  model must be provided on the right-hand side of the formula. If there
  are no covariates, `1` should be specified on the right hand side of
  the formula. By default, covariates are placed on the "location"
  parameter of the distribution through a linear model. Random effects
  can also be specified for each model parameter. See details for how
  these can be specified.

- model.cost:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be a health economic cost outcome ('c') whose name must
  correspond to that used in `data`. Any covariates in the model must be
  provided on the right-hand side of the formula. If there are no
  covariates, `1` should be specified on the right hand side of the
  formula. By default, covariates are placed on the "location" parameter
  of the distribution through a linear model. A joint bivariate
  distribution for effects and costs can be specified by including 'e'
  on the right-hand side of the formula for the costs model. Random
  effects can also be specified for each model parameter. See details
  for how these can be specified.

- model.me:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be indicated with the term 'me'(missing effects) and any
  covariates must be provided on the right-hand side of the formula. If
  there are no covariates, `1` should be specified on the right hand
  side of the formula. By default, covariates are placed on the
  "probability" parameter for the missing effects through a
  logistic-linear model. Random effects can also be specified for each
  model parameter. See details for how these can be specified.

- model.mc:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be indicated with the term 'mc'(missing costs) and any
  covariates must be provided on the right-hand side of the formula. If
  there are no covariates, `1` should be specified on the right hand
  side of the formula. By default, covariates are placed on the
  "probability" parameter for the missing costs through a
  logistic-linear model. Random effects can also be specified for each
  model parameter. See details for how these can be specified.

- dist_e:

  Distribution assumed for the effects. Current available choices are:
  Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
  Weibull ('weib'), Logistic ('logis'), Poisson ('pois'), Negative
  Binomial ('negbin') or Bernoulli ('bern').

- dist_c:

  Distribution assumed for the costs. Current available choices are:
  Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').

- type:

  Type of missingness mechanism assumed. Choices are Missing At Random
  (MAR) and Missing Not At Random (MNAR).

- prob:

  A numeric vector of probabilities within the range (0, 1),
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
  `prior = list('sigma.prior.e' = c("unif", 0, 10))`. For more
  information about how to provide prior hypervalues for different types
  of parameters and models see details. If `prior` is set to 'default',
  the default values will be used.

- center:

  Logical. If `center` is `TRUE`, all covariates in both the effect and
  cost models will be centred.

- ...:

  Additional arguments that can be provided by the user. Examples are
  the additional arguments that can be provided to the function
  [`bcea`](https://n8thangreen.github.io/BCEA/reference/bcea.html) to
  summarise the health economic evaluation results.

## Value

An object of the class 'missingHE' containing the following elements

- data_set:

  A list containing the original data set provided in `data` (see
  Arguments). Additional information is also included about, among
  others, the number of observed and missing individuals, the total
  number of individuals by treatment arm and the indicator vectors for
  the missing values

- model_output:

  A list containing the output of a `JAGS` model generated from the
  functions [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html), and the
  posterior samples for the main parameters of the model

- cea:

  A list containing the output of the economic evaluation performed
  using the function
  [`bcea`](https://n8thangreen.github.io/BCEA/reference/bcea.html)

- type:

  A character variable that indicate which type of missing value
  mechanism used to run the model, either `MAR` or `MNAR` (see details)

- data_format:

  A character variable that indicates which type of analysis was
  conducted, either using a `wide` or `long` dataset

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
as: \$\$m_i ~ Bernoulli(\pi_i)\$\$ \$\$logit(\pi_i) = \sum\gamma_j X_j +
\delta y_i\$\$ where

- \\\pi_i\\ is the individual probability of a missing value in \\y\\.

- \\\gamma_j\\ represents the impact on the probability of a missing
  value in \\y\\ of the covariate \\X_j\\.

- \\\delta\\ represents the impact on the probability of a missing value
  in \\y\\ of the outcome values.

When \\\delta = 0\\ the model assumes a 'MAR' mechanism, while when
\\\delta != 0\\ the mechanism is 'MNAR'. For the parameters indexing the
missingness model, the default prior distributions assumed are:

- \\\gamma_j ~ Normal(0, 0.01)\\

- \\\delta ~ Normal(0, 1)\\

When user-defined hyperprior values are supplied via the argument
`prior` in the function `selection`, the elements of this list (see
Arguments) must be vectors containing the user-provided distribution
name and hyperprior values and must take specific names according to the
parameters they are associated with. Specifically, the names for the
parameters indexing the model which are accepted by **missingHE** are
the following:

- auxiliary parameters \\\sigma\\: "sigma.prior.e"(effects) and/or
  "sigma.prior.c"(costs)

- covariate parameters \\\alpha_j\\ and \\\beta_j\\:
  "alpha.prior"(effects) and/or "beta.prior"(costs)

- covariate parameters in the missingness model \\\gamma_j\\ (if
  covariate data provided): "gamma.prior.e"(effects) and/or
  "gamma.prior.c"(costs)

- mnar parameter \\\delta\\: "delta.prior.e"(effects) and/or
  "delta.prior.c"(costs)

For simplicity, here we have assumed that the set of covariates \\X_j\\
used in the models for the effects/costs and in the model of the missing
effect/cost values is the same. However, it is possible to specify
different sets of covariates for each model using the arguments in the
function `selection` (see Arguments).

For each model, random effects can also be specified for each parameter
by adding the term + (x \| z) to each model formula, where x is the
fixed regression coefficient for which also the random effects are
desired and z is the clustering variable across which the random effects
are specified (must be the name of a factor variable in the dataset).
Multiple random effects can be specified using the notation + (x1 + x2
\| site) for each covariate that was included in the fixed effects
formula.

## References

Daniels, MJ. Hogan, JW. *Missing Data in Longitudinal Studies:
strategies for Bayesian modelling and sensitivity analysis*, CRC/Chapman
Hall.

Baio, G.(2012). *Bayesian Methods in Health Economics*. CRC/Chapman
Hall, London.

Gelman, A. Carlin, JB., Stern, HS. Rubin, DB.(2003). *Bayesian Data
Analysis, 2nd edition*, CRC Press.

Plummer, M. *JAGS: A program for analysis of Bayesian graphical models
using Gibbs sampling.* (2003).

## See also

[`jags`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`bcea`](https://n8thangreen.github.io/BCEA/reference/bcea.html)

## Author

Andrea Gabrio

## Examples

``` r
# Quick example to run using subset of MenSS dataset
MenSS.subset <- MenSS[50:100, ]

# Run the model using the selection function assuming a SCAR mechanism
# Use only 100 iterations to run a quick check
model.selection <- selection(data = MenSS.subset, model.eff = e ~ trt, model.cost = c ~ trt,
   model.me = me ~ 1, model.mc = mc ~ 1, dist_e = "norm", dist_c = "norm",
   type = "MAR", n.chains = 2, n.iter = 100)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 132
#>    Unobserved stochastic nodes: 80
#>    Total graph size: 755
#> 
#> Initializing model
#> 

# Print the results of the JAGS model
print(model.selection)
#> Call:
#> selection(model.eff = e ~ trt, model.cost = c ~ trt, model.me = me ~ 
#>     1, model.mc = mc ~ 1, data = MenSS.subset, dist_e = "norm", 
#>     dist_c = "norm", type = "MAR", n.chains = 2, n.iter = 100)
#> 
#> Coefficients:
#> $model.eff
#> (Intercept)        trt2 
#>       0.904       0.084 
#> 
#> $model.cost
#> (Intercept)        trt2 
#>      48.878      14.709 
#> 
#> $model.me
#> (Intercept) 
#>       0.861 
#> 
#> $model.mc
#> (Intercept) 
#>       0.889 
#> 
#

# Use dic information criterion to assess model fit
pic.dic <- pic(model.selection, criterion = "dic", cases = "cc")
pic.dic
#> $d_bar
#> [1] 249.4578
#> 
#> $pD
#> [1] 2.475565
#> 
#> $dic
#> [1] 251.9334
#> 
#> $d_hat
#> [1] 246.9823
#> 
#

# Extract regression coefficient estimates
coef(model.selection)
#> $Effects
#>              Mean    SD     QL    QU
#> (Intercept) 0.904 0.038  0.827 0.980
#> trt2        0.084 0.067 -0.043 0.201
#> 
#> $Costs
#>               Mean     SD      QL      QU
#> (Intercept) 48.878 31.705 -11.684 110.215
#> trt2        14.709 29.964 -37.039  64.457
#> 
#


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
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under MAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.904 0.041 0.824 0.985
#> 2 0.987 0.059 0.873 1.103
#> 
#>  
#>  Mean costs by intervention 
#>   Mean   SD     QL  QU
#> 1   47 37.9 -19.06 115
#> 2   67 44.0  -9.38 154
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 45161 2071 41149 49211
#> 2 49295 2947 43593 55130

# \donttest{
# Further examples which take longer to run
model.selection <- selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt + e,
   model.se = me ~ age, model.mc = mc ~ age, dist_e = "norm", dist_c = "norm",
   type = "MAR", n.chains = 2, n.iter = 500)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 410
#>    Unobserved stochastic nodes: 236
#>    Total graph size: 3109
#> 
#> Initializing model
#> 
#
# Print results for all imputed values
print(model.selection)
#> Call:
#> selection(model.eff = e ~ trt, model.cost = c ~ trt + e, model.mc = mc ~ 
#>     age, data = MenSS, dist_e = "norm", dist_c = "norm", type = "MAR", 
#>     n.chains = 2, n.iter = 500)
#> 
#> Coefficients:
#> $model.eff
#> (Intercept)        trt2 
#>       0.902       0.000 
#> 
#> $model.cost
#> (Intercept)        trt2           e 
#>      85.113      26.296      -7.200 
#> 
#> $model.me
#> (Intercept) 
#>       0.901 
#> 
#> $model.mc
#> (Intercept)         age 
#>       1.163      -0.009 
#> 

# Use looic to assess model fit
pic.looic <- pic(model.selection, criterion = "looic", cases = "cc")
pic.looic
#> $elpd
#> [1] -403.6339
#> 
#> $elpd_se
#> [1] -403.6339
#> 
#> $p
#> [1] 8.956858
#> 
#> $p_se
#> [1] 8.956858
#> 
#> $looic
#> [1] 807.2679
#> 
#> $looic_se
#> [1] 807.2679
#> 
#> $pointwise
#>         elpd_loo mcse_elpd_loo      p_loo    looic influence_pareto_k
#>  [1,]  -7.825187    0.01084903 0.05605215 15.65037        0.086099763
#>  [2,]  -7.858576    0.01077115 0.05490270 15.71715        0.110453554
#>  [3,] -12.493000    0.04912260 0.57195573 24.98600        0.316528483
#>  [4,]  -8.788244    0.01332154 0.07742027 17.57649        0.141788319
#>  [5,]  -7.725899    0.01183095 0.06284891 15.45180        0.172268281
#>  [6,]  -7.834046    0.01331965 0.07601311 15.66809        0.191233680
#>  [7,] -12.398947    0.05526224 0.66473251 24.79789        0.321208512
#>  [8,] -12.379724    0.06630884 0.71888272 24.75945        0.588482370
#>  [9,]  -8.163500    0.01286260 0.07551407 16.32700        0.020272569
#> [10,]  -7.838328    0.01125006 0.06059511 15.67666        0.062940439
#> [11,]  -8.169266    0.01197165 0.06618442 16.33853        0.055122534
#> [12,]  -8.188012    0.01293412 0.07546762 16.37602       -0.006994716
#> [13,]  -7.744277    0.01124773 0.05782128 15.48855        0.149155230
#> [14,]  -8.143820    0.01264923 0.07323945 16.28764        0.032584511
#> [15,] -17.735473    0.21601046 2.73470239 35.47095        0.996619249
#> [16,]  -8.193593    0.01217213 0.06831791 16.38719        0.109353539
#> [17,]  -8.069692    0.01889814 0.13142893 16.13938        0.164236568
#> [18,]  -7.945904    0.01140006 0.06014884 15.89181        0.043717578
#> [19,]  -7.781479    0.01172601 0.06473184 15.56296        0.079289231
#> [20,]  -8.941887    0.01497480 0.09158008 17.88377        0.122967140
#> [21,]  -8.908918    0.02069401 0.16492814 17.81784        0.177014504
#> [22,]  -7.844166    0.01094485 0.05710022 15.68833        0.096293516
#> [23,]  -9.342772    0.01824613 0.14627524 18.68554        0.232388442
#> [24,]  -8.197209    0.01170397 0.06289414 16.39442        0.184629074
#> [25,]  -8.162817    0.01197596 0.06619461 16.32563        0.078915925
#> [26,]  -8.115835    0.01260951 0.07248086 16.23167       -0.030231155
#> [27,]  -8.005677    0.01566634 0.10300759 16.01135        0.130348810
#> [28,]  -7.915311    0.01279257 0.07293672 15.83062        0.325508159
#> [29,]  -8.209591    0.01366196 0.08342116 16.41918        0.145373751
#> [30,]  -8.131228    0.01426611 0.08871708 16.26246        0.125425976
#> [31,]  -9.223388    0.01824655 0.13220601 18.44678        0.262940342
#> [32,]  -8.158981    0.01517186 0.09843337 16.31796        0.129923124
#> [33,]  -9.607338    0.02975666 0.30526974 19.21468        0.100179661
#> [34,]  -7.827069    0.01183252 0.06368886 15.65414        0.072849432
#> [35,]  -8.196958    0.02782186 0.22327056 16.39392        0.224674986
#> [36,]  -7.968912    0.01226385 0.06831006 15.93782        0.052673669
#> [37,]  -9.858171    0.03335463 0.37111768 19.71634        0.174013817
#> [38,]  -8.186854    0.01396085 0.08647280 16.37371        0.152042710
#> [39,]  -8.159541    0.01449480 0.09105226 16.31908        0.183271473
#> [40,]  -8.285482    0.01382556 0.08107018 16.57096        0.209816564
#> [41,]  -8.249914    0.01398655 0.08655397 16.49983        0.132899162
#> [42,]  -8.283421    0.01788697 0.11566132 16.56684        0.264833315
#> [43,]  -8.793350    0.01582738 0.10791532 17.58670        0.239223387
#> [44,]  -8.262877    0.01638371 0.10895738 16.52575        0.402241894
#> [45,]  -8.031047    0.01332111 0.07688793 16.06209        0.220000965
#> [46,]  -9.488250    0.02140834 0.17949463 18.97650        0.151977252
#> 
#> $diagnostics
#> $diagnostics$pareto_k
#>  [1]  0.086099763  0.110453554  0.316528483  0.141788319  0.172268281
#>  [6]  0.191233680  0.321208512  0.588482370  0.020272569  0.062940439
#> [11]  0.055122534 -0.006994716  0.149155230  0.032584511  0.996619249
#> [16]  0.109353539  0.164236568  0.043717578  0.079289231  0.122967140
#> [21]  0.177014504  0.096293516  0.232388442  0.184629074  0.078915925
#> [26] -0.030231155  0.130348810  0.325508159  0.145373751  0.125425976
#> [31]  0.262940342  0.129923124  0.100179661  0.072849432  0.224674986
#> [36]  0.052673669  0.174013817  0.152042710  0.183271473  0.209816564
#> [41]  0.132899162  0.264833315  0.239223387  0.402241894  0.220000965
#> [46]  0.151977252
#> 
#> $diagnostics$n_eff
#>  [1] 471.67347 472.02329 226.00082 458.31695 466.62248 458.62870 197.57093
#>  [8] 155.78956 461.32518 469.86457 466.11839 460.69744 469.58971 462.45032
#> [15]  20.18279 465.32518 422.98997 468.85731 467.38013 448.08461 410.53153
#> [22] 471.29278 428.85384 467.77501 466.11219 462.46079 444.73115 462.50807
#> [29] 456.88386 452.82471 428.26385 447.60327 345.14350 466.49948 358.87847
#> [36] 464.33380 319.95940 455.10392 451.64699 456.18543 454.50240 429.82377
#> [43] 443.47072 440.78445 459.00952 406.19272
#> 
#> $diagnostics$r_eff
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [39] 1 1 1 1 1 1 1 1
#> 
#> 
#> $psis_object
#> NULL
#> 

# Show density plots for mean costs parameters
diag.den <- diagnostic(model.selection, type = "denplot", param = "mu.c")

# Plots of imputations for all data
p1 <- plot(model.selection, class = "scatter", outcome = "both")

# Summarise the CEA results
summary(model.selection)
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under MAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.902 0.023 0.860 0.945
#> 2 0.902 0.027 0.849 0.956
#> 
#>  
#>  Mean costs by intervention 
#>   Mean   SD   QL  QU
#> 1   85 29.1 28.1 142
#> 2  112 36.3 36.0 185
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 45017 1132 42910 47192
#> 2 44971 1346 42349 47716

# }
#
#
```
