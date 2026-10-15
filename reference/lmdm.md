# Bayesian Longitudinal Missing Data Models

Full Bayesian cost-effectiveness models to handle missing data in
longitudinal outcomes under different missing data mechanism
assumptions, using alternative parametric distributions for the effect
and cost variables. The analysis is performed using the `BUGS` language,
which is implemented in the software `JAGS` using the function
[`jags`](https://rdrr.io/pkg/R2jags/man/jags.html) The output is stored
in an object of class 'missingHE'.

## Usage

``` r
lmdm(
  data,
  model.eff,
  model.cost,
  model.me = me ~ 1,
  model.mc = mc ~ 1,
  dist_e,
  dist_c,
  type,
  time_dep = "AR1",
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

  A data frame in which to find the longitudinal variables supplied in
  `model.eff`, `model.cost` (model formulas for effects and costs) and
  `model.me`, `model.mc` (model formulas for the missing effect and cost
  models). Among these, effectiveness, cost, time and treatment
  indicator (only two arms) variables must always be provided and named
  'e', 'c', 'time' and 'trt', respectively.

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
  response must be indicated with the term 'me' (missing effects) and
  any covariates must be provided on the right-hand side of the formula.
  If there are no covariates, `1` should be specified on the right hand
  side of the formula. By default, covariates are placed on the
  "probability" parameter for the missing effects through a
  logistic-linear model. Random effects can also be specified for each
  model parameter. See details for how these can be specified.

- model.mc:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be indicated with the term 'mc' (missing costs) and any
  covariates must be provided on the right-hand side of the formula. If
  there are no covariates, `1` should be specified on the right hand
  side of the formula. By default, covariates are placed on the
  "probability" parameter for the missing costs through a
  logistic-linear model. Random effects can also be specified for each
  model parameter. See details for how these can be specified.

- dist_e:

  Distribution assumed for the effects. Current available chocies are:
  Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
  Weibull ('weib'), Logistic ('logis'), Poisson ('pois'), Negative
  Binomial ('negbin') or Bernoulli ('bern').

- dist_c:

  Distribution assumed for the costs. Current available chocies are:
  Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').

- type:

  Type of missingness mechanism assumed. Choices are Missing At Random
  (MAR) and Missing Not At Random (MNAR).

- time_dep:

  Type of dependence structure assumed between effectiveness and cost
  outcomes. Current choices include: autoregressive structure of order
  one ('AR1') - default, bivariate at each time ('biv') and independence
  ('none').

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

- time_dep:

  A character variable that indicate which type of time dependence
  assumption was made, either `none` or `AR1`

## Details

Depending on the distributions specified for the outcome variables in
the arguments `dist_e` and `dist_c` and the type of missingness
mechanism specified in the argument `type`, different models are built
and run in the background by the function `lmdm`. These models fit
multinomial logistic regressions to estimate the probability of a given
missing data pattern `k` (1 = completers, 2 = intermittent, 3 = dropout)
in one or both longitudinal outcomes. A simple example can be used to
show how longitudinal missing data models are specified. Consider a
longitudinal data set comprising a response variable \\y\\ measured at
`S` occasions and a set of centered covariate \\X_j\\, for \\i = j, ...,
J\\. For each subject in the trial \\i = 1, ..., n\\ and time \\s = 1,
..., S\\ we define an indicator variable \\m_i\\ taking value `k = 1` if
the \\i\\-th individual is associated with no missing value (completer),
a value `k = 2` for intermittent missingness over the study period, and
a value `k = 3` for dropout. This is modelled as: \$\$m_i ~
Multinomial(\pi^k_i)\$\$ \$\$\pi^k_i = \phi^k_i/\sum\phi_i\$\$
\$\$log(\phi^k_i) = \sum\gamma^k_j X_j + \delta^k y_i\$\$ where

- \\\pi^k_i\\ is the individual probability of a missing value in \\y\\
  for pattern \\k\\ at a given time.

- \\\gamma^k_j\\ represents the impact on the missing data probability
  in \\y\\ of the covariate \\X_j\\ for pattern \\k\\ at a given time.

- \\\delta^k\\ represents the impact on the missing data probability in
  \\y\\ for pattern \\k\\ of missingness itself at a given time.

When \\\delta = 0\\ the model assumes a 'MAR' mechanism, while when
\\\delta != 0\\ the mechanism is 'MNAR'. For the parameters indexing the
missingness model, the default prior distributions assumed are:

- \\\gamma^k_j ~ Normal(0, 0.01)\\

- \\\delta^k ~ Normal(0, 1)\\

When user-defined hyperprior values are supplied via the argument
`prior` in the function `lmdm`, the elements of this list (see
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
function `lmdm` (see Arguments).

For each model, random effects can also be specified for each parameter
by adding the term + (x \| z) to each model formula, where x is the
fixed regression coefficient for which also the random effects are
desired and z is the clustering variable across which the random effects
are specified (must be the name of a factor variable in the dataset).

## References

Mason, AJ. Gomes, M. Carpenter, J. Grieve, R. (2021). *Flexible Bayesian
longitudinal models for cost‐effectiveness analyses with informative
missing data*. Health economics, 30(12), 3138-3158.

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
# Quick example to run using subset of PBS dataset

# Load longitudinal dataset

PBS.long <- PBS

# \donttest{
# Run the model using the long_miss function assuming a MAR mechanism
# Use only 100 iterations to run a quick check
model.long <- lmdm(data = PBS.long, model.eff = e ~ trt, model.cost = c ~ trt,
   model.me = me ~ 1, model.mc = mc ~ 1, dist_e = "norm", dist_c = "norm",
   type = "MAR", n.chains = 2, n.iter = 100, time_dep = "none")
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 2852
#>    Unobserved stochastic nodes: 112
#>    Total graph size: 16272
#> 
#> Initializing model
#> 

# Extract regression coefficient estimates
coef(model.long)
#> $Effects
#>                Mean    SD    QL    QU
#> (Intercept).1 0.475 0.035 0.408 0.533
#> trt2.1        0.090 0.046 0.012 0.190
#> (Intercept).2 0.501 0.031 0.443 0.556
#> trt2.2        0.138 0.047 0.039 0.230
#> (Intercept).3 0.478 0.029 0.419 0.526
#> trt2.3        0.138 0.040 0.064 0.207
#> 
#> $Costs
#>                  Mean     SD      QL      QU
#> (Intercept).1 147.799 48.826  62.582 248.354
#> trt2.1        100.570 39.491  34.235 167.534
#> (Intercept).2  73.607 33.988  17.269 139.259
#> trt2.2         45.972 28.579 -23.516  94.459
#> (Intercept).3  49.996 30.584 -21.641 104.961
#> trt2.3         30.893 31.002 -21.537  95.492
#> 
#

# Summarise the CEA information from the model
summary(model.long)
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under MAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.491 0.032 0.439 0.554
#> 2 0.615 0.031 0.554 0.673
#> 
#>  
#>  Mean costs by intervention 
#>   Mean  SD   QL   QU
#> 1  244 435 -555 1108
#> 2  451 450 -467 1340
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 24292 1669 21540 27643
#> 2 30293 1551 27081 33234

# Further examples which take longer to run
model.long <- lmdm(data = PBS.long, model.eff = e ~ trt, model.cost = c ~ trt + age,
   model.me = me ~ 1, model.mc = mc ~ 1, dist_e = "norm", dist_c = "norm",
   type = "MAR", n.chains = 2, n.iter = 500, time_dep = "none")
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 2852
#>    Unobserved stochastic nodes: 115
#>    Total graph size: 16835
#> 
#> Initializing model
#> 

# Use looic to assess model fit
pic.looic <- pic(model.long, criterion = "looic")
pic.looic
#> $elpd
#> [1] -6967.144
#> 
#> $elpd_se
#> [1] -6967.144
#> 
#> $p
#> [1] 59.2078
#> 
#> $p_se
#> [1] 59.2078
#> 
#> $looic
#> [1] 13934.29
#> 
#> $looic_se
#> [1] 13934.29
#> 
#> $pointwise
#>         elpd_loo mcse_elpd_loo        p_loo     looic influence_pareto_k
#> V1    -15.272257   0.023244564  0.209262495  30.54451       0.2805333777
#> V2     -9.311922   0.004674124  0.010396941  18.62384       0.2798825319
#> V3    -27.581939   0.126326867  1.730580967  55.16388       0.6842019030
#> V4    -14.858580   0.022239595  0.187234442  29.71716       0.3949851814
#> V5     -9.723690   0.005611725  0.015182929  19.44738       0.0804319967
#> V6    -10.117244   0.006817687  0.021809680  20.23449       0.1044244198
#> V7     -9.065612   0.003477368  0.006030832  18.13122      -0.0784228410
#> V8     -9.387119   0.004264927  0.008905620  18.77424       0.1169005633
#> V9     -9.493108   0.005009949  0.011941109  18.98622       0.1701014884
#> V10    -9.315224   0.004413534  0.009350672  18.63045       0.2632665408
#> V11    -9.363290   0.004676929  0.010686116  18.72658       0.0613053786
#> V12   -11.829521   0.011711540  0.061764458  23.65904       0.1281075094
#> V13    -9.008656   0.003636249  0.006525968  18.01731      -0.0451562902
#> V14    -9.900820   0.006668521  0.020801061  19.80164       0.1642131849
#> V15    -9.535986   0.004429394  0.009512768  19.07197       0.2337729773
#> V16    -8.961408   0.002925715  0.004301030  17.92282      -0.2180201782
#> V17   -15.074210   0.022444863  0.188790681  30.14842       0.3986629207
#> V18    -9.276591   0.004174202  0.008394571  18.55318       0.2156351498
#> V19   -10.249306   0.008652173  0.034841647  20.49861       0.1002983996
#> V20    -9.826276   0.007805444  0.028597279  19.65255      -0.1123770110
#> V21   -11.335207   0.009487803  0.041897428  22.67041       0.1484936882
#> V22    -9.307140   0.004007281  0.007752440  18.61428       0.1269803549
#> V23    -9.708290   0.007689968  0.027738650  19.41658      -0.0781622532
#> V24   -16.293034   0.029449111  0.280527374  32.58607       0.3563426528
#> V25    -9.184120   0.004074638  0.008064867  18.36824       0.2587813273
#> V26   -13.118218   0.020510179  0.167765773  26.23644       0.1668654374
#> V27    -9.513061   0.005392453  0.014185892  19.02612      -0.1340372846
#> V28    -9.543335   0.005361634  0.014053935  19.08667      -0.1051602576
#> V29    -9.518947   0.005444538  0.014402501  19.03789      -0.1323835955
#> V30    -9.808219   0.006423218  0.020181262  19.61644      -0.1285920228
#> V31    -9.596235   0.005312524  0.013798159  19.19247      -0.0941039065
#> V32    -9.775851   0.005635665  0.015276273  19.55170       0.0837300021
#> V33   -18.436826   0.030362605  0.306122171  36.87365       0.3627345285
#> V34    -9.150795   0.004014340  0.007885190  18.30159       0.0267324679
#> V35   -11.402749   0.008370438  0.033213200  22.80550       0.0622086956
#> V36    -9.194808   0.003855420  0.007373717  18.38962       0.0236363602
#> V37    -9.706027   0.004751721  0.011099563  19.41205       0.0261040656
#> V38    -9.556415   0.005366299  0.014096975  19.11283      -0.0549847466
#> V39    -9.560783   0.005343424  0.013957367  19.12157      -0.0962013255
#> V40    -8.921417   0.003162818  0.005009648  17.84283       0.0752005956
#> V41   -14.793419   0.021747545  0.185764457  29.58684       0.2924563746
#> V42    -8.857503   0.003080786  0.004760958  17.71501       0.0247981090
#> V43    -9.571344   0.005354730  0.014039069  19.14269      -0.0478509166
#> V44    -9.529011   0.005527571  0.014869566  19.05802      -0.0636718474
#> V45    -9.062423   0.003364809  0.005612545  18.12485      -0.0054117413
#> V46   -10.503639   0.008574549  0.034103685  21.00728       0.1604101631
#> V47    -9.520271   0.004633848  0.010489388  19.04054       0.0832472819
#> V48    -9.431876   0.004713757  0.010813014  18.86375       0.1582002285
#> V49   -11.526120   0.012450987  0.066402475  23.05224      -0.0311433121
#> V50   -11.753862   0.012738062  0.068994830  23.50772       0.0304137309
#> V51    -9.189620   0.003707378  0.006785240  18.37924       0.0522401237
#> V52    -8.925079   0.003115968  0.004838854  17.85016      -0.0227890113
#> V53    -9.137219   0.003782454  0.006964182  18.27444       0.1450307609
#> V54    -9.121739   0.003625545  0.006534662  18.24348       0.0484424160
#> V55    -9.050950   0.003235722  0.005233233  18.10190       0.0661998807
#> V56   -16.900087   0.030936814  0.321960659  33.80017       0.4185708151
#> V57    -9.071222   0.003962597  0.007685701  18.14244      -0.0854922568
#> V58    -8.976448   0.003170515  0.004956062  17.95290      -0.0241711073
#> V59    -9.070745   0.003767972  0.007023008  18.14149      -0.0246322754
#> V60    -9.176540   0.004042922  0.007909025  18.35308       0.0574937647
#> V61    -9.097992   0.004364791  0.009399591  18.19598      -0.1463233047
#> V62    -8.833238   0.003058599  0.004633688  17.66648       0.1207596331
#> V64    -9.061448   0.003728292  0.006885777  18.12290       0.0046194918
#> V65    -9.506463   0.004359600  0.009360933  19.01293      -0.1743627590
#> V66    -9.698015   0.006276366  0.019177915  19.39603      -0.0948302136
#> V67    -9.290186   0.004125579  0.008342343  18.58037       0.1204100158
#> V68   -14.199800   0.018259527  0.141172435  28.39960       0.1692279982
#> V69    -9.219167   0.004026741  0.007845342  18.43833       0.1519561549
#> V70    -9.553466   0.005347666  0.013967864  19.10693      -0.1054436232
#> V71    -8.862678   0.003199621  0.005076748  17.72536      -0.1675144779
#> V72    -8.961223   0.003127749  0.004892707  17.92245       0.1110184046
#> V74    -8.850532   0.003008818  0.004520829  17.70106       0.1704365351
#> V75    -8.905586   0.003274759  0.005300997  17.81117      -0.1044293986
#> V76    -9.476752   0.003651280  0.006641654  18.95350      -0.3108710849
#> V77    -9.172622   0.003914353  0.007616547  18.34524      -0.0814861310
#> V78    -9.015141   0.003498073  0.006017204  18.03028      -0.1483343527
#> V79    -9.889103   0.005223270  0.013416338  19.77821      -0.0667730513
#> V80    -8.971891   0.003107567  0.004829046  17.94378       0.1140641717
#> V81    -9.065064   0.003398896  0.005755723  18.13013       0.0693816366
#> V82    -9.048091   0.003680226  0.006723049  18.09618       0.0368815466
#> V83    -9.688239   0.006074853  0.018069270  19.37648      -0.1169148885
#> V84   -10.347788   0.008810435  0.037556401  20.69558      -0.0597296770
#> V85    -8.954553   0.002980509  0.004440708  17.90911       0.0663357224
#> V86    -9.308953   0.003431040  0.005886858  18.61791      -0.1041650067
#> V87   -14.995715   0.021958314  0.186897382  29.99143       0.3212338584
#> V88    -8.901588   0.003290125  0.005369777  17.80318      -0.1943427635
#> V89   -11.354362   0.009358015  0.042373335  22.70872       0.0029059434
#> V90    -9.273440   0.003926240  0.007659458  18.54688       0.0558260550
#> V91    -8.879990   0.003110017  0.004794676  17.75998       0.0335786921
#> V92    -9.033600   0.003685918  0.006699186  18.06720      -0.1215922457
#> V93    -8.856074   0.002963645  0.004384282  17.71215       0.0981274980
#> V94    -9.311262   0.004752215  0.011019090  18.62252      -0.0746334183
#> V95    -9.051863   0.003261859  0.005250730  18.10373      -0.0246700002
#> V96   -11.045572   0.010593641  0.054551712  22.09114       0.0460348826
#> V98   -11.371840   0.010409383  0.050348284  22.74368       0.0670681359
#> V99    -9.098039   0.003829000  0.007229978  18.19608      -0.0207851553
#> V100   -8.971335   0.003029002  0.004605910  17.94267      -0.1782799903
#> V101   -9.440476   0.005504814  0.014752615  18.88095      -0.1255035154
#> V103  -11.653750   0.012428854  0.067042537  23.30750       0.0237519147
#> V104  -11.845666   0.009813345  0.045821994  23.69133       0.0679384913
#> V105   -9.073615   0.003580565  0.006283674  18.14723       0.1168499688
#> V106   -8.934826   0.003040809  0.004604809  17.86965       0.0381957834
#> V107   -9.800341   0.006519321  0.019811549  19.60068       0.0973203980
#> V109  -15.082091   0.021553113  0.184816709  30.16418       0.3868430964
#> V110   -9.113832   0.003624205  0.006527258  18.22766       0.0479198283
#> V111   -9.952601   0.006544945  0.020445801  19.90520       0.0786809141
#> V112   -9.426978   0.004708929  0.010887561  18.85396      -0.1828303202
#> V113   -9.059797   0.003731924  0.006899462  18.11959      -0.0003820443
#> V114   -9.320956   0.005360980  0.013905657  18.64191      -0.0820497994
#> V115   -9.391559   0.005142963  0.012906926  18.78312      -0.0826846517
#> V116  -10.798562   0.006467164  0.019722768  21.59712       0.1828673186
#> V117   -9.283513   0.004456846  0.009725722  18.56703       0.1347254807
#> V118   -9.735656   0.004232982  0.008732678  19.47131      -0.1211637032
#> V119   -9.806653   0.006281756  0.018869577  19.61331       0.2476033515
#> V120  -10.644443   0.006763860  0.022327045  21.28889      -0.0765529844
#> V121  -11.014105   0.007776702  0.028368102  22.02821       0.2281717968
#> V122   -9.268861   0.004202785  0.008506237  18.53772       0.2362509445
#> V123   -9.965627   0.006606358  0.020412686  19.93125       0.1766718400
#> V124  -10.114700   0.007668917  0.027840433  20.22940       0.0849041436
#> V125  -10.034730   0.006825873  0.021852588  20.06946       0.1256285979
#> V126  -12.959871   0.013742618  0.084343810  25.91974       0.0271105966
#> V127   -9.853173   0.006496737  0.019682587  19.70635       0.1362733616
#> V128  -10.366940   0.008686441  0.035314148  20.73388       0.0792951380
#> V129   -9.789016   0.006516864  0.019788500  19.57803       0.1076717298
#> V133   -9.025372   0.002952511  0.004346014  18.05074      -0.0486681161
#> V134   -8.827664   0.003017147  0.004538785  17.65533      -0.0163138966
#> V135   -9.295705   0.004314263  0.008935228  18.59141       0.2046934068
#> V136   -9.459922   0.004686503  0.010744523  18.91984       0.0778503259
#> V137   -8.909209   0.003050925  0.004636329  17.81842       0.0320218425
#> V138   -9.219616   0.004026482  0.007844447  18.43923       0.1517138798
#> V139   -9.840758   0.006315450  0.019108742  19.68152       0.1354908627
#> V140   -9.690138   0.005961042  0.017090151  19.38028       0.2005940539
#> V141   -9.227063   0.004028629  0.007848739  18.45413       0.1729223357
#> V142  -11.074241   0.009359091  0.039825181  22.14848       0.0999930365
#> V143   -9.411552   0.004381084  0.009249720  18.82310       0.1584053306
#> V144   -8.997655   0.003262744  0.005287726  17.99531       0.0088337301
#> V145   -8.992451   0.002931204  0.004322447  17.98490      -0.2941180034
#> V146   -9.470908   0.004507461  0.009792873  18.94182       0.1908629646
#> V147  -15.899917   0.027547627  0.255242102  31.79983       0.3447135309
#> V148   -9.850749   0.006488395  0.019630412  19.70150       0.1172721549
#> V149  -14.771497   0.021837007  0.187070831  29.54299       0.2851359123
#> V150   -9.148577   0.004027675  0.007983610  18.29715      -0.0471988161
#> V151   -8.932848   0.003083880  0.004761263  17.86570       0.0911048742
#> V152   -8.846846   0.003119477  0.004843616  17.69369      -0.1552753499
#> V153  -13.622754   0.023178140  0.221603454  27.24551       0.3008410199
#> V154   -9.642254   0.005866894  0.016480636  19.28451       0.0031703482
#> V155  -11.076951   0.010170585  0.050276555  22.15390       0.0190834523
#> V156  -15.686296   0.021869693  0.194180526  31.37259       0.2240273580
#> V157   -9.840144   0.006107165  0.017811375  19.68029       0.1176278088
#> V158  -10.395129   0.007867975  0.029082520  20.79026       0.0396945883
#> V159   -9.514209   0.005580849  0.015252352  19.02842      -0.1239588789
#> V160   -9.302204   0.004507853  0.009840714  18.60441       0.1647292924
#> V161   -9.271801   0.004332300  0.009012866  18.54360       0.1857553852
#> V163   -9.837832   0.004088667  0.008130891  19.67566       0.1463164476
#> V164   -8.929271   0.003224946  0.005176938  17.85854      -0.0138781001
#> V165  -18.918724   0.032722806  0.329890034  37.83745       0.4333877144
#> V166  -10.047847   0.006504259  0.020158402  20.09569       0.0709722505
#> V167   -9.229817   0.004362204  0.009143528  18.45963       0.1173450926
#> V168   -8.905861   0.003212804  0.005142382  17.81172      -0.0001371779
#> V170   -8.995232   0.003006870  0.004524523  17.99046      -0.1568834611
#> V171   -9.117932   0.003852475  0.007329114  18.23586       0.0971300663
#> V172   -9.049077   0.003314241  0.005448766  18.09815       0.0589456646
#> V173   -8.934780   0.003078626  0.004725958  17.86956       0.0222482568
#> V174  -15.349296   0.022341552  0.191560001  30.69859       0.3378672494
#> V175   -9.275335   0.004151161  0.008300391  18.55067       0.2113421999
#> V176   -9.465805   0.003115585  0.004824397  18.93161      -0.0955657571
#> V177  -14.722162   0.022116347  0.185862099  29.44432       0.3889960179
#> V178  -11.039986   0.010343145  0.051877765  22.07997       0.0328342663
#> V179   -9.688485   0.004646313  0.010600048  19.37697      -0.1544669882
#> V180   -9.799887   0.005538955  0.015159595  19.59977      -0.1155508701
#> V181  -10.129625   0.006649671  0.020965846  20.25925       0.1485962899
#> V182   -9.056148   0.003722782  0.006868205  18.11230       0.0170324614
#> V183   -9.967865   0.006779308  0.021936197  19.93573       0.0579950015
#> V184   -8.932961   0.003127031  0.004882105  17.86592      -0.0692372189
#> V185   -9.438623   0.004841436  0.011467002  18.87725      -0.1162129358
#> V186   -8.881738   0.003145158  0.004925270  17.76348       0.0319342873
#> V187  -10.155003   0.006023625  0.017362129  20.31001       0.1894097448
#> V188  -10.436513   0.007639613  0.027431092  20.87303       0.0538888247
#> V189   -9.128188   0.003158301  0.004989949  18.25638       0.0181889305
#> V190   -9.605206   0.005858314  0.016690474  19.21041      -0.0519280013
#> V191  -11.654114   0.014178475  0.093129604  23.30823       0.1237276424
#> V194   -8.840699   0.003076917  0.004735667  17.68140      -0.0486912546
#> V195  -37.656289   0.353453346  3.916454134  75.31258       1.1120057878
#> V196  -12.244117   0.012535589  0.067533062  24.48823       0.2411138600
#> V197  -12.966115   0.020568677  0.169379314  25.93223       0.1198929207
#> V198   -8.947563   0.003205599  0.005112208  17.89513       0.0450122774
#> V199   -9.677641   0.005073232  0.012438039  19.35528       0.0489632334
#> V200  -15.503689   0.026776236  0.253138201  31.00738       0.3326752715
#> V201  -10.515755   0.007135162  0.024415314  21.03151       0.0506331409
#> V202  -11.735423   0.010155995  0.046666945  23.47085       0.0082864543
#> V204   -9.410210   0.004649183  0.010614730  18.82042      -0.1443565516
#> V205   -9.186178   0.004018354  0.008022982  18.37236       0.0220678895
#> V206   -8.955990   0.003366629  0.005660149  17.91198       0.0858328269
#> V207   -9.125297   0.003615176  0.006498698  18.25059       0.0301340470
#> V208   -9.367859   0.004367574  0.009170196  18.73572       0.2231426589
#> V209   -9.850508   0.006551872  0.020033515  19.70102       0.1688319904
#> V210   -9.605932   0.005660375  0.015428339  19.21186      -0.0505029186
#> V211   -9.554925   0.005758896  0.016177124  19.10985      -0.1033426612
#> V212   -9.651822   0.005143443  0.012795607  19.30364       0.0621011592
#> V213   -9.908641   0.006558852  0.020100590  19.81728       0.1572224772
#> V214   -9.277308   0.004297209  0.009079077  18.55462       0.0474061865
#> V215   -9.123115   0.003859605  0.007319903  18.24623      -0.0304810997
#> V216   -9.519876   0.005380044  0.014127844  19.03975      -0.1101512672
#> V217   -9.212218   0.003820297  0.007171706  18.42444      -0.0491571479
#> V218   -9.126559   0.003625376  0.006534675  18.25312       0.0500886534
#> V219   -9.594501   0.005689907  0.015577731  19.18900      -0.0475414631
#> V220   -8.898659   0.003246115  0.005208314  17.79732      -0.1172145979
#> V221   -9.202275   0.004040850  0.008033654  18.40455      -0.0484870415
#> V222   -9.675113   0.004715222  0.010737850  19.35023       0.2678563208
#> V223   -9.755397   0.006461004  0.019700613  19.51079       0.1121199592
#> V224   -9.202695   0.004145854  0.008480315  18.40539      -0.0091878531
#> V225   -9.604168   0.005364591  0.013950859  19.20834       0.0634940898
#> V226   -9.530463   0.005450331  0.014421098  19.06093      -0.1194505658
#> V227   -9.766900   0.005670579  0.015510026  19.53380       0.0688889067
#> V228   -9.842311   0.006372755  0.019875336  19.68462      -0.0880703251
#> V229   -9.518205   0.005446592  0.014412131  19.03641      -0.1296671317
#> V230  -13.622137   0.023170855  0.221499684  27.24427       0.3128650512
#> V231   -9.378474   0.004425914  0.009436925  18.75695       0.1726484374
#> V233   -9.399830   0.004438351  0.009491634  18.79966       0.1754469284
#> V234   -9.865746   0.006501263  0.019708256  19.73149       0.1715828463
#> V235   -9.121738   0.003954447  0.007672191  18.24348       0.0063436550
#> V236  -15.689894   0.027211394  0.248659734  31.37979       0.4268420884
#> V237   -8.893041   0.003159552  0.004997432  17.78608       0.1006224819
#> V238   -9.429486   0.003930796  0.007697439  18.85897      -0.0579619256
#> V239  -15.910439   0.023197148  0.203767134  31.82088       0.3225848193
#> V240   -8.945520   0.003254248  0.005263369  17.89104       0.0532717550
#> V241   -9.587513   0.005608867  0.015169785  19.17503      -0.0592390638
#> V242   -9.374549   0.004379533  0.009228125  18.74910       0.2014025505
#> V243   -8.928343   0.003261399  0.005314077  17.85669       0.0147413536
#> V244   -9.503314   0.003694490  0.006807567  19.00663      -0.0802265372
#> V245   -8.855394   0.003706368  0.006828662  17.71079      -0.1570584625
#> V246   -8.969671   0.003934778  0.007859892  17.93934      -0.0578976911
#> V247   -8.980440   0.004284672  0.009064438  17.96088      -0.1113116507
#> V249  -10.644277   0.008950450  0.037600397  21.28855       0.2547562319
#> V250   -9.062674   0.003881689  0.007562423  18.12535      -0.0588953452
#> V251   -9.029299   0.003611036  0.006532144  18.05860      -0.1538674513
#> V252   -9.039241   0.003629458  0.006621197  18.07848      -0.1042570428
#> V253   -8.756696   0.003270776  0.005391511  17.51339       0.0056075354
#> V254   -9.015455   0.003624472  0.006602675  18.03091      -0.1248211949
#> V255   -9.462572   0.005541124  0.015028239  18.92514       0.0928217518
#> V256   -9.167588   0.004257068  0.009022473  18.33518      -0.0590291043
#> V257   -9.104111   0.004246501  0.009030413  18.20822       0.1190651704
#> V258   -8.978954   0.003677897  0.006787926  17.95791      -0.0314569009
#> V259  -17.341318   0.041878677  0.529365108  34.68264       0.1954137246
#> V260   -8.795191   0.003253786  0.005304532  17.59038      -0.1426470532
#> V262   -8.782568   0.003283138  0.005398803  17.56514      -0.1114672029
#> V263   -9.879279   0.007349590  0.026128452  19.75856      -0.1245202263
#> V264  -10.147859   0.008316182  0.033107359  20.29572      -0.1364529135
#> V266   -9.129312   0.004185969  0.008784236  18.25862       0.1582425417
#> V267   -9.121465   0.004759899  0.011188878  18.24293       0.0681564583
#> V269   -9.326895   0.004679287  0.010886457  18.65379       0.0542303840
#> V270  -10.027370   0.007761906  0.028462027  20.05474       0.2433232375
#> V271   -8.933650   0.003963726  0.007818594  17.86730       0.0709814303
#> V272   -8.855845   0.003776544  0.007135357  17.71169       0.0732142630
#> V273   -9.289194   0.005535023  0.014360709  18.57839       0.1806243116
#> V274  -10.033869   0.008055566  0.030182998  20.06774       0.0457680548
#> V275   -9.340664   0.005336340  0.013536945  18.68133       0.0950844229
#> V276   -9.237151   0.004158845  0.008623595  18.47430      -0.0780743242
#> V278   -8.849889   0.003797214  0.007163168  17.69978       0.1016593364
#> V279  -10.253038   0.005769118  0.016118708  20.50608       0.0821952532
#> V280   -9.353397   0.005330722  0.013545951  18.70679       0.0655206993
#> V281   -8.848807   0.003207840  0.005165022  17.69761       0.0296837065
#> V282   -9.264607   0.005458213  0.014094313  18.52921       0.1449446555
#> V283   -9.058299   0.003531696  0.006251469  18.11660      -0.0599375196
#> V284   -9.259417   0.005478729  0.014185485  18.51883       0.1088705500
#> V286   -8.891018   0.003817743  0.007210413  17.78204       0.1034561798
#> V287   -8.740114   0.003468025  0.006074826  17.48023      -0.0097287123
#> V288   -8.917213   0.003856302  0.007346460  17.83443       0.1015928829
#> V289   -8.836924   0.003211803  0.005167373  17.67385      -0.1659320872
#> V290  -11.218121   0.011675011  0.061948294  22.43624       0.1475228442
#> V291  -11.627749   0.011176217  0.056717463  23.25550       0.1564023251
#> V292  -13.521309   0.018456330  0.133763337  27.04262       0.1567427394
#> V295   -9.328874   0.003936830  0.007621339  18.65775      -0.0531946033
#> V296   -8.997991   0.003535290  0.006155453  17.99598       0.1007910082
#> V297   -8.797196   0.003229621  0.005224666  17.59439      -0.1156180892
#> V298   -9.471523   0.005201881  0.012938504  18.94305       0.0510568621
#> V299   -8.855261   0.003768284  0.007099463  17.71052       0.0825192700
#> V300  -15.809983   0.033219503  0.381006004  31.61997       0.1115806100
#> V301   -8.878272   0.003955254  0.007753022  17.75654       0.0074122014
#> V302   -8.941100   0.003891153  0.007475650  17.88220       0.1448711303
#> V303   -8.956458   0.004018624  0.007907147  17.91292       0.1876937019
#> V304   -8.848037   0.003378166  0.005809143  17.69607      -0.0834683640
#> V305   -8.938531   0.003992375  0.007907522  17.87706       0.1292260175
#> V306   -8.816009   0.003704548  0.006883967  17.63202       0.0569002304
#> V307  -11.215273   0.009420663  0.041541843  22.43055      -0.0076087623
#> V308   -8.849460   0.003784850  0.007122288  17.69892       0.0849168704
#> V309   -8.908706   0.003997543  0.007860844  17.81741       0.1788158278
#> V310   -8.788577   0.003511394  0.006184623  17.57715      -0.0054972864
#> V311   -9.207264   0.004519256  0.010217914  18.41453       0.1739759225
#> V312  -12.131060   0.012433875  0.072923104  24.26212      -0.1120367077
#> V313   -9.130739   0.004376758  0.009619757  18.26148       0.1313413900
#> V314   -8.844801   0.003643462  0.006667637  17.68960       0.0331667401
#> V315   -8.792130   0.003629577  0.006624648  17.58426       0.0315516828
#> V316   -9.387529   0.005283267  0.013309774  18.77506       0.0702144616
#> V317  -11.678601   0.010180739  0.049469858  23.35720       0.1132014474
#> V318   -8.906458   0.003976035  0.007851922  17.81292       0.0715064739
#> V319   -8.879222   0.003983866  0.007859581  17.75844       0.0681046606
#> V320   -9.025780   0.003863661  0.007450232  18.05156      -0.0346584832
#> V321   -8.931471   0.003936221  0.007673331  17.86294      -0.1167949741
#> V322   -9.359982   0.004140558  0.008434064  18.71996       0.0236632585
#> V323   -9.230834   0.004450602  0.009678081  18.46167       0.1551879425
#> V324   -9.611876   0.005114103  0.012568336  19.22375      -0.0014398791
#> V325   -9.390878   0.005269972  0.013232887  18.78176       0.0525650847
#> V326   -9.268540   0.005508746  0.014291155  18.53708       0.1015315934
#> V327  -14.327614   0.025237857  0.233309528  28.65523       0.0998281643
#> V328   -9.873840   0.007525954  0.026584945  19.74768       0.0921396010
#> V329   -8.770001   0.003406767  0.005822760  17.54000       0.0375345713
#> V330   -8.859049   0.003771561  0.007117878  17.71810       0.0707023472
#> V332   -9.567762   0.006184307  0.018293932  19.13552      -0.0266916825
#> V333   -9.533519   0.005857277  0.016459510  19.06704       0.1273820643
#> V334   -9.543304   0.006144235  0.018053367  19.08661       0.1358734332
#> V335   -8.845841   0.003786560  0.007134646  17.69168       0.0775642789
#> V336   -8.876245   0.003795183  0.007132702  17.75249       0.1142094078
#> V337   -9.294145   0.005408325  0.013883808  18.58829       0.1021191237
#> V338   -9.221499   0.004861036  0.011551190  18.44300      -0.1999852214
#> V339   -9.794233   0.006966653  0.022950599  19.58847       0.0699775851
#> V340  -19.506225   0.052552601  0.904850875  39.01245       0.0206669551
#> V341  -11.409426   0.009640263  0.044580829  22.81885       0.0166143137
#> V342  -42.453392   0.246379000  6.685933962  84.90678       0.9776573253
#> V343   -9.321214   0.005425915  0.014050288  18.64243       0.1200897408
#> V344   -9.143773   0.004652459  0.010895879  18.28755       0.0657972685
#> V345   -9.936317   0.007434315  0.027159427  19.87263       0.1147409746
#> V346  -13.065552   0.011431375  0.060862071  26.13110      -0.0383839630
#> V349   -9.247781   0.005032421  0.012363837  18.49556       0.0435603486
#> V350   -8.984545   0.004245367  0.008896047  17.96909      -0.0762910118
#> V351  -10.017689   0.007072015  0.024236228  20.03538       0.2345514864
#> V352  -10.997149   0.009001660  0.038590230  21.99430       0.0026725119
#> V353  -14.747944   0.019412907  0.154618406  29.49589       0.1192217066
#> V354   -9.554869   0.006077327  0.017666258  19.10974       0.1251352428
#> V355   -8.811966   0.003689308  0.006825372  17.62393      -0.1303622447
#> V356   -8.873841   0.003468630  0.006023733  17.74768      -0.0660220284
#> V357   -9.245164   0.005130466  0.012798255  18.49033       0.1161042073
#> V358   -9.348067   0.005466777  0.014827396  18.69613       0.2051076517
#> V359   -8.758564   0.003519379  0.006244327  17.51713       0.0744559521
#> V360   -9.240369   0.004957341  0.012318005  18.48074       0.0738769662
#> V361   -9.720725   0.006784499  0.021971390  19.44145       0.2689155232
#> V362   -9.834872   0.007073368  0.024525043  19.66974       0.2638841057
#> V363  -10.210423   0.005807797  0.016235350  20.42085       0.0872861699
#> V364  -11.972349   0.012638978  0.074986317  23.94470      -0.0494513534
#> V365   -8.780563   0.003625168  0.006543998  17.56113      -0.0012265040
#> V366   -9.896564   0.007054692  0.024395764  19.79313       0.1831465400
#> V367   -9.828130   0.007099065  0.024762542  19.65626       0.1811045010
#> V368   -8.712514   0.003336632  0.005610806  17.42503       0.0009528022
#> V369   -9.915537   0.007178240  0.025278019  19.83107       0.2423494355
#> V370  -12.359557   0.012076642  0.069959762  24.71911       0.2352268441
#> V371   -8.795811   0.003135687  0.004943017  17.59162      -0.0660421640
#> V372  -11.750283   0.015843933  0.107276643  23.50057       0.2001702460
#> V373  -12.033191   0.013262542  0.078721180  24.06638       0.1968491073
#> V374  -11.990749   0.010861923  0.054707193  23.98150       0.0736398718
#> V375  -12.366233   0.011166548  0.057739755  24.73247       0.2503809925
#> V376  -18.394871   0.030016954  0.335390920  36.78974       0.2744156824
#> V377  -11.017542   0.009089236  0.038695608  22.03508      -0.2145835357
#> V378   -9.141049   0.003529929  0.006168261  18.28210      -0.1920492739
#> V379   -8.860828   0.003420614  0.005831861  17.72166      -0.1513661079
#> V380  -11.750280   0.012924992  0.074632865  23.50056       0.1099437721
#> V381   -9.173107   0.004567911  0.010462270  18.34621       0.1524923120
#> V382   -9.187363   0.004546512  0.010360379  18.37473       0.1313017341
#> V383   -8.917221   0.003144943  0.004955873  17.83444      -0.1576484562
#> V384   -8.999166   0.003758964  0.007115202  17.99833      -0.1367438579
#> V385   -9.207216   0.004521225  0.010223437  18.41443       0.1931438278
#> V387   -8.942230   0.003364707  0.005663643  17.88446      -0.2993539889
#> V388   -9.843757   0.007049752  0.024356837  19.68751       0.1836926476
#> V389   -8.869241   0.003121121  0.004873804  17.73848      -0.1135863116
#> V390   -9.825300   0.007092643  0.024709653  19.65060       0.1873193468
#> V392   -8.783809   0.003408374  0.005802578  17.56762      -0.0796331437
#> V394   -9.577437   0.005971905  0.017174252  19.15487       0.1039146172
#> V395   -9.264702   0.005458054  0.014095428  18.52940       0.1367063512
#> V396   -9.284590   0.005525711  0.014328374  18.56918       0.1578869870
#> V397  -12.868491   0.021257473  0.174801414  25.73698       0.1054930672
#> V398   -9.401233   0.005795593  0.015710658  18.80247       0.1717371704
#> V399   -8.814893   0.003702485  0.006874938  17.62979       0.0552662039
#> V401  -10.298610   0.008758441  0.035846122  20.59722       0.3197337529
#> V402  -10.000328   0.006856603  0.022798734  20.00066      -0.1006443333
#> V403   -9.880140   0.007574512  0.026939456  19.76028       0.0517736270
#> V404   -9.787029   0.006023661  0.017749595  19.57406      -0.0393945869
#> V405   -9.830243   0.007061064  0.024448098  19.66049       0.1823250102
#> V407  -11.620487   0.015308001  0.100587174  23.24097       0.1868590544
#> V408   -9.053571   0.004143394  0.008466569  18.10714      -0.0213122271
#> V410   -9.935177   0.007880864  0.029245650  19.87035       0.2336018814
#> V411   -9.774010   0.007124165  0.024923050  19.54802       0.1144355340
#> V412   -9.318329   0.004462901  0.009991469  18.63666       0.1548690381
#> V414   -8.828899   0.003175551  0.005048707  17.65780      -0.1291773284
#> V415   -8.983560   0.003927854  0.007779763  17.96712      -0.0987003231
#> V416   -9.084320   0.004239298  0.009050920  18.16864       0.0568135255
#> V417   -8.833340   0.003176354  0.005090814  17.66668      -0.0382827739
#> V419   -8.968795   0.003537159  0.006312896  17.93759      -0.1024059935
#> V420   -9.612110   0.003793594  0.007099788  19.22422      -0.0403553404
#> V421  -14.905464   0.018959561  0.149631850  29.81093       0.1423544454
#> V422   -9.458672   0.004984839  0.012184491  18.91734      -0.0587598412
#> V423   -8.876587   0.003962424  0.007782266  17.75317       0.0097263783
#> V424   -9.076353   0.003788458  0.007177845  18.15271      -0.1809447299
#> V425   -9.248733   0.003934107  0.007807407  18.49747      -0.0362733785
#> V426   -8.845726   0.003796440  0.007171912  17.69145       0.0716808053
#> V427  -10.338783   0.009322526  0.040134943  20.67757       0.2839127676
#> V428   -8.763330   0.003349379  0.005619172  17.52666      -0.1137876936
#> V429   -8.862628   0.003768124  0.007056980  17.72526      -0.0316493838
#> V430  -10.027190   0.007250516  0.025048818  20.05438       0.0833028923
#> V431  -10.582776   0.008838091  0.036622011  21.16555       0.2171338252
#> V432  -10.332388   0.007835703  0.029826238  20.66478       0.0670720132
#> V433   -9.521594   0.005146317  0.012679484  19.04319      -0.0013233432
#> V434   -8.843086   0.003794381  0.007165231  17.68617       0.0983013335
#> V435   -8.960956   0.003856083  0.007462193  17.92191      -0.0674521636
#> V436  -11.432174   0.009805828  0.045589593  22.86435      -0.0209751451
#> V437  -10.975827   0.008922112  0.038098687  21.95165       0.0508739569
#> V438   -8.857130   0.003626221  0.006522890  17.71426       0.0479799661
#> V439   -8.858053   0.003618162  0.006525966  17.71611       0.0288221153
#> V440  -11.473210   0.011380307  0.061859168  22.94642       0.0006679855
#> V441   -8.923391   0.003850663  0.007362282  17.84678      -0.1946406420
#> V442   -8.771932   0.003239514  0.005263323  17.54386      -0.0919918690
#> V443   -9.130839   0.004645732  0.010859470  18.26168       0.1162073396
#> V445  -13.014394   0.019484159  0.160395816  26.02879       0.0452538418
#> V447  -11.728911   0.009512390  0.042986433  23.45782       0.1349479450
#> V448   -9.484086   0.005006442  0.012414317  18.96817       0.0342864105
#> V449   -9.825933   0.007434551  0.025958113  19.65187      -0.0088882709
#> V450   -9.039608   0.004420710  0.009612875  18.07922      -0.1112715602
#> V451   -8.930723   0.003900864  0.007514898  17.86145       0.1443411647
#> V452   -9.993710   0.007045047  0.024266852  19.98742       0.1918237272
#> V453  -10.029447   0.007182791  0.025249720  20.05889       0.1787660055
#> V454   -9.336826   0.005603315  0.014664239  18.67365       0.2358305844
#> V455   -8.894867   0.003993853  0.007901812  17.78973       0.0707271333
#> V456   -9.095305   0.004660940  0.010649572  18.19061       0.0428504316
#> V457   -9.213380   0.004045758  0.008237017  18.42676       0.0241951667
#> V458   -9.741152   0.005676693  0.015822167  19.48230      -0.0130871949
#> V459  -10.028798   0.007777526  0.028313491  20.05760       0.0648240863
#> V460   -8.877981   0.003977483  0.007835196  17.75596       0.0609796460
#> V461   -9.821278   0.004741188  0.011146309  19.64256      -0.0478392554
#> V462   -9.263612   0.003977028  0.007815079  18.52722      -0.0673069177
#> V463   -8.791136   0.003448721  0.005937695  17.58227      -0.0067395588
#> V464   -8.782255   0.003484267  0.006130970  17.56451       0.0395423664
#> V465  -15.384520   0.025032200  0.259068376  30.76904      -0.0005102754
#> V466   -9.306704   0.005401345  0.014174879  18.61341       0.0803043654
#> V467   -9.439455   0.006042866  0.017059234  18.87891       0.2443770451
#> V468   -9.284879   0.004464800  0.010022102  18.56976       0.0305322852
#> V469   -9.397603   0.004744899  0.011141206  18.79521      -0.0652615347
#> V470   -9.526741   0.005388403  0.014016282  19.05348       0.0870622976
#> V471   -9.809564   0.006599719  0.020970428  19.61913       0.3140729831
#> V472   -9.029985   0.003388548  0.005765621  18.05997      -0.0416206634
#> V473   -8.865521   0.003799076  0.007155036  17.73104       0.1056485624
#> V474   -8.903964   0.003966274  0.007769019  17.80793       0.1768318406
#> V475   -8.951862   0.003179737  0.005069387  17.90372      -0.1035285292
#> V476  -11.315504   0.009241980  0.040578931  22.63101       0.0847414961
#> V477   -9.357917   0.004528757  0.010203231  18.71583       0.1168810137
#> V478   -9.898103   0.007041975  0.024295204  19.79621       0.1884031312
#> V479   -8.916849   0.004026144  0.008064794  17.83370       0.0363794076
#> V481   -8.878964   0.003840896  0.007408638  17.75793      -0.0389227394
#> V482   -9.026557   0.003615868  0.006560318  18.05311      -0.0728630292
#> V484   -9.129953   0.004654284  0.010901563  18.25991       0.1101451980
#> V485   -8.912579   0.003825946  0.007237474  17.82516       0.0875427947
#> V486   -9.207531   0.004152106  0.008617290  18.41506       0.0977772123
#> V487   -8.750112   0.003339920  0.005590106  17.50022       0.0047488296
#> V488   -9.780010   0.004582213  0.010407869  19.56002       0.0017236815
#> V489   -9.075786   0.003095573  0.004740892  18.15157       0.0255844892
#> V490   -9.377086   0.004044757  0.007903545  18.75417       0.1917546804
#> V491   -9.590043   0.004701303  0.011090690  19.18009       0.0239409002
#> V493  -11.850747   0.013418669  0.086225656  23.70149       0.0818965626
#> V494   -9.300614   0.003730784  0.006773200  18.60123       0.1834706794
#> V495   -9.767585   0.004960139  0.012120351  19.53517       0.1966633445
#> V496   -9.417914   0.004080905  0.008332896  18.83583       0.1358583408
#> V497   -9.844694   0.005712470  0.015940417  19.68939      -0.0211216831
#> V498   -9.680369   0.004734657  0.011163506  19.36074       0.0317539283
#> V499   -9.287272   0.003771182  0.007024518  18.57454      -0.0374943623
#> V500   -9.465244   0.004189717  0.008651339  18.93049       0.0219274979
#> V501  -10.468563   0.007779676  0.029194582  20.93713      -0.0527895004
#> V502   -9.641820   0.004785590  0.011459176  19.28364       0.0884042768
#> V503   -9.225727   0.003008559  0.004456253  18.45145       0.1567468332
#> V504   -9.071708   0.003079118  0.004687547  18.14342       0.0490337367
#> V506   -9.131077   0.002980272  0.004398842  18.26215      -0.0877071662
#> V507  -10.592027   0.008386966  0.033489679  21.18405       0.0962806912
#> V508  -10.041619   0.006794515  0.022091435  20.08324       0.0606304796
#> V509  -11.703846   0.010158764  0.050357532  23.40769       0.2532546804
#> V510   -9.272870   0.003459683  0.006031721  18.54574      -0.0260121304
#> V511   -9.570176   0.006044240  0.017777750  19.14035       0.0253567877
#> V512  -16.036010   0.027633607  0.261805555  32.07202       0.2799851944
#> V513  -10.325076   0.007406456  0.026649702  20.65015       0.0230966399
#> V514  -10.371950   0.007757031  0.029301587  20.74390      -0.0758955035
#> V515   -9.248213   0.003874737  0.007371454  18.49643       0.0975562356
#> V516   -9.206071   0.003199281  0.005093514  18.41214       0.0958252622
#> V517   -9.755122   0.006197853  0.018384431  19.51024       0.0902180794
#> V518  -10.332578   0.007786726  0.029086796  20.66516      -0.0365578189
#> V519   -9.417449   0.004596252  0.010340555  18.83490      -0.0492255372
#> V520   -9.740393   0.005031515  0.012439635  19.48079       0.2496865543
#> V521  -15.255205   0.025800631  0.251188069  30.51041       0.1306017214
#> V522   -9.219633   0.003682081  0.006746912  18.43927       0.1268731590
#> V523   -9.329001   0.003764727  0.006967253  18.65800       0.0604282853
#> V524   -9.777667   0.006184049  0.018302955  19.55533       0.0548046803
#> V525  -10.301713   0.007796068  0.029177875  20.60343      -0.0540005279
#> V526   -9.754438   0.006193964  0.018371173  19.50888       0.0644601393
#> V527   -9.776814   0.004343127  0.009245167  19.55363       0.0402520407
#> V528   -9.211126   0.003774856  0.007026461  18.42225       0.0149750377
#> V530   -9.212541   0.003649383  0.006623953  18.42508       0.1360177932
#> V531   -9.072049   0.003106278  0.004808298  18.14410       0.0584201196
#> V532  -15.076827   0.025731327  0.242477862  30.15365       0.2623127098
#> V533   -9.324618   0.003773833  0.006970128  18.64924       0.0471544967
#> V534  -11.181612   0.010640238  0.054623748  22.36322       0.0796577098
#> V535   -9.680681   0.004960044  0.011932096  19.36136       0.0521471920
#> V536  -10.076278   0.006412702  0.020077035  20.15256      -0.0160800650
#> V539   -9.082035   0.003023792  0.004521871  18.16407       0.0023643332
#> V540   -9.738160   0.005215666  0.013299347  19.47632       0.0665125509
#> V541   -9.410217   0.003966610  0.007919289  18.82043       0.0565005216
#> V542   -9.842177   0.006161360  0.018156165  19.68435       0.0572654776
#> V543   -9.096768   0.003216493  0.005132648  18.19354      -0.0541963231
#> V544   -9.326077   0.003770810  0.006961088  18.65215       0.0408412552
#> V545   -9.353329   0.004197418  0.008734431  18.70666       0.0231682920
#> V546   -9.257141   0.003940339  0.007641216  18.51428      -0.0106426632
#> V547   -9.349985   0.004202174  0.008752402  18.69997       0.0166338224
#> V548   -9.175688   0.003057877  0.004607479  18.35138       0.0885781980
#> V549   -9.804123   0.006412858  0.019612648  19.60825       0.0361669436
#> V550   -9.074565   0.003150081  0.004911070  18.14913       0.0045874138
#> V551  -11.368742   0.008935997  0.040209352  22.73748       0.0206731910
#> V552   -9.079460   0.003136605  0.004865244  18.15892       0.0147596186
#> V553   -9.700810   0.004884751  0.011619596  19.40162       0.0513287576
#> V554   -9.141026   0.003393137  0.005737921  18.28205       0.0949721384
#> V555   -9.529001   0.004456777  0.009977697  19.05800       0.0273464275
#> V557  -10.351031   0.007376122  0.026463267  20.70206      -0.0539178482
#> V558   -9.265219   0.003935896  0.007621461  18.53044       0.0050731009
#> V559   -9.181834   0.003654736  0.006587452  18.36367       0.0011007875
#> V560   -9.412192   0.003918246  0.007534717  18.82438       0.0784435173
#> V561  -11.243061   0.009154577  0.041663271  22.48612       0.0841488970
#> V562   -9.132478   0.003439523  0.005841723  18.26496       0.0287904448
#> V563   -9.313532   0.003772916  0.007091461  18.62706       0.1016662582
#> V564   -9.754811   0.006189461  0.018346905  19.50962       0.0570516086
#> V565   -9.815059   0.006155571  0.018491732  19.63012      -0.0205918215
#> V566   -9.371305   0.004189465  0.008706140  18.74261       0.0120819066
#> V567   -9.637580   0.003880384  0.007405623  19.27516       0.0308534476
#> V568   -9.866571   0.006136516  0.018009040  19.73314       0.0645405026
#> V569   -9.803624   0.006156793  0.018147705  19.60725       0.0395143078
#> V570   -9.753544   0.006191442  0.018362346  19.50709       0.0602810530
#> V571  -12.784631   0.018630691  0.158091816  25.56926      -0.0202914095
#> V572   -9.446086   0.004587321  0.010415135  18.89217      -0.0632564460
#> V573   -9.263190   0.003960541  0.007710875  18.52638       0.0309515163
#> V574   -9.055957   0.003068969  0.004681783  18.11191       0.0472534937
#> V576   -9.216500   0.003623290  0.006504370  18.43300       0.1761840888
#> V577  -10.124848   0.007313527  0.025809568  20.24970      -0.0495062862
#> V578   -9.618921   0.005332153  0.013972538  19.23784      -0.0177567638
#> V579   -9.243342   0.003932641  0.007610587  18.48668       0.0544221975
#> V580   -9.218996   0.003628043  0.006516787  18.43799       0.1806173966
#> V581   -9.790004   0.006176827  0.018265822  19.58001       0.0319183980
#> V582   -9.462299   0.004645305  0.010674103  18.92460      -0.0548556497
#> V583   -9.293554   0.003980977  0.007799738  18.58711      -0.0343397696
#> V584  -10.228054   0.007883898  0.030002887  20.45611      -0.0573985304
#> V585  -11.529480   0.009861498  0.047919489  23.05896       0.1231458810
#> V586  -14.824193   0.020622569  0.168453115  29.64839       0.2593658265
#> V587   -9.763941   0.006188681  0.018356091  19.52788       0.0823615340
#> V588  -10.381795   0.007401938  0.026617778  20.76359      -0.0043545195
#> V589  -10.400179   0.007700043  0.028634850  20.80036      -0.0381333696
#> V590  -11.987629   0.010243240  0.052418371  23.97526       0.0078374902
#> V592  -11.925007   0.010514958  0.053931060  23.85001       0.0303000998
#> V593  -10.340757   0.007384330  0.026523796  20.68151      -0.0414300280
#> V594   -9.425791   0.003965823  0.007909925  18.85158       0.0480584021
#> V595  -10.325742   0.007377586  0.026452518  20.65148       0.0102176407
#> V596  -11.369514   0.009084834  0.041818574  22.73903      -0.0310014931
#> V598   -9.442798   0.003386569  0.005636539  18.88560       0.1715213378
#> V599   -9.832051   0.005458692  0.014501259  19.66410       0.0898559299
#> V600   -9.797542   0.005136050  0.012916191  19.59508       0.2836080826
#> V601  -10.213531   0.007837810  0.029600382  20.42706      -0.0698239945
#> V602   -9.472162   0.005056119  0.012531798  18.94432       0.0373018653
#> V603   -9.616807   0.004450943  0.009756983  19.23361       0.1318176632
#> V604   -9.315694   0.003649144  0.006521805  18.63139       0.0510951469
#> V605  -11.027717   0.009679265  0.045245171  22.05543       0.0782674920
#> V606   -9.114619   0.003118509  0.004795955  18.22924      -0.0380959536
#> V607  -11.947655   0.014099639  0.095076023  23.89531       0.0241620974
#> V608   -9.896141   0.006278299  0.018801313  19.79228       0.1072134310
#> V609   -9.422661   0.003550629  0.006298585  18.84532       0.1474831506
#> V610   -9.115382   0.003031623  0.004572243  18.23076      -0.0247288534
#> V611   -9.601089   0.004533364  0.010259956  19.20218       0.0310083792
#> V612   -9.774397   0.005018415  0.012259187  19.54879       0.0555952392
#> V613   -9.151733   0.003203023  0.005078048  18.30347      -0.1050398389
#> V615   -9.062953   0.003066536  0.004648153  18.12591       0.0523319933
#> V616 -109.017432   0.730915232 25.958255093 218.03486       3.5985125206
#> V617  -10.755465   0.008922124  0.038655342  21.51093       0.0144059395
#> V618  -12.220237   0.010285768  0.053592556  24.44047      -0.1479506624
#> V619  -11.704796   0.009272273  0.042545933  23.40959      -0.0679327169
#> V620  -13.137253   0.012049728  0.069479120  26.27451       0.0404916171
#> V621   -9.165950   0.003031013  0.004533968  18.33190       0.1013864989
#> V622   -9.138130   0.003143226  0.004857983  18.27626       0.0172245494
#> V623   -9.089358   0.003008351  0.004477138  18.17872      -0.0023521873
#> V624   -9.694683   0.005042573  0.012442124  19.38937       0.1082480113
#> V625   -9.079186   0.003025380  0.004525450  18.15837       0.0208274033
#> V626   -9.449291   0.004124001  0.008331244  18.89858       0.0293355700
#> V627   -9.128333   0.003018627  0.004511630  18.25667      -0.0989044881
#> V628   -9.140457   0.003271682  0.005252158  18.28091       0.0696459100
#> V629   -9.101449   0.003070055  0.004649972  18.20290       0.0019646970
#> V630  -11.496155   0.009029365  0.041440715  22.99231      -0.0903588649
#> V631   -9.118190   0.003170247  0.004942362  18.23638       0.0289314570
#> V632   -9.523487   0.004460202  0.009997279  19.04697       0.0315041309
#> V633   -9.077969   0.003031127  0.004543293  18.15594       0.0067125557
#> V634  -10.745033   0.007406679  0.026450552  21.49007       0.0556888904
#> V635  -16.560808   0.027443523  0.278915926  33.12162       0.0801174457
#> V636   -9.177943   0.003271151  0.005246811  18.35589       0.0378714044
#> V638   -9.777556   0.005940795  0.017272442  19.55511      -0.0283405789
#> V639   -9.243155   0.003944706  0.007661088  18.48631      -0.0119622920
#> V640   -9.258219   0.003924654  0.007587204  18.51644       0.0618535395
#> V641  -12.421822   0.017067006  0.134033446  24.84364       0.0185048640
#> V642   -9.816033   0.006288212  0.018932401  19.63207       0.1679002613
#> V643   -9.322031   0.004245884  0.008866744  18.64406      -0.0102809605
#> V645   -9.687053   0.005119085  0.012769583  19.37411       0.0286007206
#> V646   -9.703545   0.005126965  0.012808178  19.40709       0.1108540654
#> V647  -10.854807   0.010040744  0.047625645  21.70961      -0.0320614544
#> V648   -9.444294   0.004248490  0.008874195  18.88859       0.1291145670
#> V649   -9.605007   0.004631772  0.010749593  19.21001       0.0255901383
#> V650  -12.527396   0.011290465  0.060717108  25.05479       0.0112472669
#> V651   -9.239560   0.003417591  0.005726082  18.47912       0.0280093553
#> V652   -9.094570   0.003007870  0.004476415  18.18914      -0.0245681589
#> V653  -14.939596   0.018534383  0.161415636  29.87919       0.0217158963
#> V654  -10.085354   0.006664651  0.021670278  20.17071      -0.0783585237
#> V655   -9.287592   0.003372723  0.005621584  18.57518      -0.1661132284
#> V656   -9.216221   0.003497783  0.006085017  18.43244       0.1678534793
#> V657  -11.624240   0.009093797  0.041203489  23.24848      -0.0345391137
#> V658   -9.872191   0.005619220  0.015427929  19.74438       0.0726699064
#> V659   -9.399826   0.003978850  0.007844236  18.79965       0.2258333319
#> V660   -9.184129   0.003192023  0.005002262  18.36826       0.0149584886
#> V661   -9.132718   0.003082023  0.004690354  18.26544      -0.0754960332
#> V663  -10.349058   0.007378733  0.026482240  20.69812      -0.0516687719
#> V664   -9.205856   0.003383929  0.005596335  18.41171       0.0864921586
#> V666   -9.437533   0.004351861  0.009385396  18.87507      -0.0040684833
#> V667   -9.650595   0.004207646  0.008706131  19.30119       0.1247835457
#> V668   -9.820766   0.006280513  0.018823976  19.64153       0.1032327165
#> V669   -9.277934   0.003402228  0.005738806  18.55587      -0.0630564619
#> V670   -9.135754   0.003446018  0.005849279  18.27151       0.0901196396
#> V671   -9.119398   0.002988923  0.004423145  18.23880      -0.0712885780
#> V672   -9.102232   0.003035604  0.004561711  18.20446      -0.1016852553
#> V673  -10.307132   0.006793730  0.022005260  20.61426       0.0701130700
#> V674   -9.492860   0.004428185  0.009426112  18.98572       0.2014309241
#> V675   -9.148127   0.003343319  0.005592900  18.29625       0.1018330704
#> V676  -10.538063   0.007775602  0.029148371  21.07613      -0.0240197097
#> V677   -9.874851   0.006127531  0.017956848  19.74970       0.0651785399
#> V678   -9.171980   0.003497729  0.006044987  18.34396      -0.0083736026
#> V679   -9.500402   0.004833215  0.011536854  19.00080      -0.0546318905
#> V680  -11.510026   0.009349559  0.043836875  23.02005       0.0990088448
#> V681  -11.932636   0.010212134  0.052067162  23.86527       0.0394441638
#> V682   -9.144523   0.003310396  0.005449089  18.28905       0.0159530265
#> V683  -10.800275   0.008194476  0.031994445  21.60055      -0.0444065368
#> V685  -10.157213   0.006949387  0.023547770  20.31443      -0.0545121209
#> V686   -9.231949   0.003358520  0.005553996  18.46390      -0.0803152117
#> V687   -9.548783   0.004422390  0.009744334  19.09757       0.1062550266
#> V688  -15.992995   0.027518755  0.261801009  31.98599       0.2868910958
#> V689  -11.254664   0.009363987  0.041303625  22.50933       0.0872893433
#> V690  -11.385829   0.009149014  0.041314486  22.77166      -0.0348376283
#> V691  -11.696760   0.009178610  0.042183619  23.39352      -0.0922964869
#> V692   -9.861387   0.005282064  0.013623759  19.72277       0.2617788674
#> V693  -10.240729   0.007895680  0.030100077  20.48146      -0.0468132509
#> V694   -9.814500   0.006051536  0.017764228  19.62900       0.0414431419
#> V695   -9.462619   0.004595481  0.010445410  18.92524      -0.0237765363
#> V696   -9.618616   0.004659263  0.010861015  19.23723       0.0136727337
#> V697  -10.508255   0.007568272  0.027706265  21.01651      -0.0217919433
#> V698   -9.788935   0.006211861  0.018487595  19.57787       0.1286374138
#> V699   -9.993681   0.006795116  0.021907299  19.98736       0.0931520558
#> V700   -9.468291   0.004213823  0.008707927  18.93658       0.0432567591
#> V701   -9.493909   0.004287549  0.008942369  18.98782       0.1099418847
#> V702   -9.325251   0.003794156  0.007114697  18.65050       0.0008456072
#> V703  -10.098472   0.005726638  0.015936749  20.19694       0.0541882566
#> V704   -9.247124   0.003951195  0.007680477  18.49425      -0.0014023686
#> V705   -9.218904   0.003633222  0.006575844  18.43781       0.0958348313
#> V706   -9.366857   0.004197747  0.008740593  18.73371       0.0140677666
#> V707   -9.375815   0.004345595  0.009281373  18.75163       0.0730588924
#> V708   -9.254007   0.003652525  0.006658529  18.50801       0.0396061489
#> V709   -9.661734   0.004870449  0.011591227  19.32347       0.1318333824
#> V710  -10.038076   0.005183114  0.012963038  20.07615       0.1334485792
#> V711   -9.780076   0.006251249  0.018716286  19.56015       0.1281612094
#> V712   -9.411167   0.004033073  0.008055604  18.82233      -0.0042950371
#> V713   -9.378407   0.004055443  0.008113017  18.75681       0.0611111155
#> V714   -9.253485   0.003779291  0.007084969  18.50697       0.0874575798
#> V715   -9.398970   0.003949909  0.007725322  18.79794       0.0049375207
#> V716   -9.889138   0.006132943  0.017974967  19.77828       0.0932705484
#> V717   -9.150743   0.003439541  0.005830071  18.30149       0.1019835184
#> V718   -9.311663   0.003669995  0.006717990  18.62333       0.0760789603
#> V719   -9.812535   0.005089642  0.012851250  19.62507       0.0400634186
#> V720  -12.541648   0.011993234  0.070034215  25.08330       0.0727543081
#> V721   -9.588351   0.004671993  0.010927951  19.17670       0.0316386461
#> V722  -10.389906   0.007387055  0.026516813  20.77981      -0.0322295022
#> V723   -9.405536   0.004676089  0.010668702  18.81107       0.0603058471
#> V724  -15.788953   0.025794650  0.253715769  31.57791       0.0556468289
#> V725   -9.210792   0.003809105  0.007138803  18.42158       0.0420736985
#> V726   -9.171368   0.003374178  0.005685735  18.34274       0.1107678702
#> V728   -9.671168   0.004671350  0.010878285  19.34234       0.0725775988
#> V729   -9.209259   0.003554367  0.006259224  18.41852       0.1685780730
#> V730   -9.906565   0.005491057  0.014919611  19.81313       0.0409590668
#> V731   -9.243541   0.003944974  0.007654963  18.48708       0.0552252359
#> V732   -9.246973   0.003945578  0.007662544  18.49395      -0.0208705820
#> 
#> $diagnostics
#> $diagnostics$pareto_k
#>   [1]  0.2805333777  0.2798825319  0.6842019030  0.3949851814  0.0804319967
#>   [6]  0.1044244198 -0.0784228410  0.1169005633  0.1701014884  0.2632665408
#>  [11]  0.0613053786  0.1281075094 -0.0451562902  0.1642131849  0.2337729773
#>  [16] -0.2180201782  0.3986629207  0.2156351498  0.1002983996 -0.1123770110
#>  [21]  0.1484936882  0.1269803549 -0.0781622532  0.3563426528  0.2587813273
#>  [26]  0.1668654374 -0.1340372846 -0.1051602576 -0.1323835955 -0.1285920228
#>  [31] -0.0941039065  0.0837300021  0.3627345285  0.0267324679  0.0622086956
#>  [36]  0.0236363602  0.0261040656 -0.0549847466 -0.0962013255  0.0752005956
#>  [41]  0.2924563746  0.0247981090 -0.0478509166 -0.0636718474 -0.0054117413
#>  [46]  0.1604101631  0.0832472819  0.1582002285 -0.0311433121  0.0304137309
#>  [51]  0.0522401237 -0.0227890113  0.1450307609  0.0484424160  0.0661998807
#>  [56]  0.4185708151 -0.0854922568 -0.0241711073 -0.0246322754  0.0574937647
#>  [61] -0.1463233047  0.1207596331  0.0046194918 -0.1743627590 -0.0948302136
#>  [66]  0.1204100158  0.1692279982  0.1519561549 -0.1054436232 -0.1675144779
#>  [71]  0.1110184046  0.1704365351 -0.1044293986 -0.3108710849 -0.0814861310
#>  [76] -0.1483343527 -0.0667730513  0.1140641717  0.0693816366  0.0368815466
#>  [81] -0.1169148885 -0.0597296770  0.0663357224 -0.1041650067  0.3212338584
#>  [86] -0.1943427635  0.0029059434  0.0558260550  0.0335786921 -0.1215922457
#>  [91]  0.0981274980 -0.0746334183 -0.0246700002  0.0460348826  0.0670681359
#>  [96] -0.0207851553 -0.1782799903 -0.1255035154  0.0237519147  0.0679384913
#> [101]  0.1168499688  0.0381957834  0.0973203980  0.3868430964  0.0479198283
#> [106]  0.0786809141 -0.1828303202 -0.0003820443 -0.0820497994 -0.0826846517
#> [111]  0.1828673186  0.1347254807 -0.1211637032  0.2476033515 -0.0765529844
#> [116]  0.2281717968  0.2362509445  0.1766718400  0.0849041436  0.1256285979
#> [121]  0.0271105966  0.1362733616  0.0792951380  0.1076717298 -0.0486681161
#> [126] -0.0163138966  0.2046934068  0.0778503259  0.0320218425  0.1517138798
#> [131]  0.1354908627  0.2005940539  0.1729223357  0.0999930365  0.1584053306
#> [136]  0.0088337301 -0.2941180034  0.1908629646  0.3447135309  0.1172721549
#> [141]  0.2851359123 -0.0471988161  0.0911048742 -0.1552753499  0.3008410199
#> [146]  0.0031703482  0.0190834523  0.2240273580  0.1176278088  0.0396945883
#> [151] -0.1239588789  0.1647292924  0.1857553852  0.1463164476 -0.0138781001
#> [156]  0.4333877144  0.0709722505  0.1173450926 -0.0001371779 -0.1568834611
#> [161]  0.0971300663  0.0589456646  0.0222482568  0.3378672494  0.2113421999
#> [166] -0.0955657571  0.3889960179  0.0328342663 -0.1544669882 -0.1155508701
#> [171]  0.1485962899  0.0170324614  0.0579950015 -0.0692372189 -0.1162129358
#> [176]  0.0319342873  0.1894097448  0.0538888247  0.0181889305 -0.0519280013
#> [181]  0.1237276424 -0.0486912546  1.1120057878  0.2411138600  0.1198929207
#> [186]  0.0450122774  0.0489632334  0.3326752715  0.0506331409  0.0082864543
#> [191] -0.1443565516  0.0220678895  0.0858328269  0.0301340470  0.2231426589
#> [196]  0.1688319904 -0.0505029186 -0.1033426612  0.0621011592  0.1572224772
#> [201]  0.0474061865 -0.0304810997 -0.1101512672 -0.0491571479  0.0500886534
#> [206] -0.0475414631 -0.1172145979 -0.0484870415  0.2678563208  0.1121199592
#> [211] -0.0091878531  0.0634940898 -0.1194505658  0.0688889067 -0.0880703251
#> [216] -0.1296671317  0.3128650512  0.1726484374  0.1754469284  0.1715828463
#> [221]  0.0063436550  0.4268420884  0.1006224819 -0.0579619256  0.3225848193
#> [226]  0.0532717550 -0.0592390638  0.2014025505  0.0147413536 -0.0802265372
#> [231] -0.1570584625 -0.0578976911 -0.1113116507  0.2547562319 -0.0588953452
#> [236] -0.1538674513 -0.1042570428  0.0056075354 -0.1248211949  0.0928217518
#> [241] -0.0590291043  0.1190651704 -0.0314569009  0.1954137246 -0.1426470532
#> [246] -0.1114672029 -0.1245202263 -0.1364529135  0.1582425417  0.0681564583
#> [251]  0.0542303840  0.2433232375  0.0709814303  0.0732142630  0.1806243116
#> [256]  0.0457680548  0.0950844229 -0.0780743242  0.1016593364  0.0821952532
#> [261]  0.0655206993  0.0296837065  0.1449446555 -0.0599375196  0.1088705500
#> [266]  0.1034561798 -0.0097287123  0.1015928829 -0.1659320872  0.1475228442
#> [271]  0.1564023251  0.1567427394 -0.0531946033  0.1007910082 -0.1156180892
#> [276]  0.0510568621  0.0825192700  0.1115806100  0.0074122014  0.1448711303
#> [281]  0.1876937019 -0.0834683640  0.1292260175  0.0569002304 -0.0076087623
#> [286]  0.0849168704  0.1788158278 -0.0054972864  0.1739759225 -0.1120367077
#> [291]  0.1313413900  0.0331667401  0.0315516828  0.0702144616  0.1132014474
#> [296]  0.0715064739  0.0681046606 -0.0346584832 -0.1167949741  0.0236632585
#> [301]  0.1551879425 -0.0014398791  0.0525650847  0.1015315934  0.0998281643
#> [306]  0.0921396010  0.0375345713  0.0707023472 -0.0266916825  0.1273820643
#> [311]  0.1358734332  0.0775642789  0.1142094078  0.1021191237 -0.1999852214
#> [316]  0.0699775851  0.0206669551  0.0166143137  0.9776573253  0.1200897408
#> [321]  0.0657972685  0.1147409746 -0.0383839630  0.0435603486 -0.0762910118
#> [326]  0.2345514864  0.0026725119  0.1192217066  0.1251352428 -0.1303622447
#> [331] -0.0660220284  0.1161042073  0.2051076517  0.0744559521  0.0738769662
#> [336]  0.2689155232  0.2638841057  0.0872861699 -0.0494513534 -0.0012265040
#> [341]  0.1831465400  0.1811045010  0.0009528022  0.2423494355  0.2352268441
#> [346] -0.0660421640  0.2001702460  0.1968491073  0.0736398718  0.2503809925
#> [351]  0.2744156824 -0.2145835357 -0.1920492739 -0.1513661079  0.1099437721
#> [356]  0.1524923120  0.1313017341 -0.1576484562 -0.1367438579  0.1931438278
#> [361] -0.2993539889  0.1836926476 -0.1135863116  0.1873193468 -0.0796331437
#> [366]  0.1039146172  0.1367063512  0.1578869870  0.1054930672  0.1717371704
#> [371]  0.0552662039  0.3197337529 -0.1006443333  0.0517736270 -0.0393945869
#> [376]  0.1823250102  0.1868590544 -0.0213122271  0.2336018814  0.1144355340
#> [381]  0.1548690381 -0.1291773284 -0.0987003231  0.0568135255 -0.0382827739
#> [386] -0.1024059935 -0.0403553404  0.1423544454 -0.0587598412  0.0097263783
#> [391] -0.1809447299 -0.0362733785  0.0716808053  0.2839127676 -0.1137876936
#> [396] -0.0316493838  0.0833028923  0.2171338252  0.0670720132 -0.0013233432
#> [401]  0.0983013335 -0.0674521636 -0.0209751451  0.0508739569  0.0479799661
#> [406]  0.0288221153  0.0006679855 -0.1946406420 -0.0919918690  0.1162073396
#> [411]  0.0452538418  0.1349479450  0.0342864105 -0.0088882709 -0.1112715602
#> [416]  0.1443411647  0.1918237272  0.1787660055  0.2358305844  0.0707271333
#> [421]  0.0428504316  0.0241951667 -0.0130871949  0.0648240863  0.0609796460
#> [426] -0.0478392554 -0.0673069177 -0.0067395588  0.0395423664 -0.0005102754
#> [431]  0.0803043654  0.2443770451  0.0305322852 -0.0652615347  0.0870622976
#> [436]  0.3140729831 -0.0416206634  0.1056485624  0.1768318406 -0.1035285292
#> [441]  0.0847414961  0.1168810137  0.1884031312  0.0363794076 -0.0389227394
#> [446] -0.0728630292  0.1101451980  0.0875427947  0.0977772123  0.0047488296
#> [451]  0.0017236815  0.0255844892  0.1917546804  0.0239409002  0.0818965626
#> [456]  0.1834706794  0.1966633445  0.1358583408 -0.0211216831  0.0317539283
#> [461] -0.0374943623  0.0219274979 -0.0527895004  0.0884042768  0.1567468332
#> [466]  0.0490337367 -0.0877071662  0.0962806912  0.0606304796  0.2532546804
#> [471] -0.0260121304  0.0253567877  0.2799851944  0.0230966399 -0.0758955035
#> [476]  0.0975562356  0.0958252622  0.0902180794 -0.0365578189 -0.0492255372
#> [481]  0.2496865543  0.1306017214  0.1268731590  0.0604282853  0.0548046803
#> [486] -0.0540005279  0.0644601393  0.0402520407  0.0149750377  0.1360177932
#> [491]  0.0584201196  0.2623127098  0.0471544967  0.0796577098  0.0521471920
#> [496] -0.0160800650  0.0023643332  0.0665125509  0.0565005216  0.0572654776
#> [501] -0.0541963231  0.0408412552  0.0231682920 -0.0106426632  0.0166338224
#> [506]  0.0885781980  0.0361669436  0.0045874138  0.0206731910  0.0147596186
#> [511]  0.0513287576  0.0949721384  0.0273464275 -0.0539178482  0.0050731009
#> [516]  0.0011007875  0.0784435173  0.0841488970  0.0287904448  0.1016662582
#> [521]  0.0570516086 -0.0205918215  0.0120819066  0.0308534476  0.0645405026
#> [526]  0.0395143078  0.0602810530 -0.0202914095 -0.0632564460  0.0309515163
#> [531]  0.0472534937  0.1761840888 -0.0495062862 -0.0177567638  0.0544221975
#> [536]  0.1806173966  0.0319183980 -0.0548556497 -0.0343397696 -0.0573985304
#> [541]  0.1231458810  0.2593658265  0.0823615340 -0.0043545195 -0.0381333696
#> [546]  0.0078374902  0.0303000998 -0.0414300280  0.0480584021  0.0102176407
#> [551] -0.0310014931  0.1715213378  0.0898559299  0.2836080826 -0.0698239945
#> [556]  0.0373018653  0.1318176632  0.0510951469  0.0782674920 -0.0380959536
#> [561]  0.0241620974  0.1072134310  0.1474831506 -0.0247288534  0.0310083792
#> [566]  0.0555952392 -0.1050398389  0.0523319933  3.5985125206  0.0144059395
#> [571] -0.1479506624 -0.0679327169  0.0404916171  0.1013864989  0.0172245494
#> [576] -0.0023521873  0.1082480113  0.0208274033  0.0293355700 -0.0989044881
#> [581]  0.0696459100  0.0019646970 -0.0903588649  0.0289314570  0.0315041309
#> [586]  0.0067125557  0.0556888904  0.0801174457  0.0378714044 -0.0283405789
#> [591] -0.0119622920  0.0618535395  0.0185048640  0.1679002613 -0.0102809605
#> [596]  0.0286007206  0.1108540654 -0.0320614544  0.1291145670  0.0255901383
#> [601]  0.0112472669  0.0280093553 -0.0245681589  0.0217158963 -0.0783585237
#> [606] -0.1661132284  0.1678534793 -0.0345391137  0.0726699064  0.2258333319
#> [611]  0.0149584886 -0.0754960332 -0.0516687719  0.0864921586 -0.0040684833
#> [616]  0.1247835457  0.1032327165 -0.0630564619  0.0901196396 -0.0712885780
#> [621] -0.1016852553  0.0701130700  0.2014309241  0.1018330704 -0.0240197097
#> [626]  0.0651785399 -0.0083736026 -0.0546318905  0.0990088448  0.0394441638
#> [631]  0.0159530265 -0.0444065368 -0.0545121209 -0.0803152117  0.1062550266
#> [636]  0.2868910958  0.0872893433 -0.0348376283 -0.0922964869  0.2617788674
#> [641] -0.0468132509  0.0414431419 -0.0237765363  0.0136727337 -0.0217919433
#> [646]  0.1286374138  0.0931520558  0.0432567591  0.1099418847  0.0008456072
#> [651]  0.0541882566 -0.0014023686  0.0958348313  0.0140677666  0.0730588924
#> [656]  0.0396061489  0.1318333824  0.1334485792  0.1281612094 -0.0042950371
#> [661]  0.0611111155  0.0874575798  0.0049375207  0.0932705484  0.1019835184
#> [666]  0.0760789603  0.0400634186  0.0727543081  0.0316386461 -0.0322295022
#> [671]  0.0603058471  0.0556468289  0.0420736985  0.1107678702  0.0725775988
#> [676]  0.1685780730  0.0409590668  0.0552252359 -0.0208705820
#> 
#> $diagnostics$n_eff
#>   [1] 393.093837 494.573270  55.801272 400.238078 492.051777 488.335067
#>   [7] 496.946812 495.433060 493.693768 495.140890 494.468953 466.927609
#>  [13] 496.640966 488.890770 495.070771 497.822323 398.813926 495.627561
#>  [19] 481.730376 485.033070 478.285233 495.950776 485.488291 347.915347
#>  [25] 495.846952 412.040174 492.636720 492.734449 492.521554 489.652263
#>  [31] 492.878894 491.993910 341.922875 495.939655 482.764283 496.199402
#>  [37] 494.326021 492.752989 492.790466 497.449114 403.503453 497.582978
#>  [43] 492.788741 492.277882 497.132267 482.002658 494.533972 494.434866
#>  [49] 463.103190 461.690077 496.514330 497.535620 496.372185 496.640152
#>  [55] 497.321840 338.401331 496.032153 497.469758 496.357478 495.824845
#>  [61] 495.186701 497.668237 496.437102 495.185871 490.085030 495.714589
#>  [67] 427.519980 495.901057 492.773449 497.404806 497.508107 497.721058
#>  [73] 497.283516 496.633088 496.122384 496.894674 493.155039 497.540349
#>  [79] 497.049815 496.534091 490.730857 480.896090 497.749670 496.994153
#>  [85] 402.594835 497.259134 478.778020 496.087864 497.573488 496.569472
#>  [91] 497.783478 494.277946 497.321064 473.007482 473.903301 496.258623
#>  [97] 497.664593 492.355104 463.481018 476.783482 496.738189 497.680267
#> [103] 489.266088 405.739793 496.641029 489.370992 494.401913 496.427922
#> [109] 492.819006 493.314785 489.647780 495.002800 495.454509 490.251184
#> [115] 488.597705 485.018771 495.577665 489.098887 485.499751 488.336743
#> [121] 456.145821 489.391818 481.522613 489.290806 497.783763 497.701993
#> [127] 495.327424 494.413109 497.660506 495.901448 490.059222 491.204807
#> [133] 495.906763 478.616952 495.161569 497.312094 497.807374 494.903666
#> [139] 362.086435 489.395610 402.773300 495.849573 497.580596 497.536773
#> [145] 393.359834 491.352575 474.978979 402.240732 490.651337 484.642192
#> [151] 492.135818 494.902871 495.281893 495.760007 497.361677 325.629088
#> [157] 489.419833 495.189346 497.403362 497.705152 496.246691 497.228002
#> [163] 497.588362 399.326497 495.672239 497.544652 401.182864 474.153710
#> [169] 494.528177 492.304704 488.992986 496.450189 488.566830 497.532745
#> [175] 494.064070 497.518203 490.909594 485.502360 497.439343 491.351943
#> [181] 453.566233 497.606305   7.468004 463.428600 411.261230 497.417088
#> [187] 493.514481 367.695024 487.381307 474.910731 494.543874 495.884039
#> [193] 497.109002 496.650108 495.221371 489.259747 491.944802 491.627914
#> [199] 493.362836 489.226210 495.331122 496.247705 492.682766 496.304555
#> [205] 496.642367 491.857265 497.331177 495.891256 494.407130 489.594169
#> [211] 495.656846 492.810574 492.521654 491.979159 489.844534 492.517942
#> [217] 393.459640 495.073460 495.048709 489.417455 496.040332 364.376569
#> [223] 497.459885 496.064024 393.242813 497.335499 492.101645 495.186845
#> [229] 497.290645 496.531716 496.540312 496.117891 495.366439 480.571097
#> [235] 496.176689 496.687836 496.654760 497.309590 496.666770 492.300157
#> [241] 495.401765 495.476330 496.570714 265.425412 497.313567 497.270427
#> [247] 486.580254 482.990164 495.608726 494.300115 494.464666 485.252418
#> [253] 496.028899 496.385799 492.372121 484.013301 492.940742 495.608606
#> [259] 496.345409 491.723416 492.944088 497.425572 492.635385 496.835159
#> [265] 492.537947 496.308343 496.942244 496.235016 497.384306 467.599867
#> [271] 470.257344 426.584428 496.098451 496.907866 497.355248 493.302420
#> [277] 496.402560 321.048632 496.059895 496.233380 495.947588 497.111918
#> [283] 496.020942 496.516161 478.379745 496.364162 495.977230 496.853126
#> [289] 494.894484 463.496304 495.206657 496.625362 496.648274 493.077724
#> [295] 475.185775 495.994829 496.023691 496.206433 496.072668 495.707855
#> [301] 495.089940 493.461986 493.099362 492.409755 377.046544 486.013274
#> [307] 497.097062 496.394964 490.470683 491.500985 490.683510 496.358739
#> [313] 496.353426 492.745303 494.008222 487.988122 208.954557 477.465036
#> [319]  15.549749 492.733224 494.588942 486.331042 468.750266 493.645218
#> [325] 495.460952 487.770611 480.210579 419.041102 490.866348 496.567637
#> [331] 496.953877 493.477477 492.600174 496.868897 493.837580 488.728662
#> [337] 487.710945 491.491588 462.537879 496.701605 487.738494 487.589686
#> [343] 497.204643 487.318219 466.114807 497.513872 443.398975 459.176983
#> [349] 471.443319 470.289477 345.418209 479.735666 496.851316 497.040174
#> [355] 461.197438 494.784037 494.826049 497.498407 496.408603 494.893386
#> [361] 497.126710 487.773404 497.537363 487.614667 497.064860 491.150050
#> [367] 492.628306 492.388134 406.412316 491.602144 496.518694 481.508613
#> [373] 488.235934 485.782019 490.887757 487.734849 446.716878 495.699409
#> [379] 484.813423 487.511723 495.023545 497.443948 496.088526 495.484943
#> [385] 497.446039 496.842985 496.361899 422.183587 493.787590 496.046875
#> [391] 496.366156 496.086357 496.337713 479.059773 497.161358 496.451445
#> [397] 486.938191 481.076080 484.795550 493.418527 496.351007 496.237096
#> [403] 476.624361 480.544354 496.700964 496.662505 468.808901 496.256448
#> [409] 497.354118 494.605260 419.167516 477.973242 493.655773 486.345695
#> [415] 495.059147 496.211537 487.762195 487.266000 492.197542 495.953296
#> [421] 494.565943 495.875031 491.867008 485.034934 496.033022 494.362604
#> [427] 496.019852 497.000132 496.928947 379.665013 492.701895 490.923436
#> [433] 494.985184 494.340034 492.730300 489.297204 497.087393 496.342445
#> [439] 496.088888 497.445017 479.187572 494.834200 487.787770 495.890713
#> [445] 496.243454 496.685715 494.588868 496.287666 495.654959 497.197606
#> [451] 494.750529 497.546992 495.869244 494.423805 458.244699 496.543596
#> [457] 493.870646 495.826794 491.825991 494.353288 496.376297 495.603875
#> [463] 485.071290 494.243351 497.695075 497.583927 497.727998 482.940299
#> [469] 488.546506 475.430828 496.970052 490.839541 361.366273 486.533522
#> [475] 485.137647 496.236241 497.383598 490.431174 485.018510 494.642439
#> [481] 493.709042 373.803816 496.584638 496.437001 490.473705 484.969052
#> [487] 490.425382 495.220287 496.401998 496.655054 497.534858 375.516351
#> [493] 496.425024 472.968840 493.841933 489.740534 497.673104 493.197140
#> [499] 496.043886 490.530223 497.362263 496.428241 495.587359 496.088672
#> [505] 495.578728 497.610332 489.725667 497.490759 480.618729 497.520025
#> [511] 493.987143 497.064347 495.001131 486.553215 496.100742 496.635121
#> [517] 496.137365 479.704730 497.018199 496.392629 490.437381 490.517175
#> [523] 495.591607 496.184447 490.609330 490.542449 490.423483 424.840061
#> [529] 494.687059 496.057725 497.604516 496.741596 486.709732 492.850836
#> [535] 496.105891 496.737697 490.474802 494.554184 495.986963 484.605596
#> [541] 476.576015 411.909233 490.417920 486.499508 485.398476 474.665623
#> [547] 473.456138 486.538883 496.038460 486.610829 479.909727 497.129187
#> [553] 492.539827 493.459849 484.779876 493.566658 494.966815 496.660122
#> [559] 477.456607 497.524394 454.263245 490.190608 496.829119 497.654994
#> [565] 494.817105 493.701647 497.362227 497.609199   1.411825 480.636236
#> [571] 474.425999 479.143405 465.538909 497.653515 497.504479 497.694474
#> [577] 493.644182 497.673064 495.721047 497.670070 497.278436 497.603593
#> [583] 480.146228 497.454997 494.993446 497.662616 486.475835 362.002432
#> [589] 497.310422 491.142736 496.076094 496.107042 435.526304 490.120138
#> [595] 495.415660 493.429553 493.417313 475.730880 495.469132 494.609754
#> [601] 469.420900 497.064742 497.693779 425.807552 488.918186 497.106472
#> [607] 496.934875 479.938240 492.118246 496.050075 497.436275 497.574195
#> [613] 486.546443 497.130143 495.225793 495.489426 490.191837 497.044092
#> [619] 497.023201 497.716882 497.642503 488.508610 495.081118 497.206554
#> [625] 485.111363 490.636373 496.895667 494.105010 478.877271 474.867576
#> [631] 497.205991 483.481339 488.002540 497.138958 495.090274 362.032640
#> [637] 478.877078 479.636898 479.519400 493.083228 484.569199 490.817440
#> [643] 494.676307 494.539902 485.870515 490.354825 488.519214 495.532862
#> [649] 495.405803 496.350954 491.795972 496.071775 496.679997 495.575390
#> [655] 495.198889 496.617770 493.983384 493.341292 490.222839 495.889148
#> [661] 495.862221 496.445258 496.061610 490.638478 497.031735 496.584208
#> [667] 493.501048 466.177097 494.521611 486.523362 494.507221 373.370646
#> [673] 496.357417 497.111270 494.524276 496.852605 492.482207 496.085134
#> [679] 496.077552
#> 
#> $diagnostics$r_eff
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [112] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [149] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [186] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [223] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [260] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [297] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [334] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [371] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [408] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [445] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [482] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [519] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [556] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [593] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [630] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [667] 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
#> 
#> $psis_object
#> NULL
#> 

# Show density plots for all parameters
diag.den <- diagnostic(model.long, type = "denplot", param = "alpha")

# Plots of imputations for all effect data
p1 <- plot(model.long, class = "scatter", outcome = "effects", time.plot = "all")

# Summarise the CEA results
summary(model.long)
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under MAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.489 0.022 0.449 0.531
#> 2 0.612 0.023 0.564 0.656
#> 
#>  
#>  Mean costs by intervention 
#>   Mean  SD   QL   QU
#> 1 5395 343 4746 6078
#> 2 5834 352 5173 6526
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 19042 1129 16809 21224
#> 2 24779 1207 22298 27224

# }
#
#
```
