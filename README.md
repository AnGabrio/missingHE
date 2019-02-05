
missingHE [![Travis-CI Build Status](https://travis-ci.org/AnGabrio/missingHE.svg?branch=master)](https://travis-ci.org/AnGabrio/missingHE)[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/AnGabrio/missingHE?branch=master&svg=true)](https://ci.appveyor.com/project/AnGabrio/missingHE)
===========================================================================================================================================

Missing outcome data in health economic evaluation
--------------------------------------------------

Contains a suite of functions to perform economic evaluations with missing outcome (either costs, effects or both) under a range of alternative assumptions about the missing data mechanisms. Estimation of the key parameters of interest and imputation of the missing data are carried out through Markov Chain Monte Carlo (MCMC) methods and a set of pre-defined Bayesian parametric models written in the BUGS language using the software [JAGS](http://mcmc-jags.sourceforge.net/), which is called via the R package [R2jags](https://cran.r-project.org/web/packages/R2jags). In addition, missingHE provides a set of specialised functions to assess model convergence, fit to the data, plotting of observed and imputed data, and compute different measures to summarise the cost-effectvieness results.

Installation
------------

It is possible to install `missingHE` using the "development" version - this will usually be updated frequently and may be continuously tested. On Windows machines, you need to install a few dependencies, including [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first, e.g. by running

``` r
pkgs <- c("R2jags","ggplot2","gridExtra","BCEA","ggmcmc","loo","Rtools","devtools")
repos <- c("https://cran.rstudio.com") 
install.packages(pkgs,repos=repos,dependencies = "Depends")
```

before installing the package using `devtools`:

``` r
devtools::install_github("AnGabrio/missingHE")
```
