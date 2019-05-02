
missingHE [![Travis-CI Build Status](https://travis-ci.org/AnGabrio/missingHE.svg?branch=master)](https://travis-ci.org/AnGabrio/missingHE)[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/AnGabrio/missingHE?branch=master&svg=true)](https://ci.appveyor.com/project/AnGabrio/missingHE)
===========================================================================================================================================

Missing Outcome Data in Health Economic Evaluation
--------------------------------------------------

Contains a suite of functions for health economic evaluations with missing outcome data. The package can fit different types of statistical models under a fully Bayesian approach using Markov Chain Monte Carlo (MCMC) methods. Three classes of models can be fitted under a variety of missing data assumptions: selection models, pattern mixture models and hurdle models. In addition to model fitting, `missingHE` provides a set of specialised functions to assess model convergence and summarise the statistical and economic results using different types of measures and graphs. 

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
All models implemented in `missingHE` are written in the BUGS language using the software [JAGS](http://mcmc-jags.sourceforge.net/), which needs to be installed from its own repository and instructions for installations under different OS can be found online. Once installed, the software is called in `missingHE` via the R package [R2jags](https://cran.r-project.org/package=R2jags).
Note that the `missingHE` package is currently under active development and therefore it is advisable to reinstall the package directly from GitHub before each use to ensure that you are using the most updated version.