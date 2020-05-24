
missingHE [![Travis-CI Build Status](https://travis-ci.org/AnGabrio/missingHE.svg?branch=master)](https://travis-ci.org/AnGabrio/missingHE)[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/AnGabrio/missingHE?branch=master&svg=true)](https://ci.appveyor.com/project/AnGabrio/missingHE)[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/missingHE)](https://cran.r-project.org/package=missingHE)[![CRAN_Download_Badge](http://cranlogs.r-pkg.org/badges/missingHE)](https://cran.r-project.org/package=missingHE)[![CRAN_Download_Badge](http://cranlogs.r-pkg.org/badges/grand-total/missingHE?color=orange)](http://cranlogs.r-pkg.org/badges/grand-total/missingHE?color=orange)
===========================================================================================================================================

Missing Outcome Data in Health Economic Evaluation
--------------------------------------------------

Contains a suite of functions for health economic evaluations with missing outcome data. The package can fit different types of statistical models under a fully Bayesian approach using Markov Chain Monte Carlo (MCMC) methods. Three classes of models can be fitted under a variety of missing data assumptions: selection models, pattern mixture models and hurdle models. In addition to model fitting, `missingHE` provides a set of specialised functions to assess model convergence and summarise the statistical and economic results using different types of measures and graphs. 

Installation
------------

There are two ways of installing `missingHE`. A "stable" version is packaged and binary files are available for Windows and as source. To install the stable version on a Windows machine, run the following command
```R
install.packages("missingHE")
```
which installs the package from a [CRAN mirror](https://cran.r-project.org/index.html) and ensures that `install.packages()` can also install the "dependencies" (e.g. other packages that are required for `missingHE` to work).

It is also possible to install `missingHE` using the "development" version - this will usually be updated frequently and may be continuously tested. On Windows machines, you need to install a few dependencies, including [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first, e.g. by running

``` r
pkgs <- c("R2jags","ggplot2","gridExtra","BCEA","ggmcmc","loo","Rtools","devtools", "utils")
repos <- c("https://cran.rstudio.com") 
install.packages(pkgs,repos=repos,dependencies = "Depends")
```

before installing the package using `devtools`:

``` r
devtools::install_github("AnGabrio/missingHE", build_vignettes = TRUE)
```
The optional argument `build_vignettes = TRUE` allows to install the vignettes of the package locally on your computer. These consist in brief tutorials to guid the user on how to use and customise the models in `missingHE` using different functions of the package. Once the package is installed, they can be accessed by using the command

``` r
utils::browseVignettes(package = "missingHE")
```
which shows all the vignettes available for the package.

All models implemented in `missingHE` are written in the BUGS language using the software [JAGS](http://mcmc-jags.sourceforge.net/), which needs to be installed from its own repository and instructions for installations under different OS can be found online. Once installed, the software is called in `missingHE` via the R package [R2jags](https://cran.r-project.org/package=R2jags).
Note that the `missingHE` package is currently under active development and therefore it is advisable to reinstall the package directly from GitHub before each use to ensure that you are using the most updated version.