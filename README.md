## Health Economic Evaluations with Missing Data
Health Economic Evaluations with Missing Data using a set of pre-defined Bayesian models written in the BUGS language using either [JAGS](http://mcmc-jags.sourceforge.net/) or [OpenBUGS] (http://www.openbugs.net/w/FrontPage). 

## Installation
It is possible to install `missingHE` using the "development" version - this will usually be updated frequently and may be continuously tested. On Windows machines, you need to install a few dependencies, including [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first, e.g. by running
```R
pkgs <- c("R2jags","R2OpenBUGS","ggplot2","gridExtra","BCEA","ggmcmc","Rtools","devtools")
repos <- c("https://cran.rstudio.com") 
install.packages(pkgs,repos=repos,dependencies = "Depends")
```
before installing the package using `devtools`:
```R
devtools::install_github("AnGabrio/missingHE")
```
