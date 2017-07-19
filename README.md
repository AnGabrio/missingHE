# missingHE
Health Economic Evaluations with Missing Data using a set of pre-defined Bayesian models written in BUGS 

#Installation
Contains a suite of functions to systematise the workflow involving survival analysis in health economic evaluation. survHE can fit a large range of survival models using both a frequentist approach (by calling the R package flexsurv) and a Bayesian perspective. For a selected range of models, both Integrated Nested Laplace Integration (via the R package INLA) and Hamiltonian Monte Carlo (via the R package rstan) are possible. HMC models are pre-compiled so that they can run in a very efficient and fast way. In addition to model fitting, survHE provides a set of specialised functions, for example to perform Probabilistic Sensitivity Analysis, export the results of the modelling to a spreadsheet, plotting survival curves and uncertainty around the mean estimates
