# missingHE 1.6.2

* Creation of the new S3 method function `postsum` to extract and display posterior summaries about model parameters filtered based on some optional arguments.

* Update of the `print.missingHE` function to return more standardised output about model coefficient estimates following usual practice for `print` functions in `R`.

* Code for some functions slightly updated to remove specification of deprecated arguments no longer valid in the most recent versions of some dependence packages' functions. 

* Update of the output generated using the `print` function.

* The package vignettes have been updated.

* Updated the version of the package in the `DESCRIPTION` file to `1.6.2`.

# missingHE 1.6.1
1 June 2026

* Bug fixes that prevented the correct implementation of the models depending on specific combination of distributional assumptions.

* Code for some functions slightly updated to ensure coherence among prior values assumed across different model specifications.

* Update of some default prior values to cover a less broad range of values.

* The package vignettes have been updated.

* Updated the version of the package in the `DESCRIPTION` file to `1.6.1`.

# missingHE 1.6
19 March 2026

* A few package functions have been updated to allow specification of models for multiple treatment arm trials and produce output which allows a higher degree of customisation, such as generating graphical and numeric results by treatment, time point and for different types of estimates of interest. 

* Key functions that were updated include: `selection`, `pattern`, `hurdle` and `lmdm` (previously called `long_miss`).

* Code for many functions dedicated to the fitting of models and extracting model results has been updated and more efficient. 

* Key functions that were updated include: `ppc`, `pic`, `plot.missingHE`, `coef.missingHE` and `summary.missingHE` 

* New customisation options have been incorporated into the package functions to allow a wider choice  in terms of output generation and display.

* The package vignettes have been updated and extended to provide a more general and detailed presentation of the main package functions using the built-in data sets as motivating cases.

* The package built-in data sets (`MenSS.RData` and `PBS.RData`) have been updated to align variable names and types to match the requirements of the inputs accepted by the updated functions.

* Function used to fit longitudinal model has been  removed and replaced with a new one (`lmdm`) to provide higher degree of flexibility and more customisation options for the user when fitting the model, including assumptions about time dependence, missingness patterns and specification of priors.

* Available choices in terms of prior distribution customisation and modification have been increased  for all models through a more intuitive approach to prior specification with examples provided in the package vignettes.

* New types of posterior predictive checks (`ppc`), information criteria (`pic`), imputed data graphs (`plot.missingHE`) and diagnostics (`diagnostic`) have been added to the pool of available choices in the corresponding functions. 

* The generation of the quantities/graphs has been updated through new function arguments and customised based on the type of model fitted.

* Updated the version of the package in the `DESCRIPTION` file to `1.6.0`.

# missingHE 1.5.1
3 July 2025

# missingHE 1.5
21 March 2023

# missingHE 1.4.1
25 June 2020

# missingHE 1.4
29 April 2020

# missingHE 1.3.2
8 January 2020

# missingHE 1.2.1
21 September 2019

# missingHE 1.1.1
23 June 2019

# missingHE 1.0.1
6 May 2019

# missingHE 1.0
3 May 2019
