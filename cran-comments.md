## Test environments
* Windows 10 64-bit install, R 3.5.2 and R 3.5.3 (local)
* Ubuntu 14.04, R 3.5.2 (on travis-ci)
* Ubuntu Linux 16.04 LTS, R 3.5.3 (r-hub)
* Win-builder (devel and release)
* Windows Server 2012 R2 (x64/x64, x86/i386), R 3.5.2 (AppVeyor)
* Windows Server 2008 R2 SP1, R 3.5.3, 32/64 bit (r-hub)
* MacOS 10.11, R 3.5.3 (local)
* check_for_cran (r-hub)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies
* There are currently no downstream dependencies for this package 

* This is my first submission.

* As per email correspondence with Martina Schmirl, I have updated the package DESCRIPTION and write package, software and API names in single quotes in title and description.

* Using check_for_cran (r-hub) I still get one note for possibly mis-spelled words in DESCRIPTION which is only related to two proper names of authors listed in the description field (Gabrio and Molenberghs). 
* I have checked that these are correctly spelled.

* I have included examples for the main functions in the package. All functions without an example are only for internal use and therefore no example is applicable. 
* I have used quick examples to ensure that they run in <5 sec. Additional and longer examples are included using \donttest. If these are run as well then the total running time can be >5 sec.

* I have clarified in the description field that the software JAGS needs to be installed as all models are fitted using this Bayesian software. 
* If JAGS is installed, then the software is automatically loaded via the R2jags package (included in imports)



