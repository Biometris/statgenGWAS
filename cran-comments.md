Release with some extra options, small bug fixes and significantly improved documentation
* New options in plot functions.
* Small bug fixes.
* Improved vignette and function documentation + extra examples.

## Test environments
* local Windows 10 install, R 4.1.0
* Ubuntu (on github actions, devel and release)
* macOS (on github actions, release)
* R-hub (devel and release)

## R CMD check results

0 errors | 0 warnings | 2 notes

installed size is  8.9Mb
  sub-directories of 1Mb or more:
    data   2.7Mb
    libs   5.2Mb

* Sizes vary per platform tested on. 

On some of the r-hub platforms an additional note:
Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1534/genetics.107.080101
    From: man/runSingleTraitGwas.Rd
          inst/doc/GWAS.html
    Status: 403
    Message: Forbidden
  URL: https://doi.org/10.1534/genetics.113.159731
    From: man/runSingleTraitGwas.Rd
          inst/doc/GWAS.html
    Status: 403
    Message: Forbidden
  URL: https://doi.org/10.1534/genetics.116.193987
    From: man/runSingleTraitGwas.Rd
          inst/doc/GWAS.html
    Status: 403
    Message: Forbidden

    
* All three URLs work fine in a web browser.

