Release with some extra options, small bug fixes and significantly improved documentation. 
* New options in plot functions.
* Small bug fixes, including error on cran check results.
* Improved vignette and function documentation + extra examples.

After previous attempt to release, fixed duplicated doi.org in link

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
  URL: https://doi.org/10.1214/09-sts307
    From: man/kinship.Rd
          man/runSingleTraitGwas.Rd
          inst/doc/GWAS.html

    
* This URL works fine in a web browser.

