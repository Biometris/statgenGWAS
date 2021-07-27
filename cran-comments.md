Patch release fixing test failures on some platforms. 

## Test environments
* local Windows 10 install, R 4.1.0
* Ubuntu (on github actions, devel and release)
* macOS (on github actions, release)
* R-hub (devel and release)
* R-hub patched-solaris
* R-hub noLD

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

