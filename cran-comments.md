## Test environments
* local Windows 10, R 3.6.2
* debian (on gitlab-ci), R 3.6.1
* r-hub (devel and release)

## R CMD check results

0 errors | 0 warnings | 3 notes

New submission

* This is a new release.

Possibly mis-spelled words in DESCRIPTION:
  Biometris (50:49)
  GWAS (52:72)
  statgen (50:18)
  statgenGWAS (51:68)

* These are spelled correctly.

installed size is  6.0Mb
  sub-directories of 1Mb or more:
    data   3.0Mb
    libs   2.1Mb

* Sizes vary per platform tested on. 

On some of the r-hub platforms an additional note:
Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1534/genetics.107.080101 (moved to https://www.genetics.org/lookup/doi/10.1534/genetics.107.080101)
    From: inst/doc/GWAS.html
    Status: 403
    Message: Forbidden  
  URL: https://doi.org/10.1534/genetics.113.159731 (moved to https://www.genetics.org/lookup/doi/10.1534/genetics.113.159731)
    From: inst/doc/GWAS.html
    Status: 403
    Message: Forbidden
    
* Both URLs work fine in a web browser.

