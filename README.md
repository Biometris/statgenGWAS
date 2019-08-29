# statgenGWAS

[![pipeline status](https://git.wur.nl/statistical-genetic-pipeline/statgenGWAS/badges/master/pipeline.svg)](https://git.wur.nl/statistical-genetic-pipeline/statgenGWAS/commits/master)
[![coverage report](https://git.wur.nl/statistical-genetic-pipeline/statgenGWAS/badges/master/coverage.svg)](https://git.wur.nl/statistical-genetic-pipeline/statgenGWAS/commits/master)

R Package
* recoding and imputing markers
* single trait GWAS

# Installation

For direct installation from gitlab use the following code:

``` r
## Replace the location for public and privatekey with your own.
creds <- git2r::cred_ssh_key(publickey = "M:\\.ssh\\id_rsa.pub",
                             privatekey = "M:\\.ssh\\id_rsa")
devtools::install_git(url = "git@git.wur.nl:statistical-genetic-pipeline/statgenGWAS/gwas.git",
                      credentials = creds)

```
