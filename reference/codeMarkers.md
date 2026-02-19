# Code and impute markers

`codeMarkers` codes markers in a `gData` object and optionally performs
imputation of missing values as well.  
The function performs the following steps:  

1.  replace strings in `naStrings` by `NA`.

2.  remove genotypes with a fraction of missing values higher than
    `nMissGeno`.

3.  remove SNPs with a fraction of missing values higher than `nMiss`.

4.  recode SNPs to numerical values.

5.  remove SNPs with a minor allele frequency lower than `MAF`.

6.  optionally remove duplicate SNPs.

7.  optionally impute missing values.

8.  repeat steps 5. and 6. if missing values are imputed.

## Usage

``` r
codeMarkers(
  gData,
  refAll = "minor",
  nMissGeno = 1,
  nMiss = 1,
  MAF = NULL,
  MAC = NULL,
  removeDuplicates = TRUE,
  keep = NULL,
  impute = TRUE,
  imputeType = c("random", "fixed", "beagle"),
  fixedValue = NULL,
  naStrings = NA,
  verbose = FALSE
)
```

## Arguments

- gData:

  An object of class `gData` containing at least `markers`.

- refAll:

  A character string indicating the reference allele used when recoding
  markers.  
  If "minor", then the recoding is done using the minor allele as
  reference allele. Alternatively a single character can be supplied as
  a reference allele for the whole set of SNPs, or a character vector
  with a reference allele per SNP.

- nMissGeno:

  A numerical value between 0 and 1. Genotypes with a fraction of
  missing values higher than `nMissGeno` will be removed. Genotypes with
  only missing values will always be removed.

- nMiss:

  A numerical value between 0 and 1. SNPs with a fraction of missing
  values higher than `nMiss` will be removed. SNPs with only missing
  values will always be removed.

- MAF:

  A numerical value between 0 and 1. SNPs with a Minor Allele Frequency
  (MAF) below this value will be removed. Only one of `MAF` and `MAC`
  may be specified.

- MAC:

  A numerical value. SNPs with Minor Allele Count (MAC) below this value
  will be removed. Only one of `MAF` and `MAC` may be specified.

- removeDuplicates:

  Should duplicate SNPs be removed?

- keep:

  A vector of SNPs that should never be removed in the whole process.

- impute:

  Should imputation of missing values be done?

- imputeType:

  A character string indicating what kind of imputation of values should
  be done.  

  - fixed - missing values will be replaced by a given fixed value.

  - random - missing values will be replaced by a random value
    calculated using allele frequencies per SNP.

  - beagle - missing values will be imputed using beagle software,
    version 5.2. Beagle only accepts integers as map positions. If you
    use this option, please cite the original papers in your publication
    (see references).

- fixedValue:

  A numerical value used for replacing missing values in case
  `inputType` is fixed.

- naStrings:

  A character vector of strings to be treated as NA.

- verbose:

  Should a summary of the performed steps be printed?

## Value

A copy of the input `gData` object with markers replaced by coded and
imputed markers.

## References

S R Browning and B L Browning (2007) Rapid and accurate haplotype
phasing and missing data inference for whole genome association studies
by use of localized haplotype clustering. Am J Hum Genet 81:1084-1097.
[doi:10.1086/521987](https://doi.org/10.1086/521987)

## Examples

``` r
## Create markers
markers <- matrix(c(
"AA",   "AB",   "AA",   "BB",   "BA",   "AB",   "AA",   "AA",   NA,  "AA",
"AA",   "AA",   "BB",   "BB",   "AA",   "AA",   "BB",   "AA",   NA,  "AA",
"AA",   "BA",   "AB",   "BB",   "AB",   "AB",   "AA",   "BB",   NA,  "AA",
"AA",   "AA",   "BB",   "BB",   "AA",   "AA",   "AA",   "AA",   NA,  "AA",
"AA",   "AA",   "BB",   "BB",   "AA",   "BB",   "BB",   "BB",  "AB", "AA",
"AA",   "AA",   "BB",   "BB",   "AA",    NA,    "BB",   "AA",   NA,  "AA",
"AB",   "AB",   "BB",   "BB",   "BB",   "AA",   "BB",   "BB",   NA,  "AB",
"AA",   "AA",    NA,    "BB",    NA,    "AA",   "AA",   "AA",  "AA", "AA",
"AA",    NA,     NA,    "BB",   "BB",   "BB",   "BB",   "BB",  "AA", "AA",
"AA",    NA,    "AA",   "BB",   "BB",   "BB",   "AA",   "AA",   NA,  "AA"),
ncol = 10, byrow = TRUE, dimnames = list(paste0("IND", 1:10),
paste0("SNP", 1:10)))

## create object of class 'gData'.
gData <- createGData(geno = markers)

## Code markers by minor allele, no imputation.
gDataCoded1 <- codeMarkers(gData = gData, impute = FALSE)

## Code markers by reference alleles, impute missings by fixed value.
gDataCoded2 <- codeMarkers(gData = gData,
                           refAll = rep(x = c("A", "B"), times =  5),
                           impute = TRUE, imputeType = "fixed",
                           fixedValue = 1)

## Code markers by minor allele, impute by random value.
gDataCoded3 <- codeMarkers(gData = gData, impute = TRUE,
                           imputeType = "random")
```
