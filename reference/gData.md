# S3 Class gData

`createGData` creates an object of S3 class gData with genotypic and
phenotypic data for usage in further analysis. All input to the function
is optional, however at least one input should be provided. It is
possible to provide an existing `gData` object as additional input in
which case data is added to this object. Existing data will be
overwritten with a warning.

## Usage

``` r
createGData(
  gData = NULL,
  geno = NULL,
  map = NULL,
  kin = NULL,
  pheno = NULL,
  covar = NULL
)
```

## Arguments

- gData:

  An optional gData object to be modified. If `NULL`, a new gData object
  is created.

- geno:

  A matrix or data.frame with genotypes in the rows and markers in the
  columns. A matrix from the `matrix` in the base package may be
  provided as well as as matrix from the Matrix package.  
  If no row names are provided, they are taken from `pheno` (if supplied
  and dimension matches). If no column names are provided, the row names
  from `map` are used (if supplied and dimension matches).

- map:

  A data.frame with columns `chr` for chromosome and `pos` for position.
  Positions can be in base pair (bp) or centimorgan (cM). They should
  not be cumulative over the chromosomes. Other columns are ignored.
  Marker names should be in the row names. These should match the marker
  names in `geno` (if supplied).

- kin:

  A kinship matrix or list of kinship matrices with genotype in rows and
  colums. These matrices can be from the `matrix` class, as defined in
  the base package, or from the `dsyMatrix` class, the class of
  symmetric matrices in the Matrix package.  
  The genotypes should be identical to the genotypes in `geno`.  
  If a list of kinship matrices is provided these are supposed to be
  chromosome specific matrices. In that case their names should match
  the names of the chromosomes in `map`. If no names are provided, the
  number of matrices should match the number of chromosomes in `map` in
  which case default names are provided.

- pheno:

  A data.frame or a list of data.frames with phenotypic data, with
  genotypes in the first column `genotype` and traits in the following
  columns. The trait columns should be numerical columns only. A list of
  data.frames can be used for replications, i.e. different trials.

- covar:

  A data.frame with extra covariates per genotype. Genotypes should be
  in the rows.

## Value

An object of class `gData` with the following components:

- `map`:

  a data.frame containing map data. Map is sorted by chromosome and
  position.

- `markers`:

  a matrix containing marker information.

- `pheno`:

  a list of data.frames containing phenotypic data.

- `kinship`:

  a kinship matrix.

- `covar`:

  a data.frame with extra covariates.

## See also

[summary.gData](summary.gData.md)

## Author

Bart-Jan van Rossum

## Examples

``` r
set.seed(1234)
## Create genotypic data.
geno <- matrix(sample(x = c(0, 1, 2), size = 15, replace = TRUE), nrow = 3)
dimnames(geno) <- list(paste0("G", 1:3), paste0("M", 1:5))

## Construct map.
map <- data.frame(chr = c(1, 1, 2, 2, 2), pos = 1:5,
                  row.names = paste0("M", 1:5))

## Compute kinship matrix.
kin <- kinship(X = geno, method = "IBS")

## Create phenotypic data.
pheno <- data.frame(paste0("G", 1:3),
                    matrix(rnorm(n = 12, mean = 50, sd = 5), nrow = 3),
                    stringsAsFactors = FALSE)
dimnames(pheno) = list(paste0("G", 1:3), c("genotype", paste0("T", 1:4)))

## Combine all data in gData object.
gData <- createGData(geno = geno, map = map, kin = kin, pheno = pheno)
summary(gData)
#> map
#>  Number of markers: 5 
#>  Number of chromosomes: 2 
#> 
#> markers
#>  Number of markers: 5 
#>  Number of genotypes: 3 
#>  Content:
#>      0    1    2 <NA>  
#>    0.2  0.6  0.2  0.0  
#> 
#> pheno
#>  Number of trials: 1 
#> 
#>  pheno:
#>      Number of traits: 4 
#>      Number of genotypes: 3 
#> 
#>               T1       T2       T3       T4
#> Min.    41.24633 41.56337 46.76490 51.55131
#> 1st Qu. 44.33149 44.21309 48.42824 51.71474
#> Median  47.41665 46.86282 50.09158 51.87818
#> Mean    47.68783 48.42541 50.12757 52.59013
#> 3rd Qu. 50.90859 51.85644 51.80890 53.10954
#> Max.    54.40052 56.85005 53.52622 54.34090
#> NA's     0.00000  0.00000  0.00000  0.00000
#> 

## Construct covariate.
covar <- data.frame(C1 = c("a", "a", "b"), row.names = paste0("G", 1:3))

## Compute alternative kinship matrix.
kin2 <- kinship(X = geno, method = "astle")

## Add covariates to previously created gData object and overwrite
## current kinship matrix by newly computed one.
gData2 <- createGData(gData = gData, kin = kin2, covar = covar)
#> Warning: Existing kinship will be overwritten.
```
