# Functions for calculating kinship matrices

A collection of functions for calculating kinship matrices using
different algorithms. The following algorithms are included: astle
(Astle and Balding, 2009), Identity By State (IBS) and VanRaden
(VanRaden, 2008) for marker matrices. For method identity an identity
kinship matrix is returned.

## Usage

``` r
kinship(
  X,
  method = c("astle", "IBS", "vanRaden", "identity"),
  MAF = NULL,
  denominator = NULL
)
```

## Arguments

- X:

  An n x m marker matrix with genotypes in the rows (n) and markers in
  the columns (m).

- method:

  The method used for computing the kinship matrix.

- MAF:

  The minor allele frequency (MAF) threshold used in kinship
  computation. A numerical value between 0 and 1. SNPs with MAF below
  this value are not taken into account when computing the kinship. If
  NULL all markers are used regardless of their allele frequency.

- denominator:

  A numerical value. See details.

## Value

An n x n kinship matrix.

## Marker matrices

In all algorithms the input matrix `X` is first cleaned, i.e. markers
with a variance of 0 are excluded from the calculation of the kinship
matrix. Then some form of scaling is done which differs per algorithm.
This gives a scaled matrix `Z`. The matrix \\ZZ^t / denominator\\ is
returned. By default the denominator is equal to the number of columns
in `Z` for `astle` and `IBS` and \\2 \* p \* (1-p)\\ where \\p =
colSums(X) / (2 \* nrow(X))\\ for `vanRaden`. This denominator can be
overwritten by the user, e.g. when computing kinship matrices by
splitting `X` in smaller matrices and then adding the results together
in the end.

## References

Astle, William, and David J. Balding. 2009. “Population Structure and
Cryptic Relatedness in Genetic Association Studies.” Statistical Science
24 (4): 451–71.
[doi:10.1214/09-sts307](https://doi.org/10.1214/09-sts307) .

VanRaden P.M. (2008) Efficient methods to compute genomic predictions.
Journal of Dairy Science 91 (11): 4414–23.
[doi:10.3168/jds.2007-0980](https://doi.org/10.3168/jds.2007-0980) .

## Examples

``` r
## Create example matrix.
M <- matrix(c(1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1), nrow = 4)

## Compute kinship matrices using different methods.
kinship(M, method = "astle")
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  0.3111111  0.1333333 -0.1333333 -0.3111111
#> [2,]  0.1333333  0.6666667 -0.3111111 -0.4888889
#> [3,] -0.1333333 -0.3111111  0.3111111  0.1333333
#> [4,] -0.3111111 -0.4888889  0.1333333  0.6666667
kinship(M, method = "IBS")
#>           [,1]      [,2]      [,3]      [,4]
#> [1,] 1.0000000 0.6666667 0.6666667 0.3333333
#> [2,] 0.6666667 1.0000000 0.3333333 0.0000000
#> [3,] 0.6666667 0.3333333 1.0000000 0.6666667
#> [4,] 0.3333333 0.0000000 0.6666667 1.0000000
kinship(M, method = "vanRaden")
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  0.2857143  0.0952381 -0.0952381 -0.2857143
#> [2,]  0.0952381  0.6666667 -0.2857143 -0.4761905
#> [3,] -0.0952381 -0.2857143  0.2857143  0.0952381
#> [4,] -0.2857143 -0.4761905  0.0952381  0.6666667

## Only use markers with a Minor Allele Frequency of 0.3 or more.
kinship(M, method = "astle", MAF = 0.3)
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  0.6666667  0.6666667 -0.6666667 -0.6666667
#> [2,]  0.6666667  0.6666667 -0.6666667 -0.6666667
#> [3,] -0.6666667 -0.6666667  0.6666667  0.6666667
#> [4,] -0.6666667 -0.6666667  0.6666667  0.6666667

## Compute kinship matrix using astle and balding method with denominator 2.
kinship(M, method = "astle", denominator = 2)
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  0.4666667  0.2000000 -0.2000000 -0.4666667
#> [2,]  0.2000000  1.0000000 -0.4666667 -0.7333333
#> [3,] -0.2000000 -0.4666667  0.4666667  0.2000000
#> [4,] -0.4666667 -0.7333333  0.2000000  1.0000000
```
