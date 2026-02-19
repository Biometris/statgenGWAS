# Summary function for the class `GWAS`

Gives a summary for an object of S3 class `GWAS`.

## Usage

``` r
# S3 method for class 'GWAS'
summary(object, ..., trials = NULL, traits = NULL)
```

## Arguments

- object:

  An object of class `GWAS`.

- ...:

  Not used.

- trials:

  A vector of strings or numeric indices indicating for which trials the
  summary should be made. If `NULL`, a summary is made for all trials.

- traits:

  A vector of strings indicating which traits to include in the summary.
  If `NULL`, all traits are included in the summary.
