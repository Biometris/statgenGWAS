# Summary function for the class `gData`

Gives a summary for an object of S3 class `gData`.

## Usage

``` r
# S3 method for class 'gData'
summary(object, ..., trials = NULL)
```

## Arguments

- object:

  An object of class `gData`.

- ...:

  Not used.

- trials:

  A vector of trials to include in the summary. These can be either
  numeric indices or character names of list items in `pheno`. If
  `NULL`, all trials are included.

## Value

A list with a most four components:

- mapSum:

  A list with number of markers and number of chromosomes in the map.

- markerSum:

  A list with number of markers, number of genotypes and the
  distribution of the values within the markers.

- phenoSum:

  A list of data.frames, one per trial with a summary of all traits
  within the trial.

- covarSum:

  A list of data.frames, one per trial with a summary of all covariates
  within the trial.

All components are only present in the output if the corresponding
content is present in the gData object.
