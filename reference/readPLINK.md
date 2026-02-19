# Read PLINK binary data

Read PLINK binary data and save in gData format. This is a wrapper
around
[snpStats::read.plink](https://rdrr.io/pkg/snpStats/man/read.plink.html)
in the Bioconductor package `snpStats`. This package needs to be
installed for the function to work.

## Usage

``` r
readPLINK(bed, bim, fam, ...)
```

## Arguments

- bed:

  The name of the file containing the packed binary SNP genotype data.
  It should have the extension .bed; If it doesn't, then this extension
  will be appended.

- bim:

  The file containing the SNP descriptions. If not specified `bed` is
  used with its file extension replaced by bim.

- fam:

  The file containing subject (and, possibly, family) identifiers. This
  is basically a tab-delimited "pedfile". If not specified `bed` is used
  with its file extension replaced by fam.

- ...:

  Further arguments passed to
  [snpStats::read.plink](https://rdrr.io/pkg/snpStats/man/read.plink.html).

## Value

An object of class `gData`.
