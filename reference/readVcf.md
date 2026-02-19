# Read variant call format data

Read variant call format (VCF) data and save in gData format. This is a
wrapper around vcfR::read.vcfR in the package `vcfR`. This package needs
to be installed for the function to work.

## Usage

``` r
readVcf(vcfFile, ...)
```

## Arguments

- vcfFile:

  The name of the vcf file. This can either be a plain text file with
  extension (.vcf), or a gzipped file with extencion (.vcf.gz).

- ...:

  Further arguments passed to vcfR::read.vcfR.

## Value

An object of class `gData`.

## References

Knaus BJ, Grünwald NJ (2017). “VCFR: a package to manipulate and
visualizevariant call format data in R.” *Molecular Ecology Resources*,
*17*(1), 44-53. ISSN 757,
[doi:10.1111/1755-0998.12549](https://doi.org/10.1111/1755-0998.12549) .
