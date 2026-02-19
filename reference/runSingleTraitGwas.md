# Perform single-trait GWAS

`runSingleTraitGwas` performs a single-trait Genome Wide Association
Study (GWAS) on phenotypic and genotypic data contained in a `gData`
object. A covariance matrix is computed using the EMMA algorithm (Kang
et al., 2008), the Newton-Raphson algorithm (Tunnicliffe, 1989) in the
`sommer` package, or the modified Henderson algorithm in the `LMMsolver`
package. Then a Generalized Least Squares (GLS) method is used for
estimating the marker effects and corresponding p-values. This is done
using either one kinship matrix for all chromosomes or different
chromosome-specific kinship matrices for each chromosome. Significant
SNPs are selected based on a user defined threshold.

## Usage

``` r
runSingleTraitGwas(
  gData,
  traits = NULL,
  trials = NULL,
  covar = NULL,
  snpCov = NULL,
  kin = NULL,
  kinshipMethod = c("astle", "IBS", "vanRaden"),
  remlAlgo = c("EMMA", "NR", "Hend"),
  GLSMethod = c("single", "multi"),
  useMAF = TRUE,
  MAF = 0.01,
  MAC = 10,
  genomicControl = FALSE,
  thrType = c("bonf", "fixed", "small", "fdr"),
  alpha = 0.05,
  LODThr = 4,
  nSnpLOD = 10,
  pThr = 0.05,
  rho = 0.5,
  sizeInclRegion = 0,
  minR2 = 0.5,
  nCores = NULL
)
```

## Arguments

- gData:

  An object of class `gData` containing at least `map`, `markers` and
  `pheno`.

- traits:

  A vector of traits on which to run GWAS. These can be either numeric
  indices or character names of columns in `pheno`. If `NULL`, GWAS is
  run on all traits.

- trials:

  A vector of trials on which to run GWAS. These can be either numeric
  indices or character names of list items in `pheno`. If `NULL`, GWAS
  is run for all trials. GWAS is run for the selected trials in
  sequential order.

- covar:

  An optional vector of covariates taken into account when running GWAS.
  These can be either numeric indices or character names of columns in
  `covar` in `gData`. If `NULL` no covariates are used.

- snpCov:

  An optional character vector of snps to be included as covariates.

- kin:

  An optional kinship matrix or list of kinship matrices. These matrices
  can be from the `matrix` class as defined in the base package or from
  the `dsyMatrix` class, the class of symmetric matrices in the Matrix
  package.  
  If `GLSMethod` = "single" then one matrix should be provided, if
  `GLSMethod` = "multi", a list of chromosome specific matrices of
  length equal to the number of chromosomes in `map` in `gData`.  
  If `NULL` then matrix `kinship` in `gData` is used.  
  If both `kin` is provided and `gData` contains a matrix `kinship` then
  `kin` is used.

- kinshipMethod:

  An optional character indicating the method used for calculating the
  kinship matrix(ces). Currently "astle" (Astle and Balding, 2009),
  "IBS" and "vanRaden" (VanRaden, 2008) are supported. If a kinship
  matrix is supplied either in `gData` or in parameter `kin`,
  `kinshipMethod` is ignored.

- remlAlgo:

  A character string indicating the algorithm used to estimate the
  variance components. Either `EMMA`, for the EMMA algorithm, `NR`, for
  the Newton-Raphson algorithm (using
  [sommer::mmes](https://rdrr.io/pkg/sommer/man/mmes.html)), or `Hend`
  for the modified Henderson algorithm (using
  [LMMsolver::LMMsolve](https://rdrr.io/pkg/LMMsolver/man/LMMsolve.html)).

- GLSMethod:

  A character string indicating the method used to estimate the marker
  effects. Either `single` for using a single kinship matrix, or `multi`
  for using chromosome specific kinship matrices.

- useMAF:

  Should the minor allele frequency be used for selecting SNPs for the
  analysis. If `FALSE`, the minor allele count is used instead.

- MAF:

  The minor allele frequency (MAF) threshold used in GWAS. A numerical
  value between 0 and 1. SNPs with MAF below this value are not taken
  into account in the analysis, i.e. p-values and effect sizes are put
  to missing (`NA`). Ignored if `useMAF` is `FALSE`.

- MAC:

  A numerical value. SNPs with minor allele count below this value are
  not taken into account for the analysis, i.e. p-values and effect
  sizes are set to missing (`NA`). Ignored if `useMAF` is `TRUE`.

- genomicControl:

  Should genomic control correction as in Devlin and Roeder (1999) be
  applied?

- thrType:

  A character string indicating the type of threshold used for the
  selection of candidate loci. See Details.

- alpha:

  A numerical value used for calculating the LOD-threshold for `thrType`
  = "bonf" and the significant p-Values for `thrType` = "fdr".

- LODThr:

  A numerical value used as a LOD-threshold when `thrType` = "fixed".

- nSnpLOD:

  A numerical value indicating the number of SNPs with the smallest
  p-values that are selected when `thrType` = "small".

- pThr:

  A numerical value just as the cut off value for p-Values for `thrType`
  = "fdr".

- rho:

  A numerical value used a the minimum value for SNPs to be considered
  correlated when using `thrType` = "fdr".

- sizeInclRegion:

  An integer. Should the results for SNPs close to significant SNPs be
  included? If so, the size of the region in centimorgan or base pairs.
  Otherwise 0.

- minR2:

  A numerical value between 0 and 1. Restricts the SNPs included in the
  region close to significant SNPs to only those SNPs that are in
  sufficient Linkage Disequilibrium (LD) with the significant snp, where
  LD is measured in terms of \\R^2\\. If for example `sizeInclRegion` =
  200000 and `minR2` = 0.5, then for every significant SNP also those
  SNPs whose LD (\\R^2\\) with the significant SNP is at least 0.5 AND
  which are at most 200000 away from this significant snp are included.
  Ignored if `sizeInclRegion` = 0.

- nCores:

  A numerical value indicating the number of cores to be used by the
  parallel part of the algorithm. If `NULL` the number of cores used
  will be equal to the number of cores available on the machine - 1.

## Value

An object of class `GWAS`.

## Threshold types for the selection of candidate loci

For the selection of candidate loci from the GWAS output four different
methods can be used. The method used can be specified in the function
parameter `thrType`. Further parameters can be used to fine tune the
method.

- bonf:

  The Bonferroni threshold, a LOD-threshold of \\-log10(alpha / m)\\,
  where m is the number of SNPs and alpha can be specified by the
  function parameter `alpha`

- fixed:

  A fixed LOD-threshold, specified by the function parameter `LODThr`

- small:

  The n SNPs with with the smallest p-Values are selected. n can be
  specified in `nSnpLOD`

- fdr:

  Following the algorithm proposed by Brzyski D. et al. (2017), SNPs are
  selected in such a way that the False Discovery Rate (fdr) is
  minimized. To do this, first the SNPs are restricted to the SNPs with
  a p-Value below `pThr`. Then clusters of SNPs are created using a two
  step iterative process in which first the SNP with the lowest p-Value
  is selected as cluster representative. This SNP and all SNPs that have
  a correlation with this SNP of \\\rho\\ or higher will form a cluster.
  \\\rho\\ can be specified as an argument in the function and has a
  default value of 0.5, which is a recommended starting value in
  practice. The selected SNPs are removed from the list and the
  procedure repeated until no SNPs are left. Finally to determine the
  number of significant clusters the first cluster is determined for
  which the p-Value of the cluster representative is not smaller than
  \\cluster\_{number} \* \alpha / m\\ where m is the number of SNPs and
  alpha can be specified by the function parameter `alpha`. All previous
  clusters are selected as significant.

## References

Astle, William, and David J. Balding. 2009. Population Structure and
Cryptic Relatedness in Genetic Association Studies. Statistical Science
24 (4): 451–71.
[doi:10.1214/09-sts307](https://doi.org/10.1214/09-sts307) .

Brzyski D. et al. (2017) Controlling the Rate of GWAS False Discoveries.
Genetics 205 (1): 61-75.
[doi:10.1534/genetics.116.193987](https://doi.org/10.1534/genetics.116.193987)

Devlin, B., and Kathryn Roeder. 1999. Genomic Control for Association
Studies. Biometrics 55 (4): 997–1004.
[doi:10.1111/j.0006-341x.1999.00997.x](https://doi.org/10.1111/j.0006-341x.1999.00997.x)
.

Kang et al. (2008) Efficient Control of Population Structure in Model
Organism Association Mapping. Genetics 178 (3): 1709–23.
[doi:10.1534/genetics.107.080101](https://doi.org/10.1534/genetics.107.080101)
.

Millet, E. J., Pommier, C., et al. (2019). A multi-site experiment in a
network of European fields for assessing the maize yield response to
environmental scenarios - Data set.
[doi:10.15454/IASSTN](https://doi.org/10.15454/IASSTN)

Rincent et al. (2014) Recovering power in association mapping panels
with variable levels of linkage disequilibrium. Genetics 197 (1):
375–87.
[doi:10.1534/genetics.113.159731](https://doi.org/10.1534/genetics.113.159731)
.

Segura et al. (2012) An efficient multi-locus mixed-model approach for
genome-wide association studies in structured populations. Nature
Genetics 44 (7): 825–30.
[doi:10.1038/ng.2314](https://doi.org/10.1038/ng.2314) .

Sun et al. (2010) Variation explained in mixed-model association
mapping. Heredity 105 (4): 333–40.
[doi:10.1038/hdy.2010.11](https://doi.org/10.1038/hdy.2010.11) .

Tunnicliffe W. (1989) On the use of marginal likelihood in time series
model estimation. JRSS 51 (1): 15–27.

VanRaden P.M. (2008) Efficient methods to compute genomic predictions.
Journal of Dairy Science 91 (11): 4414–23.
[doi:10.3168/jds.2007-0980](https://doi.org/10.3168/jds.2007-0980) .

## See also

[summary.GWAS](summary.GWAS.md), [plot.GWAS](plot.GWAS.md)

## Examples

``` r
## Create a gData object Using the data from the DROPS project.
## See the included vignette for a more extensive description on the steps.
data(dropsMarkers)
data(dropsMap)
data(dropsPheno)

## Add genotypes as row names of dropsMarkers and drop Ind column.
rownames(dropsMarkers) <- dropsMarkers[["Ind"]]
dropsMarkers <- dropsMarkers[colnames(dropsMarkers) != "Ind"]

## Add genotypes as row names of dropsMap.
rownames(dropsMap) <- dropsMap[["SNP.names"]]

## Rename Chomosome and Position columns.
colnames(dropsMap)[match(c("Chromosome", "Position"), 
                   colnames(dropsMap))] <- c("chr", "pos")
                   
## Convert phenotypic data to a list.
dropsPhenoList <- split(x = dropsPheno, f = dropsPheno[["Experiment"]])

## Rename Variety_ID to genotype and select relevant columns.
dropsPhenoList <- lapply(X = dropsPhenoList, FUN = function(trial) {
  colnames(trial)[colnames(trial) == "Variety_ID"] <- "genotype"
  trial <- trial[c("genotype", "grain.yield", "grain.number", "seed.size",
                 "anthesis", "silking", "plant.height", "tassel.height",
                 "ear.height")]
return(trial)
}) 

## Create gData object.
gDataDrops <- createGData(geno = dropsMarkers, map = dropsMap, 
                          pheno = dropsPhenoList)
                          
## Run single trait GWAS for trait 'grain.yield' for trial Mur13W.
# \donttest{
GWASDrops <- runSingleTraitGwas(gData = gDataDrops,
                               trials = "Mur13W",
                               traits = "grain.yield")
# }
                               
## Run single trait GWAS for trait 'grain.yield' for trial Mur13W.
## Use chromosome specific kinship matrices calculated using vanRaden method.
# \donttest{
GWASDropsMult <- runSingleTraitGwas(gData = gDataDrops,
                                    trials = "Mur13W",
                                    traits = "grain.yield",
                                    kinshipMethod = "vanRaden",
                                    GLSMethod = "multi")  
# }
```
