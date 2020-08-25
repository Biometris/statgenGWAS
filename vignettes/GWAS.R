## ----setup, include = FALSE-----------------------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>",
fig.dim = c(7, 4)
)
library(statgenGWAS)
options(width = 100, digits = 2)

## ----loadData-------------------------------------------------------------------------------------
data(dropsMarkers)
data(dropsMap)
data(dropsPheno)

## ----convertMarkers-------------------------------------------------------------------------------
## Add genotypes as row names of dropsMarkers and drop Ind column.
rownames(dropsMarkers) <- dropsMarkers[["Ind"]]
dropsMarkers <- dropsMarkers[colnames(dropsMarkers) != "Ind"]

## ----convertMap-----------------------------------------------------------------------------------
## Add genotypes as row names of dropsMap.
rownames(dropsMap) <- dropsMap[["SNP.names"]]
## Rename Chomosome and Position columns.
colnames(dropsMap)[match(c("Chromosome", "Position"), colnames(dropsMap))] <- c("chr", "pos")

## ----createGdata----------------------------------------------------------------------------------
## Create a gData object containing map and marker information.
gDataDrops <- createGData(geno = dropsMarkers, map = dropsMap)

## ----addPheno-------------------------------------------------------------------------------------
## Rename Variety_ID to genotype.
colnames(dropsPheno)[colnames(dropsPheno) == "Variety_ID"] <- "genotype"
## Select relevant columns and convert data to a list.
dropsPhenoList <- split(x = dropsPheno[c("genotype", "grain.yield",
                                         "grain.number", "seed.size",
                                         "anthesis", "silking", "plant.height",
                                         "tassel.height", "ear.height")], 
                        f = dropsPheno[["Experiment"]])
## Add phenotypic data to gDataDrops.
gDataDrops <- createGData(gData = gDataDrops, pheno = dropsPhenoList)

## ----sumGData-------------------------------------------------------------------------------------
## Summarize gDataDrops.
summary(gDataDrops, trials = "Mur13W")

## ----removeDupMarkers-----------------------------------------------------------------------------
## Remove duplicate SNPs from gDataDrops.
gDataDropsDedup <- codeMarkers(gDataDrops, impute = FALSE, verbose = TRUE) 

## ----addMissings----------------------------------------------------------------------------------
## Copy gData object.
gDataDropsMiss <- gDataDrops
## Add random missing values to 1% of the values in the marker matrix.
set.seed(1)
nVal <- nrow(gDataDropsMiss$markers) * ncol(gDataDropsMiss$markers)
gDataDropsMiss$markers[sample(x = 1:nVal, size = nVal / 100)] <- NA

## ----imputeMissings-------------------------------------------------------------------------------
## Impute missing values with random value.
## Remove SNPs and genotypes with proportion of NA larger than 0.01.
gDataDropsImputed <- codeMarkers(gData = gDataDropsMiss, 
                                 nMissGeno = 0.01, 
                                 nMiss = 0.01, 
                                 impute = TRUE, 
                                 imputeType = "random", 
                                 verbose = TRUE)

## ----imputeMissingsBeagle, eval=FALSE-------------------------------------------------------------
#  ## Impute missing values using beagle software.
#  gDataDropsImputedBeagle <- codeMarkers(gData = gDataDropsMiss,
#                                         impute = TRUE,
#                                         imputeType = "beagle",
#                                         verbose = TRUE)

## ----stg------------------------------------------------------------------------------------------
## Run single trait GWAS for traits 'grain.yield' and 'anthesis' for trial Mur13W.
GWASDrops <- runSingleTraitGwas(gData = gDataDropsDedup,
                                trials = "Mur13W",
                                traits = c("grain.yield", "anthesis"))

## ----gwaRes---------------------------------------------------------------------------------------
print(head(GWASDrops$GWAResult$Mur13W), row.names = FALSE)

## ----signSnp--------------------------------------------------------------------------------------
print(GWASDrops$signSnp$Mur13W, row.names = FALSE)

## ----sumStg---------------------------------------------------------------------------------------
## Create summary of GWASDrops.
summary(GWASDrops)

## ----qqStg----------------------------------------------------------------------------------------
## Plot a QQ-plot of GWAS Drops.
plot(GWASDrops, plotType = "qq", trait = "grain.yield")

## ----manhattanStg---------------------------------------------------------------------------------
## Plot a manhattan plot of GWAS Drops.
plot(GWASDrops, plotType = "manhattan", trait = "grain.yield")

## ----manhattanStgThr------------------------------------------------------------------------------
## Plot a manhattan plot of GWAS Drops.
## Set significance threshold to 4 and only plot chromosomes 6 to 8.
plot(GWASDrops, plotType = "manhattan", trait = "grain.yield", yThr = 4, chr = 6:8)

## ----manhattanLod---------------------------------------------------------------------------------
## Plot a manhattan plot of GWAS Drops.
## Plot only 5% of SNPs with a LOD below 3.
set.seed(1)
plot(GWASDrops, plotType = "manhattan", trait = "grain.yield", lod = 3)

## ----qtlStg---------------------------------------------------------------------------------------
## Plot a qtl plot of GWAS Drops for Mur13W.
plot(GWASDrops, plotType = "qtl")

## ----qtlStgThr------------------------------------------------------------------------------------
## Plot a qtl plot of GWAS Drops for Mur13W.
## Set significance threshold to 4.
plot(GWASDrops, plotType = "qtl", yThr = 4)

## ----qtlStgNorm-----------------------------------------------------------------------------------
## Plot a qtl plot of GWAS Drops for Mur13W.
## Set significance threshold to 4 and normalize effect estimates.
plot(GWASDrops, plotType = "qtl", yThr = 4, normalize = TRUE)

## ----stgChrSpec-----------------------------------------------------------------------------------
## Run single trait GWAS for trial 'Mur13W' and trait 'grain.yield'
## Use chromosome specific kinship matrices computed using method of van Raden.
GWASDropsChrSpec <- runSingleTraitGwas(gData = gDataDropsDedup, 
                                       traits = "grain.yield",
                                       trials = "Mur13W",
                                       GLSMethod = "multi",
                                       kinshipMethod = "vanRaden")

## ----stgSNPFixThr---------------------------------------------------------------------------------
## Run single trait GWAS for trait 'grain.yield' for Mur13W.
## Use a fixed significance threshold of 4.
GWASDropsFixThr <- runSingleTraitGwas(gData = gDataDropsDedup,
                                      trials = "Mur13W",
                                      traits = "grain.yield",
                                      thrType = "fixed",
                                      LODThr = 4)

## ----stgSNPNR-------------------------------------------------------------------------------------
## Run single trait GWAS for trait 'grain.yield' for Mur13W.
## Use the Newton Raphson algorithm for computing the variance components.
GWASDropsNR <- runSingleTraitGwas(gData = gDataDropsDedup,
                                  trials = "Mur13W",
                                  traits = "grain.yield",
                                  remlAlgo = "NR")

## ----inflation------------------------------------------------------------------------------------
GWASDrops$GWASInfo$inflationFactor$Mur13W

## ----stgSNPGenomicCorrection----------------------------------------------------------------------
## Run single trait GWAS for trait 'grain.yield' for Mur13W.
## Perform genomic correction on the p-Values.
GWASDropsGenControl <- runSingleTraitGwas(gData = gDataDropsDedup,
                                          trials = "Mur13W",
                                          traits = "grain.yield",
                                          genomicControl = TRUE)

## ----stgSNPCovar----------------------------------------------------------------------------------
## Run single trait GWAS for trait 'grain.yield' for Mur13W.
## Use PZE-106021410, the most significant SNP, a SNP covariate.
GWASDropsSnpCov <- runSingleTraitGwas(gData = gDataDropsDedup,
                                      trials = "Mur13W",
                                      traits = "grain.yield",
                                      snpCov = "PZE-106021410")

## ----stgMAC---------------------------------------------------------------------------------------
## Run single trait GWAS for trait 'grain.yield' for Mur13W.
## Only include SNPs that have a MAC of at least 20
GWASDropsMAC <- runSingleTraitGwas(gData = gDataDropsDedup,
                                   trials = "Mur13W",
                                   traits = "grain.yield",
                                   useMAF = FALSE,
                                   MAC = 20)

## ----stgInclReg-----------------------------------------------------------------------------------
## Run single trait GWAS for trait 'grain.yield' for Mur13W.
## Include SNPs within 200000 centimorgan of significant SNPs with a minimum LD of 0.1.
GWASDropsInclClose <- runSingleTraitGwas(gData = gDataDropsDedup,
                                         trials = "Mur13W",
                                         traits = "grain.yield",
                                         sizeInclRegion = 200000,
                                         minR2 = 0.1)
## Check signSnp in output.
print(head(GWASDropsInclClose$signSnp$Mur13W), row.names = FALSE)

