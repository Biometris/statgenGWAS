## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>",
fig.dim = c(7, 4)
)
library(statgenGWAS)

## ----loadData-----------------------------------------------------------------
data(dropsMarkers)
data(dropsMap)
data(dropsPheno)

## ----convertMarkers-----------------------------------------------------------
## Add genotypes as row names of dropsMarkers and drop Ind column.
rownames(dropsMarkers) <- dropsMarkers$Ind
dropsMarkers <- dropsMarkers[colnames(dropsMarkers) != "Ind"]

## ----convertMap---------------------------------------------------------------
## Add genotypes as row names of dropsMap.
rownames(dropsMap) <- dropsMap$SNP.names
## Rename Chomosome and Position columns.
colnames(dropsMap)[match(c("Chromosome", "Position"), colnames(dropsMap))] <-
  c("chr", "pos")

## ----createGdata--------------------------------------------------------------
## Create a gData object containing map and marker information.
gDataDrops <- createGData(geno = dropsMarkers, map = dropsMap)

## ----addPheno-----------------------------------------------------------------
## Convert phenotypic data to a list
dropsPhenoList <- split(x = dropsPheno, f = dropsPheno$Experiment)
## Rename Variety_ID to genotype and select relevant columns.
dropsPhenoList <- lapply(X = dropsPhenoList, FUN = function(trial) {
  colnames(trial)[colnames(trial) == "Variety_ID"] <- "genotype"
  trial <- trial[c("genotype", "grain.yield", "grain.number", "seed.size",
                   "anthesis", "silking", "plant.height", "tassel.height",
                   "ear.height")]
  return(trial)
})
## Add phenotypic data to gDataDrops
gDataDrops <- createGData(gData = gDataDrops, pheno = dropsPhenoList)

## ----sumGData-----------------------------------------------------------------
## Summarize gDataDrops
summary(gDataDrops, trials = "Bol12R")

## ----removeDupMarkers---------------------------------------------------------
## Remove duplicate markers from gDataDrops
gDataDropsDedup <- codeMarkers(gDataDrops, impute = FALSE, verbose = TRUE) 

## ----addMissings--------------------------------------------------------------
## Copy gData object.
gDataDropsMiss <- gDataDrops
## Add random missing values to 1% of the values in the marker matrix.
set.seed(1)
nVal <- nrow(gDataDropsMiss$markers) * ncol(gDataDropsMiss$markers)
gDataDropsMiss$markers[sample(x = 1:nVal, size = nVal / 100)] <- NA

## ----imputeMissings-----------------------------------------------------------
## Impute missing values with random value.
## Remove SNPs and genotypes with proportion of NA larger than 0.01
gDataDropsImputed <- codeMarkers(gData = gDataDropsMiss, nMissGeno = 0.01, 
                                 nMiss = 0.01, impute = TRUE, 
                                 imputeType = "random", verbose = TRUE)

## ----imputeMissingsBeagle, eval=FALSE-----------------------------------------
#  ## Impute missing values using beagle software.
#  gDataDropsImputedBeagle <- codeMarkers(gData = gDataDropsMiss, impute = TRUE,
#                                         imputeType = "beagle", verbose = TRUE)

## ----stg----------------------------------------------------------------------
## Run single trait GWAS for trial 'Bol12W' and trait 'grain.yield'
GWASDrops <- runSingleTraitGwas(gData = gDataDropsDedup, traits = "grain.yield",                                 trials = "Bol12W")

## ----sumStg-------------------------------------------------------------------
## Create summary of GWASDrops.
summary(GWASDrops)

## ----qqStg--------------------------------------------------------------------
## Plot a qq plot of GWAS Drops.
plot(GWASDrops, plotType = "qq")

## ----manhattanStg-------------------------------------------------------------
## Plot a manhattan plot of GWAS Drops.
plot(GWASDrops, plotType = "manhattan")

## ----eval=FALSE---------------------------------------------------------------
#  ## Run single trait GWAS for trial 'Bol12W' and trait 'grain.yield'
#  ## Use chromosome specific kinship matrices computed using method of van Raden.
#  GWASDropsChrSpec <- runSingleTraitGwas(gData = gDataDropsDedup,
#                                         traits = "grain.yield",
#                                         trials = "Bol12W",
#                                         GLSMethod = "multi",
#                                         kinshipMethod = "vanRaden")

