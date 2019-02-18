## ----eval=FALSE----------------------------------------------------------
#  markersRaw <- read.csv(gzcon(
#    url("http://ricediversity.org/data/sets/44kgwas/RiceDiversity.44K.MSU6.Genotypes.csv.gz"),
#    text = TRUE), stringsAsFactors = FALSE)
#  phenoRaw <- read.delim(
#    "http://ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt",
#    stringsAsFactors = FALSE)
#  covarRaw <- read.csv(
#    "http://ricediversity.org/data/sets/44kgwas/RiceDiversity.44K.germplasm.csv",
#    skip = 1,
#    stringsAsFactors = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  ## Create mapRice with columns chr and pos and SNP names as row names.
#  mapRice <- data.frame(chr = markersRaw$chr,
#                        pos = markersRaw$position,
#                        row.names = markersRaw$id)
#  ## Create markersRice by removing the 'map' columns from the raw data and then transposing the result
#  markersRice <- t(markersRaw[, 4:ncol(markersRaw)])
#  ## Add SNP names as column names
#  colnames(markersRice) <- markersRaw$id
#  ## Missing values are coded in more than one way. Set everything that is not A, C, T or G to NA.
#  markersRice[!markersRice %in% c("A", "C", "T", "G")] <- NA

## ----eval=FALSE----------------------------------------------------------
#  library(genStatPipeline)
#  gDataRice <- createGData(geno = markersRice, map = mapRice)

## ----eval=FALSE----------------------------------------------------------
#  phenoRice <- cbind(genotype = paste0("NSFTV_", phenoRaw$NSFTVID),
#                     phenoRaw[, 3:(ncol(phenoRaw) - 2)],
#                     stringsAsFactors = FALSE)
#  covarRice <- covarRaw[!is.na(covarRaw$NSFTV.ID), c("PC1", "PC2", "PC3", "PC4")]
#  rownames(covarRice) <- paste0("NSFTV_", covarRaw[!is.na(covarRaw$NSFTV.ID), "NSFTV.ID"])

## ----eval=FALSE----------------------------------------------------------
#  gDataRice <- createGData(gData = gDataRice, pheno = phenoRice, covar = covarRice)

## ----eval=FALSE----------------------------------------------------------
#  summary(gDataRice)

## ----eval=FALSE----------------------------------------------------------
#  gDataRiceCoded <- codeMarkers(gData = gDataRice)

## ----eval=FALSE----------------------------------------------------------
#  ## Perform imputation using beagle.
#  gDataRiceCodedBeagle <- codeMarkers(gdata = gDataRice,
#                                      removeDuplicates = FALSE,
#                                      imputeType = "beagle")

## ----eval=FALSE----------------------------------------------------------
#  ## Load the data imputed by beagle directly.
#  data(gDataRiceCodedBeagle)

## ----eval=FALSE----------------------------------------------------------
#  summary(gDataRiceCodedBeagle)

## ----eval=FALSE----------------------------------------------------------
#  GWASRicePlantHeight0 <- runSingleTraitGwas(gData = gDataRiceCodedBeagle,
#                                             traits = "Plant.height")

## ----eval=FALSE----------------------------------------------------------
#  summary(GWASRicePlantHeight0)
#  plot(GWASRicePlantHeight0, type = "qq")
#  plot(GWASRicePlantHeight0)

## ----eval=FALSE----------------------------------------------------------
#  GWASRicePlantHeight <-
#    runSingleTraitGwas(gData = gDataRiceCodedBeagle, traits = "Plant.height",
#                       covar = c("PC1", "PC2", "PC3", "PC4"), GLSMethod = "multi",
#                       kinshipMethod = "IBS", thrType = "fixed", LODThr = 4)
#  
#  summary(GWASRicePlantHeight)
#  plot(GWASRicePlantHeight, type = "qq")
#  plot(GWASRicePlantHeight)

## ----eval=FALSE----------------------------------------------------------
#  GWASRicePlantHeightMult <-
#    runSingleTraitGwas(gData = gDataRiceCodedBeagle,
#                       covar = c("PC1", "PC2", "PC3", "PC4"), GLSMethod = "multi",
#                       kinshipMethod = "IBS", thrType = "fixed", LODThr = 4)
#  
#  summary(GWASRicePlantHeightMult)
#  plot(GWASRicePlantHeightMult, type = "qq", trait = "Plant.height")
#  plot(GWASRicePlantHeightMult, type = "qtl")

