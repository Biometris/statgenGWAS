#' Perform single-trait GWAS
#'
#' \code{runSingleTraitGwas} performs a single-trait Genome Wide Association
#' Study (GWAS) on phenotypic and genotypic data contained in a \code{gData}
#' object. A covariance matrix is computed using the EMMA algorithm (Kang et
#' al., 2008) or the Newton-Raphson algorithm (Tunnicliffe, 1989) in the
#' \code{sommer} package. Then a Generalized Least Squares (GLS) method is used
#' for estimating the marker effects and corresponding p-values. This is done
#' using either one kinship matrix for all chromosomes or different chromosome-specific
#' kinship matrices for each chromosome. Significant SNPs are selected
#' based on a user defined threshold.
#'
#' @param gData An object of class \code{gData} containing at least \code{map},
#' \code{markers} and \code{pheno}.
#' @param traits A vector of traits on which to run GWAS. These can be either
#' numeric indices or character names of columns in \code{pheno}. If \code{NULL},
#' GWAS is run on all traits.
#' @param trials A vector of trials on which to run GWAS. These can
#' be either numeric indices or character names of list items in \code{pheno}.
#' If \code{NULL}, GWAS is run for all trials. GWAS is run for the
#' selected trials in sequential order.
#' @param covar An optional vector of covariates taken into account when
#' running GWAS. These can be either numeric indices or character names of
#' columns in \code{covar} in \code{gData}. If \code{NULL} no covariates are
#' used.
#' @param snpCov An optional character vector of snps to be included as
#' covariates.
#' @param kin An optional kinship matrix or list of kinship matrices. These
#' matrices can be from the \code{matrix} class as defined in the base package
#' or from the \code{dsyMatrix} class, the class of symmetric matrices in the
#' Matrix package.\cr
#' If \code{GLSMethod} = "single" then one matrix should be provided, if
#' \code{GLSMethod} = "multi", a list of chromosome specific matrices of length
#' equal to the number of chromosomes in \code{map} in \code{gData}.\cr
#' If \code{NULL} then matrix \code{kinship} in \code{gData} is used. \cr
#' If both \code{kin} is provided and \code{gData} contains a matrix
#' \code{kinship} then \code{kin} is used.
#' @param kinshipMethod An optional character indicating the method used for
#' calculating the kinship matrix(ces). Currently "astle" (Astle and Balding,
#' 2009), "IBS" and "vanRaden" (VanRaden, 2008) are supported. If a
#' kinship matrix is supplied either in \code{gData} or in parameter \code{kin},
#' \code{kinshipMethod} is ignored.
#' @param remlAlgo A character string indicating the algorithm used to estimate
#' the variance components. Either \code{EMMA}, for the EMMA algorithm, or
#' \code{NR}, for the Newton-Raphson algorithm.
#' @param GLSMethod A character string indicating the method used to estimate
#' the marker effects. Either \code{single} for using a single kinship matrix,
#' or \code{multi} for using chromosome specific kinship matrices.
#' @param useMAF Should the minor allele frequency be used for selecting SNPs
#' for the analysis. If \code{FALSE}, the minor allele count is used instead.
#' @param MAF The minor allele frequency (MAF) threshold used in GWAS. A
#' numerical value between 0 and 1. SNPs with MAF below this value are not taken
#' into account in the analysis, i.e. p-values and effect sizes are put to
#' missing (\code{NA}). Ignored if \code{useMAF} is \code{FALSE}.
#' @param MAC A numerical value. SNPs with minor allele count below this value
#' are not taken into account for the analysis, i.e. p-values and effect sizes
#' are set to missing (\code{NA}). Ignored if \code{useMAF} is \code{TRUE}.
#' @param genomicControl Should genomic control correction as in Devlin and
#' Roeder (1999) be applied?
#' @param thrType A character string indicating the type of threshold used for
#' the selection of candidate loci. Either \code{bonf} for using the
#' Bonferroni threshold, a LOD-threshold of \eqn{-log10(alpha/p)}, where p is
#' the number of markers and alpha can be specified in \code{alpha},
#' \code{fixed} for a self-chosen fixed LOD-threshold, specified in \code{LODThr}
#' or \code{small}, the LOD-threshold is chosen such as the SNPs with the
#' \code{nSnpLOD} smallest p-values are selected. \code{nSnpLOD} can be
#' specified.
#' @param alpha A numerical value used for calculating the LOD-threshold for
#' \code{thrType} = "bonf".
#' @param LODThr A numerical value used as a LOD-threshold when
#' \code{thrType} = "fixed".
#' @param nSnpLOD A numerical value indicating the number of SNPs with the
#' smallest p-values that are selected when \code{thrType} = "small".
#' @param rho A numerical value ...
#' @param pThr A numerical value ...
#' @param sizeInclRegion An integer. Should the results for SNPs close to
#' significant SNPs be included? If so, the size of the region in centimorgan
#' or base pairs. Otherwise 0.
#' @param minR2 A numerical value between 0 and 1. Restricts the SNPs included
#' in the region close to significant SNPs to only those SNPs that are in
#' sufficient Linkage Disequilibrium (LD) with the significant snp, where LD
#' is measured in terms of \eqn{R^2}. If for example \code{sizeInclRegion} = 
#' 200000 and \code{minR2} = 0.5, then for every significant SNP also those SNPs 
#' whose LD (\eqn{R^2}) with the significant SNP is at least 0.5 AND which are 
#' at most 200000 away from this significant snp are included. Ignored if
#' \code{sizeInclRegion} = 0.
#' @param nCores A numerical value indicating the number of cores to be used by
#' the parallel part of the algorithm. If \code{NULL} the number of cores used
#' will be equal to the number of cores available on the machine - 1.
#'
#' @return An object of class \code{\link{GWAS}}.
#'
#' @references Astle, William, and David J. Balding. 2009. Population Structure
#' and Cryptic Relatedness in Genetic Association Studies. Statistical Science 
#' 24 (4): 451–71. \url{https://doi.org/10.1214/09-sts307}.
#' @references Devlin, B., and Kathryn Roeder. 1999. Genomic Control for 
#' Association Studies. Biometrics 55 (4): 997–1004. 
#' \url{https://doi.org/10.1111/j.0006-341x.1999.00997.x}.
#' @references Kang et al. (2008) Efficient Control of Population Structure in
#' Model Organism Association Mapping. Genetics 178 (3): 1709–23. 
#' \url{https://doi.org/10.1534/genetics.107.080101}.
#' @references Millet, E. J., Pommier, C., et al. (2019). A multi-site 
#' experiment in a network of European fields for assessing the maize yield 
#' response to environmental scenarios [Data set]. 
#' \url{https://doi.org/10.15454/IASSTN}
#' @references Rincent et al. (2014) Recovering power in association mapping
#' panels with variable levels of linkage disequilibrium. Genetics 197 (1): 
#' 375–87. \url{https://doi.org/10.1534/genetics.113.159731}.
#' @references Segura et al. (2012) An efficient multi-locus mixed-model
#' approach for genome-wide association studies in structured populations.
#' Nature Genetics 44 (7): 825–30. \url{https://doi.org/10.1038/ng.2314}.
#' @references Sun et al. (2010) Variation explained in mixed-model association
#' mapping. Heredity 105 (4): 333–40. \url{https://doi.org/10.1038/hdy.2010.11}.
#' @references Tunnicliffe W. (1989) On the use of marginal likelihood in time
#' series model estimation. JRSS 51 (1): 15–27.
#' @references VanRaden P.M. (2008) Efficient methods to compute genomic
#' predictions. Journal of Dairy Science 91 (11): 4414–23. 
#' \url{https://doi.org/10.3168/jds.2007-0980}.
#' @references Brzyski D. et al. (2017) Controlling the Rate of GWAS False 
#' Discoveries. Genetics 205 (1): 61-75.
#' \url{https://doi.org/10.1534/genetics.116.193987 }
#' 
#' @examples 
#' ## Create a gData object Using the data from the DROPS project.
#' ## See the included vignette for a more extensive description on the steps.
#' data(dropsMarkers)
#' data(dropsMap)
#' data(dropsPheno)
#' ## Add genotypes as row names of dropsMarkers and drop Ind column.
#' rownames(dropsMarkers) <- dropsMarkers[["Ind"]]
#' dropsMarkers <- dropsMarkers[colnames(dropsMarkers) != "Ind"]
#' ## Add genotypes as row names of dropsMap.
#' rownames(dropsMap) <- dropsMap[["SNP.names"]]
#' ## Rename Chomosome and Position columns.
#' colnames(dropsMap)[match(c("Chromosome", "Position"), 
#'                    colnames(dropsMap))] <- c("chr", "pos")
#' ## Convert phenotypic data to a list.
#' dropsPhenoList <- split(x = dropsPheno, f = dropsPheno[["Experiment"]])
#' ## Rename Variety_ID to genotype and select relevant columns.
#' dropsPhenoList <- lapply(X = dropsPhenoList, FUN = function(trial) {
#'   colnames(trial)[colnames(trial) == "Variety_ID"] <- "genotype"
#'   trial <- trial[c("genotype", "grain.yield", "grain.number", "seed.size",
#'                  "anthesis", "silking", "plant.height", "tassel.height",
#'                  "ear.height")]
#' return(trial)
#' }) 
#' gDataDrops <- createGData(geno = dropsMarkers, map = dropsMap, 
#'                           pheno = dropsPhenoList)
#'                           
#' ## Run single trait GWAS for trait 'grain.yield' for trial Mur13W.
#'  \donttest{
#' GWASDrops <- runSingleTraitGwas(gData = gDataDrops,
#'                                trials = "Mur13W",
#'                                traits = "grain.yield")
#'  }
#'                                
#' ## Run single trait GWAS for trait 'grain.yield' for trial Mur13W.
#' ## Use chromosome specific kinship matrices calculated using vanRaden method.
#' \donttest{
#' GWASDropsMult <- runSingleTraitGwas(gData = gDataDrops,
#'                                     trials = "Mur13W",
#'                                     traits = "grain.yield",
#'                                     kinshipMethod = "vanRaden",
#'                                     GLSMethod = "multi")  
#' }
#'
#' @seealso \code{\link{GWAS}}, \code{\link{kinship}}, \code{\link{dropsData}}
#'
#' @importFrom data.table :=
#' @import data.table
#' @import stats
#'
#' @export
runSingleTraitGwas <- function(gData,
                               traits = NULL,
                               trials = NULL,
                               covar = NULL,
                               snpCov = NULL,
                               kin = NULL,
                               kinshipMethod = c("astle", "IBS", "vanRaden"),
                               remlAlgo = c("EMMA", "NR"),
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
                               nCores = NULL) {
  ## Checks.
  chkGData(gData)
  chkMarkers(gData$markers)
  chkTrials(trials, gData)
  ## If trials is null set trials to all trials in pheno.
  if (is.null(trials)) {
    trials <- seq_along(gData$pheno)
  }
  chkTraits(traits, trials, gData, multi = TRUE)
  chkCovar(covar, gData)
  ## If covar is given as numeric convert to character.
  if (is.numeric(covar)) {
    covar <- colnames(gData$covar)[covar]
  }
  chkSnpCov(snpCov, gData)
  GLSMethod <- match.arg(GLSMethod)
  chkKin(kin, gData, GLSMethod)
  kinshipMethod <- match.arg(kinshipMethod)
  remlAlgo <- match.arg(remlAlgo)
  chkNum(sizeInclRegion, min = 0)
  if (sizeInclRegion > 0) {
    chkNum(minR2, min = 0, max = 1)
  }
  if (useMAF) {
    chkNum(MAF, min = 0, max = 1)
    MAF <- max(MAF, 1e-6)
  } else {
    chkNum(MAC, min = 0)
    MAC <- max(MAC, 1)
  }
  thrType <- match.arg(thrType)
  if (thrType == "bonf") {
    chkNum(alpha, min = 0)
  } else if (thrType == "fixed") {
    chkNum(LODThr, min = 0)
  } else if (thrType == "small") {
    chkNum(nSnpLOD, min = 0)
  } else if (thrType == "fdr") {
    chkNum(alpha, min = 0)
    chkNum(rho, min = 0, max = 1)
    chkNum(pThr, min = 0, max = 1)
  }
  ## Compute kinship matrix (GSLMethod single)
  ## or kinship matrices per chromosome (GLSMethod multi).
  K <- computeKin(GLSMethod = GLSMethod, kin = kin, gData = gData,
                  markers = gData$markers, map = gData$map,
                  kinshipMethod = kinshipMethod)
  ## Compute max value in markers
  maxScore <- min(max(gData$markers, na.rm = TRUE), 2)
  ## Define data.frames for total results.
  GWATot <- signSnpTot <- varCompTot <- LODThrTot <- inflationFactorTot <-
    setNames(vector(mode = "list", length = length(trials)),
             names(gData$pheno[trials]))
  for (trial in trials) {
    ## Get traits for current trial.
    if (is.numeric(traits)) {
      ## If traits is given as numeric convert to character.
      traits <- colnames(gData$pheno[[trial]])[traits]
    } else if (is.null(traits)) {
      ## If no traits supplied extract them from pheno data.
      traits <- colnames(gData$pheno[[trial]])[-1]
    }
    ## Add covariates to phenotypic data.
    phExp <- expandPheno(gData = gData, trial = trial, covar = covar,
                         snpCov = snpCov)
    ## Get phenotypic data for trial from expanded phenotypic data.
    phTr <- phExp$phTr
    ## Restrict phenotypic data to genotypes that are actually in markers.
    phTr <- phTr[phTr[["genotype"]] %in% rownames(gData$markers), ]
    ## Get covariates for trial from expanded phenotypic data.
    covTr <- phExp$covTr
    ## Create vectors and lists for storing trait specific results.
    LODThrTr <- inflationFactorTr <-
      setNames(numeric(length = length(traits)), traits)
    GWATotTr <- signSnpTotTr <- varCompTr <-
      setNames(vector(mode = "list", length = length(traits)), traits)
    ## Perform GWAS for all traits.
    for (trait in traits) {
      ## Remove missings and select relevant columns only.
      phTrTr <- phTr[!is.na(phTr[trait]), c("genotype", trait, covTr)]
      ## If only NA values skip this trait.
      if (nrow(phTrTr) == 0) {
        next
      }
      ## Select genotypes where trait is not missing.
      nonMissRepId <- phTrTr[["genotype"]]
      nonMiss <- unique(nonMissRepId)
      ## Create a reduced map containing only markers that are in markers.
      mapRed <- gData$map[rownames(gData$map) %in% colnames(gData$markers), ]
      ## Get chromosomes from map that actually have at least one 
      ## corresponding marker in markers.
      chrs <- if (GLSMethod == "multi") {
        unique(mapRed[["chr"]])
      }
      ## Create a reduced marker file removing genotypes that are missing and
      ## markers that are not in map
      markersRed <- gData$markers[nonMiss, colnames(gData$markers) %in%
                                    rownames(mapRed), drop = FALSE]
      ## Compute allele frequencies based on genotypes for which phenotypic
      ## data is available.
      allFreq <- colMeans(markersRed, na.rm = TRUE) / maxScore
      if (!useMAF) {
        ## MAC used. Compute MAF from MAC.
        MAF <- MAC / (maxScore * length(nonMiss)) - 1e-5
      }
      ## Estimate variance components.
      vc <- estVarComp(GLSMethod = GLSMethod, remlAlgo = remlAlgo,
                       trait = trait, pheno = phTrTr, covar = covTr,
                       K = K, chrs = chrs, nonMiss = nonMiss, 
                       nonMissRepId = nonMissRepId)
      varCompTr[[trait]] <- vc$varComp
      vcovMatrix <- vc$vcovMatrix
      ## Define single column matrix with trait non missing values.
      y <- phTrTr[which(phTrTr$genotype %in% nonMiss), trait]
      ## Set up a data.table for storing results containing map info and
      ## allele frequencies.
      GWAResult <- data.table::data.table(trait = trait, snp = rownames(mapRed),
                                          mapRed, allFreq = allFreq, 
                                          key = "snp")
      if (GLSMethod == "single") {
        ## Determine segregating markers. Exclude snps used as covariates.
        segMarkers <- which(allFreq >= MAF & allFreq <= (1 - MAF))
        ## Exclude snpCovariates from segregating markers.
        exclude <- exclMarkers(snpCov = snpCov, markers = markersRed,
                               allFreq = allFreq)
        ## The following is based on the genotypes, not the replicates:
        X <- markersRed[nonMissRepId, setdiff(segMarkers, exclude)]
        Z <- if (length(covTr) > 0) {
          ## Define covariate matrix Z.
          as.matrix(phTrTr[which(phTrTr$genotype %in% nonMiss), covTr])
        }
        ## Compute pvalues and effects using fastGLS.
        GLSResult <- fastGLS(y = y, X = X, Sigma = vcovMatrix, covs = Z,
                             nCores = nCores)
        ## Merge GLSResult to GWAResult.
        GWAResult[GLSResult, names(GLSResult[, -1]) := GLSResult[, -1]]
        ## Compute p-values and effects for snpCovariates using fastGLS.
        for (snpCovariate in snpCov) {
          GLSResultSnpCov <- 
            fastGLS(y = y,
                    X = markersRed[nonMissRepId, snpCovariate, drop = FALSE],
                    Sigma = vcovMatrix,
                    covs = Z[, which(colnames(Z) != snpCovariate), 
                             drop = FALSE], nCores = nCores)
          ## Merge GLSResult for snp covariate to GWAResult.
          GWAResult[GLSResultSnpCov, 
                    names(GLSResultSnpCov[, -1]) := GLSResultSnpCov[, -1]]
        }
      } else if (GLSMethod == "multi") {
        ## Similar to GLSMethod single except using chromosome specific kinship
        ## matrices.
        for (chr in chrs) {
          mapRedChr <- mapRed[which(mapRed[["chr"]] == chr), ]
          markersRedChr <- markersRed[, colnames(markersRed) %in%
                                        rownames(mapRedChr), drop = FALSE]
          allFreqChr <- colMeans(markersRedChr, na.rm = TRUE) / maxScore
          ## Determine segregating markers. Exclude snps used as covariates.
          segMarkersChr <- which(allFreqChr >= MAF & allFreqChr <= (1 - MAF))
          ## Exclude snpCovariates from segregating markers.
          exclude <- exclMarkers(snpCov = snpCov, markers = markersRedChr,
                                 allFreq = allFreqChr)
          ## Remove excluded snps from segreg markers for current chromosome.
          segMarkersChr <- setdiff(segMarkersChr, exclude)
          ## If there are no segregating markers for current chromosome 
          ## continue with next chromosome.
          ## This is highly unlikely for real data.
          if (!length(segMarkersChr)) next
          X <- markersRedChr[nonMissRepId, segMarkersChr, drop = FALSE]
          Z <- if (length(covTr) > 0) {
            ## Define covariate matrix Z.
            as.matrix(phTrTr[which(phTrTr$genotype %in% nonMiss), covTr])
          }
          GLSResult <- fastGLS(y = y, X = X,
                               Sigma = vcovMatrix[[which(chrs == chr)]],
                               covs = Z, nCores = nCores)
          ## Merge GLSResult to GWAResult.
          GWAResult[GLSResult, names(GLSResult[, -1]) := GLSResult[, -1]]
          ## Compute pvalues and effects for snpCovariates using fastGLS.
          for (snpCovariate in intersect(snpCov, colnames(markersRedChr))) {
            GLSResultSnpCov <-
              fastGLS(y = y, X = markersRedChr[nonMissRepId, snpCovariate,
                                               drop = FALSE],
                      Sigma = vcovMatrix[[which(chrs == chr)]],
                      covs = Z[, which(colnames(Z) != snpCovariate),
                               drop = FALSE], nCores = nCores)
            ## Merge GLSResult for snp covariate to GWAResult.
            GWAResult[GLSResultSnpCov, 
                      names(GLSResultSnpCov[, -1]) := GLSResultSnpCov[, -1]]
          }
        }
      }
      ## Effects should be for a single allele, not for 2
      if (maxScore == 1) {
        GWAResult[, "effect" := GWAResult[["effect"]] / 2]
      }
      ## Calculate the genomic inflation factor.
      GC <- genCtrlPVals(pVals = GWAResult[["pValue"]], nObs = length(nonMiss),
                         nCov = length(covTr))
      inflationFactorTr[trait] <- GC$inflation
      ## Rescale p-values.
      if (genomicControl) {
        GWAResult[, "pValue" := GC$pValues]
      }
      ## Compute LOD score.
      GWAResult[, "LOD" := -log10(GWAResult[["pValue"]])]
      ## When thrType is bonferroni or small, determine the LOD threshold.
      if (thrType == "bonf") {
        ## Compute LOD threshold using Bonferroni correction.
        LODThr <- -log10(alpha / sum(!is.na(GWAResult[["pValue"]])))
      } else if (thrType == "small") {
        ## Compute LOD threshold by computing the 10log of the nSnpLOD item
        ## of ordered p values.
        LODThr <- sort(na.omit(GWAResult[["LOD"]]), decreasing = TRUE)[nSnpLOD]
      } else if (thrType == "fdr") {
        LODThr <- NA
      }
      LODThrTr[trait] <- LODThr
      ## Select the SNPs whose LOD-scores are above the threshold.
      if (thrType == "fdr") {
        signSnpTotTr[[trait]] <-
          extrSignSnpsFDR(GWAResult = GWAResult, markers = markersRed,
                          maxScore = maxScore, pheno = phTrTr, trait = trait,
                          rho = rho, pThr = pThr, alpha = alpha)
      } else {
        signSnpTotTr[[trait]] <-
          extrSignSnps(GWAResult = GWAResult, LODThr = LODThr,
                       sizeInclRegion = sizeInclRegion, minR2 = minR2,
                       map = mapRed, markers = markersRed,
                       maxScore = maxScore, pheno = phTrTr, trait = trait)
      }
      ## Sort columns.
      data.table::setkeyv(x = GWAResult, cols = c("trait", "chr", "pos"))
      GWATotTr[[trait]] <- GWAResult
    } # end for (trait in traits)
    ## Bind data together for results and significant SNPs.
    GWATot[[trial]] <- data.table::rbindlist(GWATotTr)
    signSnpTot[[trial]] <- data.table::rbindlist(signSnpTotTr)
    varCompTot[[trial]] <- varCompTr
    LODThrTot[[trial]] <- LODThrTr
    inflationFactorTot[[trial]] <- inflationFactorTr
  } # end for (trial in trials)
  ## No significant SNPs should return NULL instead of data.frame().
  signSnpTot <- lapply(signSnpTot, FUN = function(x) {
    if (is.null(x) || nrow(x) == 0) NULL else x
  })
  ## Collect info.
  GWASInfo <- list(call = match.call(),
                   remlAlgo = remlAlgo,
                   thrType = thrType,
                   MAF = MAF,
                   GLSMethod = GLSMethod,
                   varComp = varCompTot,
                   genomicControl = genomicControl,
                   inflationFactor = inflationFactorTot)
  return(createGWAS(GWAResult = GWATot,
                    signSnp = signSnpTot,
                    kin = K,
                    thr = LODThrTot,
                    GWASInfo = GWASInfo))
}

