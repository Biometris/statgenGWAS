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
#' @param sizeInclRegion An integer. Should the results for SNPs close to
#' significant SNPs be included? If so, the size of the region in centimorgan
#' or base pairs. Otherwise 0.
#' @param minR2 A numerical value between 0 and 1. Restricts the SNPs included
#' in the region close to significant SNPs to only those SNPs that are in
#' sufficient Linkage Disequilibruim (LD) with the significant snp, where LD
#' is measured in terms of \eqn{R^2}. If for example \code{sizeInclRegion} = 200000
#' and \code{minR2} = 0.5, then for every significant SNP also those SNPs whose
#' LD (\eqn{R^2}) with the significant SNP is at least 0.5 AND which are at
#' most 200kb away from this significant snp are included. Ignored if
#' \code{sizeInclRegion} = 0.
#' @param nCores A numerical value indicating the number of cores to be used by
#' the parallel part of the algorithm. If \code{NULL} the number of cores used
#' will be equal to the number of cores available on the machine - 1.
#'
#' @return An object of class \code{\link{GWAS}}.
#'
#' @references Astle W., Balding D. J. (2009) Population structure and cryptic
#' relatedness in genetic association studies, Stat. Sci., November 2009,
#' Vol. 24, no. 4, p. 451–471.
#' @references Devlin, B. and Roeder K. (1999) Genomic control for association
#' studies. Biometrics, December 1999, Vol. 55(4), p. 997-1004.
#' @references Kang et al. (2008) Efficient Control of Population Structure in
#' Model Organism Association Mapping. Genetics, March 2008, Vol. 178, no. 3,
#' p. 1709-1723.
#' @references Rincent et al. (2014) Recovering power in association mapping
#' panels with variable levels of linkage disequilibrium. Genetics, May 2014.
#' Vol. 197. p. 375–387.
#' @references Segura et al. (2012) An efficient multi-locus mixed-model
#' approach for genome-wide association studies in structured populations.
#' Nature Genetics, June 2012, Vol. 44, p. 825–830.
#' @references Sun et al. (2010) Variation explained in mixed-model association
#' mapping. Heredity, February 2010, Vol. 105, p. 333–340.
#' @references Tunnicliffe W. (1989) On the use of marginal likelihood in time
#' series model estimation. JRSS, Vol.51(1), p.15-27.
#' @references VanRaden P.M. (2008) Efficient methods to compute genomic
#' predictions. J Dairy Sci, November 2008, Vol. 91 p. 4414–4423.
#'
#' @seealso \code{\link{GWAS}}, \code{\link{kinship}}
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
                               thrType = c("bonf", "fixed", "small"),
                               alpha = 0.05 ,
                               LODThr = 4,
                               nSnpLOD = 10,
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
      }
      LODThrTr[trait] <- LODThr
      ## Select the SNPs whose LOD-scores are above the threshold.
      signSnpTotTr[[trait]] <-
        extrSignSnps(GWAResult = GWAResult, LODThr = LODThr,
                     sizeInclRegion = sizeInclRegion, minR2 = minR2,
                     map = mapRed, markers = markersRed,
                     maxScore = maxScore, pheno = phTrTr, trait = trait)
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

