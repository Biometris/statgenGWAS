#' Compute REML estimates of variance components using EMMA algorithm.
#'
#' Using the EMMA algorithm as in Kang et al. (2008) to compute REML estimates
#' of genetic and residual variance components.
#'
#' @param gData An object of class gData containing at least a data.frame
#' \code{pheno}. If \code{K} is not supplied a matrix \code{kinship} should be
#' in \code{gData}. If covariates are included then a data.frame covar is
#' needed as wel and if an extra snp is to be included as covariate (defined
#' in \code{snpName}) then a data.frame \code{markers} is also needed. Missing
#' values in \code{pheno} are allowed but will be excluded from the
#' calculations.
#' @param trait A trait for which to estimate variance components. This can be
#' either numeric index or character name of a column in \code{pheno}.
#' @param environment An environment for which to estimate variance components.
#' This can be either numeric index or character name of a list item in
#' \code{pheno}.
#' @param K An optional kinship matrix. If \code{NULL} then matrix
#' \code{kinship} in \code{gData} is used. If both \code{K} is provided and
#' \code{gData} contains a matrix \code{kinship} then \code{K} is used.
#' @param covar An optional vector of covariates taken into account when
#' estimating variance components. These can be either numeric indices or
#' character names of columns in \code{covar} in \code{gData}. If \code{NULL}
#' no covariates are used.
#' @param snpName An optional character string of a marker in \code{markers} in
#' \code{gData} to be included as covariate. If used the \code{gData} object
#' should contain a data.frame \code{markers}.
#' @param Z An optional incidence matrix mapping each observed phenotype to
#' one of inbred strains.
#' @param nGrids An integer indicating the number of intervals used for local
#' optimisation within the algorithm.
#' @param lLim A numerical value indicating the lower limit of the interval over
#' which optimisating is done.
#' @param uLim A numerical value indicating the upper limit of the interval over
#' which optimisating is done.
#' @param eps A numerical value used as computational tolerance in the
#' algorithm.
#'
#' @return A list with two components:
#' \itemize{
#' \item{\code{varcomp} a vector of genetic variance Vg and residual variance
#' Ve}
#' \item{\code{vcovMatrix} The variance covariance matrix corresponding to
#' the computed variances.}
#' }
#' @references Kang et al. (2008) (Efficient Control of Population Structure in
#' Model Organism Association Mapping. Genetics, March 2008, Vol. 178, no. 3,
#' p. 1709-1723
#'
#' @import stats
#'
#' @keywords internal
EMMA <- function(gData,
                 trait,
                 environment,
                 K = NULL,
                 covar = NULL,
                 snpName = NULL,
                 Z = NULL,
                 nGrids = 100,
                 lLim = -10,
                 uLim = 10,
                 eps = .Machine$double.eps ^ 0.25) {
  ## Checks.
  chkGData(gData, comps = "pheno")
  if (missing(environment) || length(environment) > 1 ||
      !(is.numeric(environment) || is.character(environment))) {
    stop("environment should be a single numeric or character.\n")
  }
  if ((is.character(environment) && !environment %in% names(gData$pheno)) ||
      (is.numeric(environment) && environment > length(gData$pheno))) {
    stop("environment should be a list item in pheno.\n")
  }
  chkTraits(trait, environment, gData, multi = FALSE)
  if (!is.null(K) && !(inherits(K, "Matrix") || is.matrix(K))) {
    stop("K should be a matrix.\n")
  }
  if (is.null(K) && is.null(gData$kinship)) {
    stop("gData contains no matrix kinship so K should be provided.\n")
  }
  chkCovar(covar, gData)
  if (!is.null(snpName) && (length(snpName) > 1 || !is.character(snpName))) {
    stop("snpName should be a single character.\n")
  }
  if (!is.null(Z) && !is.matrix(Z)) {
    stop("Z should be a matrix.\n")
  }
  chkNum(nGrids, min = 1)
  chkNum(lLim)
  chkNum(uLim)
  if (lLim >= uLim) {
    stop("lLim should be smaller than uLim.\n")
  }
  chkNum(eps, min = 0)
  ## Add column genotype to environment.
  phEnv <- gData$pheno[[environment]]
  ## Remove data with missings in trait or any of the covars.
  nonMiss <- phEnv$genotype[!is.na(phEnv[trait])]
  nonMissId <- which(!is.na(phEnv[trait]))
  if (!is.null(covar)) {
    misCov <- rownames(gData$covar)[rowSums(is.na(gData$covar[covar])) == 0]
    nonMiss <- nonMiss[nonMiss %in% misCov]
    nonMissId <- intersect(nonMissId, which(phEnv$genotype %in% misCov))
  }
  if (is.null(K)) {
    K <- gData$kinship[nonMiss, nonMiss]
  } else {
    K <- K[nonMiss, nonMiss]
  }
  y <- phEnv[nonMissId, trait]
  ## Define intercept.
  X <- matrix(data = 1, nrow = length(nonMiss), ncol = 1)
  if (!is.null(covar)) {
    ## Add covars to intercept.
    X <- cbind(X, as.matrix(gData$covar[nonMiss, covar, drop = FALSE]))
  }
  if (!is.null(snpName)) {
    ## Add extra snp to intercept + covars.
    X <- cbind(X, as.numeric(gData$markers[phEnv$genotype, snpName][nonMiss]))
  }
  ## Check resulting X for singularity.
  if (!is.matrix(try(solve(crossprod(X)), silent = TRUE))) {
    warning("X is singular.")
    return(list(varcomp = c(0, 0), K = K))
  }
  resEmma <- emmaCPP(y = y, k = K, x = X, eps = .Machine$double.eps ^ 0.25)
  vcovMatrix <- resEmma$vcovMatrix
  rownames(vcovMatrix) <- colnames(vcovMatrix) <- rownames(K)
  return(list(varComp = c(Vg = resEmma$maxVg, Ve = resEmma$maxVe),
              vcovMatrix = vcovMatrix))
}
