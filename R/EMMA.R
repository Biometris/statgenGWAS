#' Compute REML estimates of variance components using EMMA algorithm.
#'
#' Using the EMMA algorithm as in Kang et al. (2008) to compute REML estimates
#' of genetic and residual variance components.
#'
#' @param dat A data.frame containing the phenotypic data on which the 
#' analysis is performed.
#' @param trait A trait for which to estimate variance components. This can be
#' either numeric index or character name of a column in \code{pheno}.
#' @param K An optional kinship matrix. If \code{NULL} then matrix
#' \code{kinship} in \code{gData} is used. If both \code{K} is provided and
#' \code{gData} contains a matrix \code{kinship} then \code{K} is used.
#' @param covar A data.frame of covariates taken into account when
#' estimating variance components. If \code{NULL} no covariates are used.
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
EMMA <- function(dat, 
                 trait,
                 K,
                 covar = NULL,
                 nGrids = 100,
                 lLim = -10,
                 uLim = 10,
                 eps = .Machine$double.eps ^ 0.25) {
  ## Remove data with missings in trait or any of the covars.
  nonMiss <- dat[!is.na(dat[trait]), "genotype"]
  nonMissId <- which(!is.na(dat[trait]))
  if (!is.null(covar)) {
    nonMissCov <- rownames(covar)[rowSums(is.na(covar)) == 0]
    nonMiss <- nonMiss[nonMiss %in% nonMissCov]
    nonMissId <- intersect(nonMissId, which(dat[["genotype"]] %in% nonMissCov))
  }
  K <- K[nonMiss, nonMiss]
  y <- dat[nonMissId, trait]
  ## Define intercept.
  X <- matrix(data = 1, nrow = length(nonMiss), ncol = 1)
  if (!is.null(covar)) {
    ## Add covars to intercept.
    X <- cbind(X, as.matrix(covar[nonMiss, , drop = FALSE]))
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
