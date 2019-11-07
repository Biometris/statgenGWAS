#' Fast Generalized Least Squares algoritm
#'
#' Compute statistics for the Generalized Least Squares (GLS) F-test based on
#' the algorithm proposed by Segura (2012). Also the \eqn{R_LR^2} statistics
#' as in Sun (2010) is computed.
#'
#' @inheritParams runSingleTraitGwas
#'
#' @param y A numeric vector of length n of phenotypic scores. No missing
#' values allowed.
#' @param X An n x m matrix of marker-scores, n being the number of
#' individuals, m the number of markers. no missing values allowed.
#' @param Sigma An n x n covariance matrix. No missing values allowed.
#' @param covs An n x c matrix of covariates (NOT including an intercept).
#' No missing values allowed.
#'
#' @return A data.table with the following columns:
#' \itemize{
#' \item{\code{pValue} p-values for the GLS F-test}
#' \item{\code{beta} effect sizes}
#' \item{\code{betaSe} standard errors of the effect sizes}
#' \item{\code{RLR2} L_LR^2 statistics as defined in Sun et al.}
#' \item{\code{rn} SNP name}
#' }
#'
#' @references Segura et al. (2012) An efficient multi-locus mixed-model
#' approach for genome-wide association studies in structured populations.
#' Nature Genetics, June 2012, Vol. 44, p. 825–830.
#' @references Sun et al. (2010) Variation explained in mixed-model association
#' mapping. Heredity, February 2010, Vol. 105, p. 333–340.
#'
#' @keywords internal
fastGLS <- function(y,
                    X,
                    Sigma,
                    covs = NULL,
                    nCores = NULL) {
  resCpp <- fastGLSCPP(X, y, Sigma, covs, nCores = nCores)
  ## Set row and column names to output merging results.
  rownames(resCpp) <- colnames(X)
  colnames(resCpp) <- c("pValue", "effect", "effectSe", "RLR2")
  ## Convert output to data.table.
  resCpp <- as.data.table(resCpp, keep.rownames = TRUE)
  setkeyv(resCpp, cols = "rn")
  return(resCpp)
}
