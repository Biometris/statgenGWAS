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
#' @return A data.frame with the following columns:
#' \itemize{
#' \item{\code{pValue} p-values for the GLS F-test}
#' \item{\code{beta} effect sizes}
#' \item{\code{betaSe} standard errors of the effect sizes}
#' \item{\code{RLR2} L_LR^2 statistics as defined in Sun et al.}
#' }
#'
#' @references Segura et al. (2012) An efficient multi-locus mixed-model
#' approach for genome-wide association studies in structured populations.
#' Nature Genetics, June 2012, Vol. 44, p. 825–830.
#' @references Sun et al. (2010) Variation explained in mixed-model association
#' mapping. Heredity, February 2010, Vol. 105, p. 333–340.
#'
#' @import stats
#'
#' @keywords internal
fastGLS <- function(y,
                    X,
                    Sigma,
                    covs = NULL,
                    nCores = NULL) {
  ## Check class and missing values.
  if (missing(y) || !(inherits(y, "Matrix") || is.numeric(y)) || anyNA(y)) {
    stop("y should be a numeric vector without missing values.\n")
  }
  #  if (missing(X) || !(inherits(X, "Matrix") || is.matrix(X)) || anyNA(X))
  #    stop("X should be a matrix without missing values.")
  if (missing(Sigma) || !(inherits(Sigma, "Matrix") || is.matrix(Sigma)) ||
      anyNA(Sigma)) {
    stop("Sigma should be a matrix without missing values.\n")
  }
  if (!is.null(covs) && (!(inherits(covs, "Matrix") || is.matrix(covs)) ||
                         anyNA(covs))) {
    stop("covs should be a numeric vector without missing values.\n")
  }
  n <- length(y)
  ## Check dimensions.
  if (nrow(X) != n) {
    stop(paste("The number of elements in y should be identical to the",
               "number of rows in X.\n"))
  }
  if (nrow(Sigma) != n || ncol(Sigma) != n) {
    stop(paste("The number of elements in y should be identical to the",
               "number of rows and columns in Sigma.\n"))
  }
  if (!is.null(covs) && nrow(covs) != n) {
    stop(paste("The number of elements in y should be identical to the",
               "number of rows in covs.\n"))
  }
  ## If necessary convert input to matrix
  resCpp <- fastGLSCPP(as.matrix(X), y, as.matrix(Sigma), covs, nCores = nCores)
  ## Construct output data.frame.
  GLS <- data.frame(pValue = resCpp$pVal,
                    beta = resCpp$beta,
                    betaSe = resCpp$betaSe,
                    RLR2 = resCpp$RLR2)
  rownames(GLS) <- colnames(X)
  return(GLS)
}
