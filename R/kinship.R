#' Functions for calculating kinship matrices
#'
#' A collection of functions for calculating kinship matrices using different
#' algorithms. The following algorithms are included: astle (Astle and Balding,
#' 2009), Identity By State (IBS) and VanRaden (VanRaden, 2008) for
#' marker matrices. For method identity an identity kinship matrix is returned.
#'
#' @section Marker matrices:
#' In all algorithms the input matrix \code{X} is first cleaned, i.e. markers
#' with a variance of 0 are excluded from the calculation of the kinship matrix.
#' Then some form of scaling is done which differs per algorithm. This gives a
#' scaled matrix \code{Z}. The matrix \eqn{ZZ^t / denominator} is returned.
#' By default the denominator is equal to the number of columns in \code{Z} for
#' \code{astle} and \code{IBS} and \eqn{2 * p * (1-p)} where
#' \eqn{p = colSums(X) / (2 * nrow(X))} for \code{vanRaden}. This denominator
#' can be overwritten by the user, e.g. when computing kinship matrices by
#' splitting \code{X} in smaller matrices and then adding the results together
#' in the end.
#'
#' @param X An n x m marker matrix with genotypes in the rows (n) and markers in
#' the columns (m).
#' @param method The method used for computing the kinship matrix. 
#' @param denominator A numerical value. See details.
#'
#' @return An n x n kinship matrix.
#'
#' @references Astle, William, and David J. Balding. 2009. “Population Structure
#' and Cryptic Relatedness in Genetic Association Studies.” Statistical Science
#' 24 (4): 451–71. \url{https://doi.org/10.1214/09-sts307}.
#' @references VanRaden P.M. (2008) Efficient methods to compute genomic
#' predictions. Journal of Dairy Science 91 (11): 4414–23. 
#' \url{https://doi.org/10.3168/jds.2007-0980}.
#'
#' @examples 
#' ## Create example matrix.
#' M <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1), nrow = 3)
#' 
#' ## Compute kinship matrices using different methods.
#' kinship(M, method = "astle")
#' kinship(M, method = "IBS")
#' kinship(M, method = "vanRaden")
#' 
#' ## Compute kinship matrix using astle and balding method with denominator 2.
#' kinship(M, method = "astle", denominator = 2)
#'
#' @export
kinship <- function(X,
                    method = c("astle", "IBS", "vanRaden", "identity"),
                    denominator = NULL) {
  method = match.arg(method)
  chkMarkers(X)
  if (!is.null(denominator)) {
    chkNum(denominator, min = 0)
  }
  if (method == "identity") {
    K <- diag(nrow = nrow(X), ncol = nrow(X))
  } else {
    K <- do.call(what = paste0(method, "CPP"),
                 args = list(x = X, denom = denominator))
  }
  rownames(K) <- colnames(K) <- rownames(X)
  return(K)
}
