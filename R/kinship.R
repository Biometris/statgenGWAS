#' Functions for calculating kinship matrices
#'
#' A collection of functions for calculating kinship matrices using different
#' algorithms. The following algorithms are included: astle (Astle and Balding,
#' 2009), Identity By State (IBS) and VanRaden (VanRaden, 2008) for
#' marker matrices.
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
#' @references Astle W., Balding D. J. (2009) Population structure and cryptic
#' relatedness in genetic association studies, Stat. Sci., November 2009,
#' Vol. 24, no. 4, p. 451–471.
#' @references VanRaden P.M. (2008) Efficient methods to compute genomic
#' predictions. J Dairy Sci, November 2008, Vol. 91 p. 4414–4423.
#'
#' @examples X <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1), nrow = 3)
#' kinship(X, method = "astle")
#' kinship(X, method = "IBS")
#' kinship(X, method = "vanRaden")
#'
#' @export
kinship <- function(X,
                    method = c("astle", "IBS", "vanRaden"),
                    denominator = NULL) {
  method = match.arg(method);
  if (!is.null(denominator)) {
    chkNum(denominator, min = 0)
  }
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  K <- do.call(what = paste0(method, "CPP"),
               args = list(x = X, denom = denominator))
  rownames(K) <- colnames(K) <- rownames(X)
  return(K)
}
