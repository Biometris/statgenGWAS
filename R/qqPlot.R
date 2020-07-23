#' qq-plot of observed versus expected LOD-scores
#'
#' Given a vector of pvalues, generate a qq-plot of observed LOD-scores versus
#' expected LOD-scores. Code taken from Segura et al. and adapted
#'
#' @inheritParams manhattanPlot
#'
#' @param pValues A numeric vector of pValues. Missings are ignored when
#' plotting.
#'
#' @references Segura V, Vilhjálmsson BJ, Platt A, et al. An efficient
#' multi-locus mixed model approach for genome-wide association studies in
#' structured populations. Nature genetics. 2012;44(7):825-830.
#' doi:10.1038/ng.2314.
#'
#' @importFrom graphics plot
#'
#' @noRd
#' @keywords internal
qqPlot <- function(pValues,
                   ...,
                   output = TRUE) {
  dotArgs <- list(...)
  if (is.null(pValues) || !is.numeric(pValues) || any(pValues < 0) ||
      any(pValues > 1)) {
    stop("pValues should be an numeric vector with values between 0 and 1")
  }
  ## Remove missing values.
  pValues <- na.omit(pValues)
  ## Construct vector of expected values given the number of points.
  expected <- -log10(ppoints(n = length(pValues)))
  ## Convert p-Values to log scale.
  observed <- -log10(sort(pValues))
  ## Create a data.frame used as input for ggplot.
  plotDat <- data.frame(expected, observed)
  ## Construct title.
  if (!is.null(dotArgs$title)) {
    plotTitle <- dotArgs$title
  } else {
    plotTitle <- "QQ-plot"
  }
  p <- ggplot2::ggplot(plotDat, 
                       ggplot2::aes_string(x = "expected", y = "observed")) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "blue") +
    ## Center title.
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ## Add custom axis labels.
    ggplot2::labs(x = expression(Expected~~-log[10](p)), 
                  y = expression(Observed~~-log[10](p))) +
    ggplot2::ggtitle(plotTitle)
  if (output) {
    plot(p)
  }
  invisible(p)
}
