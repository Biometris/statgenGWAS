#' Create a manhattan plot
#'
#' Given vectors of marker positions and corresponding LOD-scores plot a
#' LOD-profile. Significant markers can be highlighted with red dots. If there
#' are previously known marker effect also false positives and true negatives
#' will be marked.
#'
#' @param xValues A vector of cumulative marker positions.
#' @param yValues A vector of LOD-scores.
#' @param map A data.frame with at least the columns chr, the number of the
#' chromosome and cumPos, the cumulative position of the SNP the cumulative
#' position of the SNP starting from the first chromosome.
#' @param xLab A character string, the x-axis label.
#' @param yLab A character string, the y-axis label.
#' @param xSig A vector of integers, indicating which components in the vectors
#' xValues and yValues are significant.
#' @param xEffects A vector of integers, indicating which components in the
#' vector xValues correspond to a real (known) effect.
#' @param colPalette A color palette used for plotting.
#' @param yThr A numerical value for the LOD-threshold.
#' @param signLwd A numerical value giving the thickness of the
#' points that are false/true positives/negatives.
#' @param ... Other graphical parameters passed on to the actual plot function.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @return A LOD-profile with LOD-scores per snip. Markers declared significant
#' get a red dot, markers with a real effect get a blue dot. If both significant
#' and real effects are given false positives get an orange dot, true negatives
#' a yellow dot and true positives a green dot.
#'
#' @importFrom graphics plot
#'
#' @noRd
#' @keywords internal
manhattanPlot <- function(xValues,
                          yValues,
                          map,
                          xLab = "Chromosomes",
                          yLab = expression(-log[10](p)),
                          xSig = integer(),
                          xEffects = integer(),
                          colPalette = rep(c("black", "grey50"), 50),
                          yThr = NULL,
                          signLwd = 0.6,
                          ...,
                          output = TRUE) {
  ## Basic argument checks
  if (is.null(xValues) || !is.numeric(xValues)) {
    stop("xValues should be a numerical vector.\n")
  }
  if (is.null(yValues) || !is.numeric(yValues)) {
    stop("yValues should be a numerical vector.\n")
  }
  ## Check correspondence xValues and yValues
  if (length(xValues) != length(yValues)) {
    stop("xValues and yValues should be of the same length.\n")
  }
  if (is.null(xSig) || !is.numeric(xSig) || !all(xSig == round(xSig))) {
    stop("xSig should be an integer vector.\n")
  }
  if (is.null(xEffects) || !is.numeric(xEffects) || 
      !all(xEffects == round(xEffects))) {
    stop("xEffects should be an integer vector.\n")
  }
  if (!is.null(yThr)) {
    chkNum(yThr, min = 0)
  }
  chkNum(signLwd, min = 0)
  ## Extract central chromosome positions from map.
  ## Differentiate cases to deal with character chromosomes.
  if (is.numeric(map$chr)) {
    chrs <- as.numeric(levels(as.factor(map$chr)))
  } else {
    chrs <- levels(as.factor(map$chr))
  }
  xMarks <- aggregate(x = map$cumPos, by = list(map$chr), FUN = function(x) {
    min(x) + (max(x) - min(x)) / 2
  })[, 2]
  ## ylim has to be set to cope with subsetting based on lod.
  plotDat <- data.frame(x = xValues, y = yValues, chr = as.factor(map$chr))
  yMax <- max(c(yValues, yThr), na.rm = TRUE)
  p <- ggplot2::ggplot(data = if (length(c(xSig, xEffects)) == 0) {
    plotDat
  } else {
    plotDat[-c(xSig, xEffects),  ]
  }, ggplot2::aes_string(x = "x", y = "y", color = "chr")) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::scale_y_continuous(limits = c(0, yMax),
                                expand = c(0, 0, 0.1, 0)) +
    ggplot2::scale_x_continuous(breaks = xMarks, labels = chrs, limits = c(0, NA), 
                                expand = c(0, 0)) +
    ggplot2::scale_color_manual(values = colPalette, labels = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = xLab, y = yLab) +
    ggplot2::theme(legend.position = "none",
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank())
  if (length(xSig) > 0 && length(xEffects) == 0) {
    p <- p + ggplot2::geom_point(data = plotDat[xSig, , drop = FALSE], 
                                 color = "red")
  } else if (length(xSig) == 0 && length(xEffects) > 0) {
    p <- p + ggplot2::geom_point(data = plotDat[xEffects, , drop = FALSE], 
                                 color = "blue")
  } else if (length(xSig) > 0 && length(xEffects) > 0) {
    falsePos <- setdiff(xSig, xEffects)
    trueNeg <- setdiff(xEffects, xSig)
    truePos <- intersect(xSig, xEffects)
    p <- p + ggplot2::geom_point(data = plotDat[falsePos, , drop = FALSE],
                                 color = "orange", na.rm = TRUE) +
      ggplot2::geom_point(data = plotDat[trueNeg, , drop = FALSE],
                          color = "yellow", na.rm = TRUE) +
      ggplot2::geom_point(data = plotDat[truePos, , drop = FALSE],
                          color = "green", na.rm = TRUE)
  }
  if (!is.null(yThr)) {
    ## na.rm = TRUE is needed for plotting the results of GWAS with 
    ## thrType = "fdr".
    ## This results in several local thresholds that cannot be displayed as 
    ## a single line.
    p <- p + ggplot2::geom_hline(yintercept = yThr, linetype = 2, na.rm = TRUE)
  }
  if (output) {
    plot(p)
  }
  invisible(p)
}
