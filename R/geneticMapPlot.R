geneticMapPlot <- function(map, 
                           highlight = NULL,
                           title = NULL,
                           ...,
                           output = TRUE) {
  ## Input checks.
  if (!is.null(highlight) && !inherits(highlight, "data.frame")) {
    stop("highlight should be a data.frame.\n")
  }
  if (!is.null(highlight)) {
    reqCols <- c("chr", "pos")
    missCols <- reqCols[!sapply(X = reqCols, FUN = function(reqCol) {
      hasName(x = highlight, name = reqCol)
    })]
    if (length(missCols) > 0) {
      stop("The following columns are missing in highlight:\n",
           paste(missCols, collapse = ", "))
    }
  }
  ## Convert chr to factor for plotting. Keep chr order.
  if (!is.factor(map[["chr"]])) {
    map[["chr"]] <- factor(map[["chr"]], levels = unique(map[["chr"]]))
  }
  ## Get chromosome lengths.
  maxPos <- tapply(X = map[["pos"]], INDEX = map[["chr"]], FUN = max)
  chrLen <- data.frame(chr = unique(map[["chr"]]), min = 0, max = maxPos)
  ## Construct title.
  if (is.null(title)) {
    title <- "Genetic map"
  }
  ## Add name to annotation points if not present in data.
  if (!is.null(highlight)) { 
    if (!hasName(x = highlight, name = "name")) {
      highlight[["name"]] <- paste0(highlight[["chr"]], "@", highlight[["pos"]])
    }
    highlight[["name"]] <- paste("\u2190", highlight[["name"]])
  }
  p <- ggplot2::ggplot(data = map,
                       ggplot2::aes_string(x = "chr", y = "pos")) +
    ggplot2::geom_segment(data = chrLen,
                          ggplot2::aes_string(x = "chr", y = "min", 
                                              yend = "max", xend = "chr"),
                          lineend = "round", size = 8, color = "grey90") +
    ggplot2::geom_segment(data = chrLen,
                          ggplot2::aes_string(x = "chr", y = "min", 
                                              yend = "max", xend = "chr")) +
    ggplot2::geom_point(pch = "_", size = 4) +
    ggplot2::labs(title = title, y = "Position", x = "Chromosome") +
    ggplot2::scale_y_reverse() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(),
      plot.title = ggplot2::element_text(hjust = 0.5))
  if (!is.null(highlight)) {
    p <- p +
      ggplot2::geom_point(data = highlight, color = "red", pch = "_", size = 5) +
      ggplot2::geom_text(data = highlight,
                         ggplot2::aes_string(x = "chr", y = "pos", 
                                             label = "name"),
                         color = "red", size = 2, nudge_x = 0.4)
  }
  if (output) {
    plot(p)
  }
  invisible(p)
}