#' Draw a qtl plot
#'
#' Visualisation of the final set of QTLs selected after the multi-trait GWAS.
#'
#' Each horizontal line contains QTLs of one trait, phenotypic trait or
#' trial. Option: Vertical white lines can indicate chromosome
#' subdivision, genes of interest, known QTLs, etc. Circle diameters are
#' proportional to the absolute value of allelic effect. Colors indicate the
#' direction of the effect: green when the allele increases the trait value,
#' and blue when it decreases the value.
#'
#' @param dat A data.frame with QTL data to be plotted.
#' @param map A data.frame with at least the columns \code{chr}, the number of
#' the chromosome and \code{pos}, the position of the SNP on the chromosome.
#' These are used for calculating the physical limits of the chromosome.
#' @param normalize Should the snpEffect be normalized?
#' @param sortData Should the data be sorted before plotting? Either
#' \code{FALSE}, if no sorting should be done, or a character string indicating
#' the data column to use for sorting. This should be a numerical column.
#' @param binPositions An optional data.frame containing at leasts two columns,
#' chr and pos. Vertical lines are plotted at those positions.
#' @param printVertGrid Should default vertical grid lines be plotted.
#' @param yLab A character string, the y-axis label.
#' @param exportPptx Should the plot be exported to a .pptx file?
#' @param pptxName A character string, the name of the .pptx file to which the
#' plot is exported. Ignored if exportPptx = \code{FALSE}.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#' @param ... Further parameter.
#'
#' @return A ggplot object is returned and plotted on screen. Optionally the
#' plot is exported to pptx as well.
#'
#' @references Millet et al. (2016) Genome-wide analysis of yield in Europe:
#' Allelic effects vary with drought and heat scenarios. Plant Physiology,
#' October 2016, Vol. 172, p. 749–764
#'
#' @import ggplot2
#' @importFrom graphics plot
#'
#' @noRd
#' @keywords internal
qtlPlot <- function(dat,
                    map,
                    normalize = FALSE,
                    sortData = FALSE,
                    binPositions = NULL,
                    printVertGrid = TRUE,
                    yLab = "Traits",
                    exportPptx = FALSE,
                    pptxName = "",
                    output = TRUE,
                    ...) {
  ## Basic argument checks
  if (is.null(dat) || !is.data.frame(dat)) {
    stop("dat should be a data.frame.\n")
  }
  if (is.null(map) || !is.data.frame(map)) {
    stop("map should be a data.frame.\n")
  }
  if (is.null(normalize) || length(normalize) > 1 || !is.logical(normalize)) {
    stop("normalize should be a single logical value.\n")
  }
  if (is.null(sortData) || (is.logical(sortData) && sortData) ||
      (is.character(sortData) && length(sortData) > 1)) {
    stop("sortData should be either FALSE or a single character string.\n")
  }
  if (!is.null(binPositions) && (!is.data.frame(binPositions))) {
    stop("binPositions should be either NULL or an data.frame.\n")
  }
  if (is.null(exportPptx) || length(exportPptx) > 1 ||
      !is.logical(exportPptx)) {
    stop("exportPptx should be a single logical value.\n")
  }
  if (exportPptx) {
    if (is.null(pptxName) || length(pptxName) > 1 || !is.character(pptxName)) {
      stop("pptxName should be a single character string.\n")
    }
    if (tools::file_ext(pptxName) != "pptx") {
      stop(pptxName, " should have '.pptx' extension.\n")
    }
  }
  ## Check that all necessary columns are in the data
  reqCols <- c("chr", "trait", "effect", "pos")
  if (is.character(sortData)) reqCols <- c(reqCols, sortData)
  reqChk <- hasName(x = dat, name = reqCols)
  if (!all(reqChk)) {
    stop("dat lacks the following columns: ",
         paste0(reqCols[!reqChk], collapse = ", "), ".\n\n")
  }
  ## Check that sortData is a numerical column.
  if (is.character(sortData) && !is.numeric(dat[[sortData]])) {
    stop("sortData should be a numerical column.\n") 
  }
  ## Check that all necessary columns are in the bin file
  if (!is.null(binPositions)) {
    reqColsBin <- c("chr", "pos")
    reqChkBin <- hasName(x = binPositions, name = reqColsBin)
    if (!all(reqChkBin)) {
      stop("binPositions lacks the following columns: ",
           paste0(reqColsBin[!reqChkBin], collapse = ", "), ".\n\n")
    }
  } else {
    binPositions <- data.frame(pos = integer())
  }
  ## Reduce the allelic effect (because of the different units).
  ## This should be done on the full data to make sense.
  ## Therefore restricting to significant SNPs only is done afterwards.
  if (normalize) {
    datLst <- lapply(X = unique(dat$trait), FUN = function(x) {
      dat0 <- dat[dat[["trait"]] == as.character(x), ]
      if (nrow(dat0[!is.na(dat0[["effect"]]), ]) <= 1) {
        dat0$eff <- dat0[["effect"]]
      } else {
        dat0[!is.na(dat0[["effect"]]), "eff"] <- 
        scale(dat0[!is.na(dat0[["effect"]]), "effect"], center = FALSE,
              scale = apply(dat0[!is.na(dat0[["effect"]]), "effect", 
                                 drop = FALSE], 2, sd, na.rm = TRUE))
      }
      return(dat0)
    })
    dat <- do.call("rbind", datLst)
  } else dat$eff <- dat[["effect"]]
  ## After normalizing restrict the data to the significant SNPs only.
  dat <- dat[dat$sign, ]
  if (is.character(sortData)) {
    dat$sort <- dat[[sortData]]
  } else {
    dat$sort <- 1
  }
  ## Add the physical limits of the chromosomes, calculated from the map file
  ## This ensures plotting of all chromosomes
  limLow <- aggregate(x = map$pos, by = list(map$chr), FUN = min)
  limHigh <- aggregate(x = map$pos, by = list(map$chr), FUN = max)
  ## Empty data.frame with 2 lines per chromosomes and as many columns as the
  ## QTL data.frame
  lim <- data.frame(matrix(ncol = ncol(dat), nrow = 2 * nrow(limLow)))
  names(lim) <- names(dat)
  ## Trait and sort have to be filled. Value is not important
  lim["trait"] <- dat[1, "trait"]
  lim["sort"] <- dat[1, "sort"]
  ## Set eff to suppress warnings in printing. Setting it to -Inf ensures
  ## nothing is plotted
  lim$eff <- -Inf
  lim[, c("chr", "pos")] <- rbind(limLow, limHigh)
  dat <- rbind(dat, lim)
  ## Select relevant columns for plotting
  plotDat <- dat[, c("trait", "chr", "effect", "pos", "sort", "eff")]
  ## Add a column with the allelic effect direction (for points color)
  plotDat$color <- ifelse(plotDat$eff != -Inf,
                          ifelse(plotDat$eff > 0, "pos", "neg"), NA)
  ## Redefine trait as factor to keep traits in order as they appeared in
  ## original data when plotting.
  plotDat$trait <- factor(plotDat$trait, levels = unique(plotDat$trait))
  plotDat <- droplevels(plotDat)
  ## Create theme for plot
  qtlPlotTheme <- theme(plot.background = element_blank(),
                        panel.grid.major.x = 
                          element_line(color = ifelse(printVertGrid,
                                                      "white", "gray")),
                        panel.grid.major.y = element_line(color = "white"),
                        panel.grid.minor = element_blank(),
                        plot.title = element_text(size = 20, face = "bold", 
                                                  vjust = 2) ,
                        axis.ticks = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_blank(),
                        legend.position = "none",
                        panel.background = element_rect(fill = "gray"),
                        axis.text.x = element_text(size = 0),
                        axis.text.y = element_text(size = 13, color = "black"),
                        axis.title = element_text(size = 20, color = "navyblue"),
                        strip.background = element_rect(fill = "gray40"),
                        strip.text = element_text(),
                        strip.text.x = element_text(size = 14),
                        strip.text.y = element_text(size = 0))
  ## Create the plot object
  p <- ggplot(data = plotDat,
              ## Y data is sorted in reverse order because of the way ggplot plots.
              ## Point size proportional to allelic effect.
              ## Point color depends on the effect direction.
              aes_string(x = "pos", y = "reorder(trait, -sort)", 
                         size = "abs(eff)", color = "factor(color)")) +
    # add vertical lines at the bin positions
    geom_vline(aes_(xintercept = ~pos), data = binPositions,
               linetype = 1, color = "white") +
    ## Add the points with a slight transparency in case of overlap.
    geom_point(alpha = I(0.7), na.rm = TRUE) +
    ## Split of the plot according to the chromosomes on the x axis.
    ## Do not resize the x axis (otherwise every chromosome has the same size.
    ## Do not add extra space between two facets.
    ## Place the chromosome labels at the bottom.
    facet_grid(". ~ chr", scales = "free", space = "free", switch = "both") +
    ## Ascribe a color to the allelic direction (column 'color').
    scale_color_manual(labels = c("neg", "pos"), 
                       values = c("darkblue", "green4")) +
    ## Turn off clipping so points are plotted outside plot area.
    coord_cartesian(clip = "off") +
    ## use custom made theme
    qtlPlotTheme +
    labs(x = "Chromosomes", y = yLab) 
  if (exportPptx) {
    ## Save figure in .pptx
    if (requireNamespace("officer", quietly = TRUE) &&
        requireNamespace("rvg", quietly = TRUE)) {
      ## Create empty .pptx file
      pptOut <- officer::read_pptx()
      ## Add new slide (always necessary)
      pptOut <- officer::add_slide(x = pptOut, layout = "Title and Content",
                                   master = "Office Theme")
      ## Add plot to the document
      pptOut <- rvg::ph_with_vg_at(x = pptOut, ggobj = p, left = 0.9,
                                   top = 0.9, width = 8, height = 6.4)
      ## Add date to slide
      pptOut <- officer::ph_with(x = pptOut, value = format(Sys.Date()), 
                                 location = officer::ph_location_type(type = "dt"))
      ##Write .pptx
      print(pptOut, target = pptxName)
    } else {
      message(paste("Package officer needs to be installed to be able",
                    "to export to .pptx"))
    }
  }
  if (output) {
    plot(p)
  }
  invisible(p)
}



