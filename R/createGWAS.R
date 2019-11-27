#' S3 Class GWAS
#'
#' \code{createGWAS} creates an object of S3 class GWAS containing the results
#' of a GWAS analysis.
#' \code{GWAResult} and \code{signSnp} are both optional, however at least one
#' of those should be provided as input.\cr
#' \code{summary} and \code{plot} functions are available.
#'
#' @param GWAResult An optional data.table or list of data.frames containing the
#' overall analysis results. Should a least contain columns \code{trait}, the
#' evaluated trait, \code{snp}, the name of the SNP,
#' \code{chr}, the chromosome number, \code{pos}, the position of the SNP on the 
#' chromosome, \code{pValue}, the p-values from the analysis and \code{LOD} 
#' the LOD-score.
#' @param signSnp An optional data.table or list of data.frames containing
#' information on the significant SNPs and optionally the SNPs close to the
#' significant SNPs. Should at least contain columns \code{trait}, the
#' evaluated trait, \code{snp}, the name of the SNP, \code{pValue}, the p-values
#' from the analysis and \code{LOD} the LOD-score.
#' @param kin An optional kinship matrix or list of chromosome-specific kinship
#' matrices.
#' @param thr An optional numerical value, the threshold used in performing the
#' GWAS analysis.
#' @param GWASInfo a list containing extra information concering the GWAS
#' analysis.
#'
#' @return An object of class GWAS, a list of the input items.
#'
#' @seealso \code{\link{summary.GWAS}}, \code{\link{plot.GWAS}}
#'
#' @noRd
#' @keywords internal
createGWAS <- function(GWAResult = NULL,
                       signSnp = NULL,
                       kin = NULL,
                       thr = NULL,
                       GWASInfo = NULL) {
  ## Create GWAS object.
  GWAS <- structure(list(GWAResult = GWAResult,
                         signSnp = signSnp,
                         kinship = kin,
                         thr = thr,
                         GWASInfo = GWASInfo),
                    class = "GWAS")
  return(GWAS)
}

#' Summary function for the class \code{GWAS}
#'
#' Gives a summary for an object of S3 class \code{GWAS}.
#'
#' @param object An object of class \code{GWAS}
#' @param ... Not used
#' @param trials A vector of strings or numeric indices indicating for
#' which trials the summary should be made. If \code{NULL}, a summary is
#' made for all trials.
#'
#' @export
summary.GWAS <- function(object, 
                         ..., 
                         trials = NULL) {
  ## Checks.
  if (!is.null(trials) && !is.character(trials) &&
      !is.numeric(trials)) {
    stop("trials should be a character or numeric vector.\n")
  }
  if ((is.character(trials) &&
       !all(trials %in% names(object$GWAResult))) ||
      (is.numeric(trials) &&
       !all(trials %in% seq_along(object$GWAResult)))) {
    stop("All trials should be in object.\n")
  }
  ## Convert character input to numeric.
  if (is.character(trials)) {
    trials <- which(names(object$GWAResult) == trials)
  }
  ## If NULL then summary of all trials.
  if (is.null(trials)) {
    trials <- seq_along(object$GWAResult)
  }
  for (trial in trials) {
    GWAResult <- object$GWAResult[[trial]]
    signSnp <- object$signSnp[[trial]]
    GWASInfo <- object$GWASInfo
    traits <- unique(GWAResult$trait)
    ## Print trial.
    cat(names(object$GWAResult)[trial], ":\n", sep = "")
    ## Print traits.
    cat("\tTraits analysed:", paste(traits, collapse = ", "), "\n\n")
    ## Print SNP numbers.
    cat("\tData are available for", length(unique(GWAResult$snp)),
        "SNPs.\n")
    if (!is.null(GWASInfo$MAF)) {
      cat("\t", length(unique(GWAResult$snp[is.na(GWAResult$pValue)])),
          "of them were not analyzed because their minor allele frequency is",
          "below", GWASInfo$MAF, "\n\n")
    }
    for (trait in traits) {
      cat("\tTrait:", trait, "\n\n")
      if (substr(GWASInfo$call[[1]], 4, 4) == "S" &&
          !is.null(GWASInfo$GLSMethod) && GWASInfo$GLSMethod == "single") {
        ## Print mixed model info.
        cat("\t\tMixed model with only polygenic effects,",
            "and no marker effects:\n")
        cat("\t\tGenetic variance:",
            GWASInfo$varComp[[trial]][[trait]][1], "\n")
        cat("\t\tResidual variance:",
            GWASInfo$varComp[[trial]][[trait]][2], "\n\n")
      }
      if (!is.null(GWASInfo$thrType)) {
        ## Print significant SNP info.
        cat("\t\tLOD-threshold:", object$thr[[trial]][trait], "\n")
        signSnpTrait <- signSnp[signSnp[["snpStatus"]] == "significant snp" &
                                  signSnp[["trait"]] == 
                                  get("trait", envir = parent.frame(3))]
        if (!is.null(signSnpTrait)) {
          nSignSnp <- nrow(signSnpTrait)
          cat("\t\tNumber of significant SNPs:" , nSignSnp, "\n")
          if (nSignSnp > 0) {
            cat("\t\tSmallest p-value among the significant SNPs:",
                min(signSnpTrait[["pValue"]]), "\n")
            cat("\t\tLargest p-value among the significant SNPs: ",
                max(signSnpTrait[["pValue"]]),
                " (LOD-score: ", min(signSnpTrait[["LOD"]]), ")\n\n", sep = "")
          } else {
            cat("\n")
          }
        } else {
          cat("\t\tNo significant SNPs found.","\n\n")
        }
      }
      if (!is.null(GWASInfo$genomicControl) && GWASInfo$genomicControl) {
        ## Print genomic control.
        cat("\t\tGenomic control correction was applied\n")
      } else {
        cat("\t\tNo genomic control correction was applied\n")
      }
      if (!is.null(GWASInfo$inflationFactor)) {
        cat("\t\tGenomic control inflation-factor:",
            round(GWASInfo$inflationFactor[[trial]][trait], 3), "\n\n")
      }
    }
  }
}

#' Plot function for the class \code{GWAS}
#'
#' Creates a plot of an object of S3 class \code{GWAS}. The following types of
#' plot can be made:
#' \itemize{
#' \item{a manhattan plot, i.e. a plot of LOD-scores per SNP}
#' \item{a qq plot of observed LOD-scores versus expected LOD-scores}
#' \item{a qtl plot of effect sizes and directions for multiple traits}
#' }
#' Manhattan plots and qq plots are made for a single trait which
#' should be indicated using the parameter \code{trait} unless the analysis was
#' done for only one trait in which case it is detected automatically. The qtl
#' plot will plot all traits analysed.\cr
#' See details for a detailed description of the plots and the plot options
#' specific to the different plots.
#'
#' @section Manhattan Plot:
#' A LOD-profile of all marker positions and corresponding LOD-scores is
#' plotted. Significant markers are highlighted with red dots. By default these
#' are taken from the result of the GWAS analysis however the LOD-threshold for
#' significant parameters may be modified using the parameter \code{yThr}. The
#' treshold is plotted as a horizontal line. If there are previously known
#' marker effect, false positives and true negatives can also be marked.\cr
#' Extra parameter options:
#' \describe{
#' \item{\code{xLab}}{A character string, the x-axis label. Default =
#' \code{"Chromosomes"}}
#' \item{\code{yLab}}{A character string, the y-axis label. Default =
#' \code{-log10(p)}}
#' \item{\code{effects}}{A character vector, indicating which SNPs correspond
#' to a real (known) effect. Used for determining true/false positives and
#' false negatives. True positives are colored green, false positives orange and
#' false negatives yellow.}
#' \item{\code{colPalette}}{A color palette used for plotting. Default
#' coloring is done by chromosome, using black and grey.}
#' \item{\code{yThr}}{A numerical value for the LOD-threshold. The value from
#' the GWAS analysis is used as default.}
#' \item{\code{signLwd}}{A numerical value giving the thickness of the
#' points that are false/true positives/negatives. Default = 0.6}
#' \item{\code{lod}}{A positive numerical value. For the SNPs with a LOD-value
#' below this value, only 5\% is plotted. The chance of a SNP being plotting is
#' proportional to its LOD-score. This option can be useful when plotting a
#' large number of SNPs. The 5\% of SNPs plotted is selected randomly. For 
#' reproducible results use set.seed before calling the function.}
#' \item{\code{chr}}{A vector of chromosomes to be plotted. By default, all
#' chromosomes are plotted. Using this option allows restricting the plot to a
#' subset of chromosomes.}
#' }
#'
#' @section QQ Plot:
#' From the LOD-scores calculated in the GWAS analysis, a qq-plot is generated with
#' observed LOD-scores versus expected LOD-scores. Code is adapted from
#' Segura et al. (2012).
#'
#' @section QTL Plot:
#' A plot of effect sizes for the significant SNPs found in the GWAS analysis
#' is created. Each horizontal line contains QTLs of one trait,
#' phenotypic trait or trial. Optionally, vertical white lines can indicate
#' chromosome subdivision, genes of interest, known QTL, etc. Circle diameters
#' are proportional to the absolute value of allelic effect. Colors indicate the
#' direction of the effect: green when the allele increases the trait value,
#' and blue when it decreases the value.\cr
#' Extra parameter options:
#' \describe{
#' \item{\code{normalize}}{Should the snpEffect be normalized? Default =
#' \code{FALSE}}
#' \item{\code{sortData}}{Should the data be sorted before plotting? Either
#' \code{FALSE}, if no sorting should be done, or a character string indicating
#' the data column to use for sorting. This should be a numerical column.
#' Default = \code{FALSE}}
#' \item{\code{binPositions}}{An optional data.frame containing at leasts two
#' columns, chr(omosome) and pos(ition). Vertical lines are plotted at those
#' positions. Default = \code{NULL}}
#' \item{\code{printVertGrid}}{Should default vertical grid lines be plotted.
#' Default = \code{TRUE}}
#' \item{\code{yLab}}{A character string, the y-axis label. Default =
#' \code{"Traits"}}
#' \item{\code{yThr}}{A numerical value for the LOD-threshold. The value from
#' the GWAS analysis is used as default.}
#' \item{\code{chr}}{A vector of chromosomes to be plotted. By default all
#' chromosomes are plotted. Using this option this can be restricted to a
#' subset of chromosomes.}
#' \item{\code{exportPptx}}{Should the plot be exported to a .pptx file?
#' Default = \code{FALSE}}
#' \item{\code{pptxName}}{A character string, the name of the .pptx file to
#' which the plot is exported. Ignored if exportPptx = \code{FALSE}.}
#' }
#'
#' @param x An object of class \code{GWAS}.
#' @param ... further arguments to be passed on to the actual plotting
#' functions.
#' @param plotType A character string indicating the type of plot to be made.
#' One of "manhattan", "qq" and "qtl".
#' @param trial A character string or numeric index indicating for which
#' trial the plot should be made. If \code{x} only contains results for
#' one trial, \code{trial} may be \code{NULL}.
#' @param trait A character string indicating for which trait the results
#' should be plotted. For \code{type} "qtl" all traits are plotted. If \code{x}
#' only contains results for one trait, \code{trait} may be \code{NULL}.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE}, only a list of ggplot objects is invisibly returned.
#'
#' @references Millet et al. (2016) Genome-wide analysis of yield in Europe:
#' Allelic effects vary with drought and heat scenarios. Plant Physiology,
#' October 2016, Vol. 172, p. 749–764
#' @references Segura et al. (2012) An efficient multi-locus mixed-model
#' approach for genome-wide association studies in structured populations.
#' Nature Genetics, June 2012, Vol. 44, p. 825–830.
#'
#' @export
plot.GWAS <- function(x,
                      ...,
                      plotType = c("manhattan", "qq", "qtl"),
                      trial = NULL,
                      trait = NULL,
                      output = TRUE) {
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
  ## Checks.
  if (!is.null(trial) && !is.character(trial) &&
      !is.numeric(trial)) {
    stop("trial should be a character or numerical value.\n")
  }
  if ((is.character(trial) && !trial %in% names(x$GWAResult)) ||
      (is.numeric(trial) && !trial %in% seq_along(x$GWAResult))) {
    stop("trial should be in x.\n")
  }
  ## Convert character input to numeric.
  if (is.character(trial)) {
    trial <- which(names(x$GWAResult) == trial)
  }
  ## If NULL then plot of all trials.
  if (is.null(trial)) {
    if (length(x$GWAResult) != 1) {
      stop("trial not supplied but multiple trials detected in data.\n")
    } else {
      trial <- 1
    }
  }
  GWAResult <- x$GWAResult[[trial]]
  signSnp <- x$signSnp[[trial]]
  if (plotType != "qtl") {
    if (is.null(trait)) {
      trait <- unique(GWAResult$trait)
    }
    if (length(trait) > 1) {
      if (substr(as.character(x$GWASInfo$call)[1], 1, 9) == "runSingle") {
        stop("Trait not supplied but multiple traits detected in data.\n")
      } else {
        ## For multi trait GWAS p-values are the same for all traits.
        trait <- trait[1]
      }
    }
    GWAResult <- GWAResult[trait == as.name(trait), ]
    signSnp <- signSnp[trait == as.name(trait), ]
  }
  if (plotType == "manhattan") {
    ## Compute chromosome boundaries.
    GWAResult <- GWAResult[!is.na(GWAResult$pos), ]
    ## Select specific chromosome(s) for plotting.
    if (!is.null(dotArgs$chr)) {
      GWAResult <- GWAResult[GWAResult$chr %in% dotArgs$chr, ]
      if (nrow(GWAResult) == 0) {
        stop("Select at least one valid chromosome for plotting.\n")
      }
    }
    ## Get the boundary for each of the chromosomes. 
    ## Has to be converted to numeric to avoid integer overflow in the next step.
    chrBnd <- aggregate(x = GWAResult$pos, by = list(GWAResult$chr), 
                        FUN = function(p) {as.numeric(max(p))})
    ## Compute cumulative positions.
    addPos <- data.frame(chr = chrBnd[, 1],
                         add = c(0, cumsum(chrBnd[, 2]))[1:nrow(chrBnd)],
                         stringsAsFactors = FALSE)
    map <- GWAResult[, c("snp", "chr", "pos", "LOD")]
    if (!is.null(dotArgs$effects)) {
      if (!all(dotArgs$effects %in% map$snp)) {
        stop("All known effects should be in the map.\n")
      }
    }
    if (!is.null(dotArgs$lod)) {
      ## Of markers below lod only 5% will be plotted.
      ## A weighted sample is taken of those markers.
      ## The weigth of LOD ^ 1.5 empirically derived.
      lod <- dotArgs$lod
      chkNum(lod, min = 0)
      ## Split the map in a part that is always shown, i.e high enough lod
      ## And the remainder that is to be sampled.
      mapShw <- map[(!is.na(map$LOD) & map$LOD >= dotArgs$lod) |
                      map$snp %in% dotArgs$effect, ]
      mapRem <- map[(!is.na(map$LOD) & map$LOD < dotArgs$lod) &
                      !map$snp %in% dotArgs$effect, ]
      ## Take a sample of the number of rows in the remainder and 
      ## use as indices for actual sampling.
      sampIdx <- sample(seq_len(nrow(mapRem)), ceiling(nrow(mapRem) / 20),
                        prob = mapRem$LOD ^ 1.5)
      map <- rbind(mapRem[sampIdx, ], mapShw)
    }
    map <- merge(map, addPos, by = "chr")
    map$cumPos <- map$pos + map$add
    ## Extract row numbers for significant SNPs.
    if (!is.null(dotArgs$yThr)) {
      chkNum(dotArgs$yThr, min = 0)
      signSnpNr <- which(map$LOD > dotArgs$yThr)
    } else if (!is.null(signSnp)) {
      signSnpNr <- which(map$snp %in%
                           signSnp$snp[as.numeric(signSnp$snpStatus) == 1])
    } else {
      signSnpNr <- integer()
    }
    if (is.null(dotArgs$yThr)) {
      yThr <- x$thr[[trial]][trait]
    } else {
      yThr <- dotArgs$yThr
    }
    xEffects <- which(map$snp %in% dotArgs$effects)
    ## Create manhattan plot.
    do.call(manhattanPlot,
            args = c(list(xValues = map$cumPos, yValues = map$LOD,
                          map = map, xSig = signSnpNr, xEffects = xEffects,
                          chrBoundaries = chrBnd[, 2], yThr = yThr,
                          output = output),
                     dotArgs[!(names(dotArgs) %in% c("yThr", "lod", "chr",
                                                     "effects"))]
            ))
  } else if (plotType == "qq") {
    ## Create qq-plot.
    qqPlot(pValues = na.omit(GWAResult$pValue), ..., output = output)
  } else if (plotType == "qtl") {
    if (!is.null(dotArgs$yThr)) {
      chkNum(dotArgs$yThr, min = 0)
      yThr <- dotArgs$yThr
    } else {
      ## Get yThr from GWAS object.
      yThr <- x$thr[[trial]][[1]]
    }
    GWAResult$sign <- !is.na(GWAResult$LOD) & GWAResult$LOD > yThr
    if (!sum(GWAResult$sign)) {
      stop("No significant SNPs found. No plot can be made.\n")
    }
    map <- GWAResult[!is.na(GWAResult$pos), c("chr", "pos")]
    if (!is.null(dotArgs$chr)) {
      GWAResult <- GWAResult[GWAResult$chr %in% dotArgs$chr, ]
      map <- map[map$chr %in% dotArgs$chr, ]
      if (nrow(GWAResult) == 0) {
        stop("Select at least one valid chromosome for plotting.\n")
      }
    }
    do.call(qtlPlot,
            args = c(list(dat = GWAResult, map = map, output = output),
                     dotArgs[!(names(dotArgs) %in% c("yThr", "chr"))]
            ))
  } 
}
