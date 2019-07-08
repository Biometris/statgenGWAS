#' S3 Class GWAS
#'
#' \code{createGWAS} creates an object of S3 class GWAS containing the results
#' of a GWAS analysis.
#' \code{GWAResult} and \code{signSnp} are both optional, however at least one
#' of those should be provided as input.\cr
#' \code{summary} and \code{plot} functions are available.
#'
#' @param GWAResult An optional data.frame or list of data.frames containing the
#' overall analysis results. Should a least contain columns \code{trait}, the
#' evaluated trait, \code{snp}, the name of the SNP,
#' \code{chr}, the chromosome number, \code{pos}, the position of the SNP on the 
#' chromosome, \code{pValue}, the p-values from the analysis and \code{LOD} 
#' the LOD-score.
#' @param signSnp An optional data.frame or list of data.frames containing
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
#' @keywords internal
createGWAS <- function(GWAResult = NULL,
                       signSnp = NULL,
                       kin = NULL,
                       thr = NULL,
                       GWASInfo = NULL) {
  ## Check that at least one input is provided.
  if (is.null(GWAResult) && is.null(signSnp)) {
    stop("At least one of GWAResult and signSnp should be provided.\n")
  }
  ## Check GWAResults
  if (!is.null(GWAResult)) {
    if (!is.data.frame(GWAResult) &&
        !(is.list(GWAResult) &&
          all(sapply(X = GWAResult, FUN = is.data.frame)))) {
      stop("GWAResult should be a data.frame or a list data.frames.\n")
    }
    if (is.data.frame(GWAResult)) {
      ## If not a list already put data.frame in a list.
      GWAResult <- list(GWAResult)
    }
    if (!all(sapply(GWAResult, FUN = function(x) {
      all(c("trait", "snp", "chr", "pos", "pValue", "LOD") %in%
          colnames(x))}))) {
      stop(paste("GWAResult should contain columns trait, snp, chr, pos,",
                 "pValue and LOD.\n"))
    }
  }
  ## Check signSnps
  if (!all(sapply(signSnp, FUN = is.null))) {
    if (!is.data.frame(signSnp) &&
        !(is.list(signSnp) && all(sapply(X = signSnp, FUN = function(x) {
          is.null(x) || is.data.frame(x)
        })))) {
      stop("signSnp should be a data.frame or a list of data.frames.\n")
    }
    if (is.data.frame(signSnp)) {
      ## If not a list already put data.frame in a list.
      signSnp <- list(signSnp)
    }
    if (!all(sapply(signSnp, FUN = function(x) {
      is.null(x) || all(c("trait", "snp", "snpStatus", "pValue", "LOD") %in%
                        colnames(x))}))) {
      stop(paste("signSnp should contain columns trait, snp, snpStatus,",
                 "pValue and LOD.\n"))
    }
  }
  ## Check kin
  if (!is.null(kin)) {
    if (!(inherits(kin, "Matrix") || is.matrix(kin)) &&
        !(is.list(kin) && all(sapply(X = kin, FUN = function(x) {
          inherits(x, "Matrix") || is.matrix(x)})))) {
      stop("kin should be a matrix or a list of matrices.\n")
    }
  }
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
#' @param environments A vector of strings or numeric indices indicating for
#' which environment the summary should be made. If \code{NULL}, a summary is
#' made for all environments.
#'
#' @export
summary.GWAS <- function(object, ..., environments = NULL) {
  ## Checks.
  if (!is.null(environments) && !is.character(environments) &&
      !is.numeric(environments)) {
    stop("environments should be a character or numeric vector.\n")
  }
  if ((is.character(environments) &&
       !all(environments %in% names(object$GWAResult))) ||
      (is.numeric(environments) &&
       !all(environments %in% 1:length(object$GWAResult)))) {
    stop("all environments should be in object.\n")
  }
  ## Convert character input to numeric.
  if (is.character(environments)) {
    environments <- which(names(object$GWAResult) == environments)
  }
  ## If NULL then summary of all environments.
  if (is.null(environments)) {
    environments <- 1:length(object$GWAResult)
  }
  for (environment in environments) {
    GWAResult <- object$GWAResult[[environment]]
    signSnp <- object$signSnp[[environment]]
    GWASInfo <- object$GWASInfo
    traits <- unique(GWAResult$trait)
    ## Print environment.
    cat(names(object$GWAResult)[environment], ":\n", sep = "")
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
            GWASInfo$varComp[[environment]][[trait]][1], "\n")
        cat("\t\tResidual variance:",
            GWASInfo$varComp[[environment]][[trait]][2], "\n\n")
      }
      if (!is.null(GWASInfo$thrType) && !is.null(GWASInfo$thrType)) {
        ## Print significant SNP info.
        cat("\t\tLOD-threshold:", object$thr[[environment]][trait], "\n")
        signSnpTrait <- signSnp[signSnp$trait == trait, ]
        if (!is.null(signSnpTrait)) {
          nSignSnp <-
            nrow(signSnpTrait[signSnpTrait$snpStatus == "significant snp", ])
          cat("\t\tNumber of significant SNPs:" , nSignSnp, "\n")
          if (nSignSnp > 0) {
            cat("\t\tSmallest p-value among the significant SNPs:",
                min(signSnpTrait[signSnpTrait$snpStatus == "significant snp",
                                 "pValue"]), "\n")
            cat("\t\tLargest p-value among the significant SNPs: ",
                max(signSnpTrait[signSnpTrait$snpStatus == "significant snp",
                                 "pValue"]),
                " (LOD-score: ",
                min(signSnpTrait[signSnpTrait$snpStatus == "significant snp",
                                 "LOD"]), ")\n\n", sep = "")
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
        cat("\t\tNo Genomic control correction was applied\n")
      }
      if (!is.null(GWASInfo$inflationFactor)) {
        cat("\t\tGenomic control inflation-factor:",
            round(GWASInfo$inflationFactor[[environment]][trait], 3), "\n\n")
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
#' large number of SNPs.}
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
#' phenotypic trait or environment. Optionally, vertical white lines can indicate
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
#' the data column to use for sorting. Default =
#' \code{FALSE}}
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
#' @param environment A character string or numeric index indicating for which
#' environment the plot should be made. If \code{x} only contains results for
#' one environment, \code{environment} may be \code{NULL}.
#' @param trait A character string indicating for which trait the results
#' should be plotted. For \code{type} "qtl" all traits are plotted. If \code{x}
#' only contains results for one trait, \code{trait} may be \code{NULL}.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE}, only a list of ggplot objects is invisibly returned.
#'
#' @references Millet et al. (2016) Genome-wide analysis of yield in Europe:
#' Allelic effects vary with drought and heat scenarios. Plant Physiology,
#' October 2016, Vol. 172, p. 749–764
#'
#' @export
plot.GWAS <- function(x,
                      ...,
                      plotType = c("manhattan", "qq", "qtl", "matrix"),
                      environment = NULL,
                      trait = NULL,
                      output = TRUE) {
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
  ## Checks.
  if (!is.null(environment) && !is.character(environment) &&
      !is.numeric(environment)) {
    stop("Environment should be a character or numerical value.\n")
  }
  if ((is.character(environment) && !environment %in% names(x$GWAResult)) ||
      (is.numeric(environment) && !environment %in% 1:length(x$GWAResult))) {
    stop("Environment should be in x.\n")
  }
  ## Convert character input to numeric.
  if (is.character(environment)) {
    environment <- which(names(x$GWAResult) == environment)
  }
  ## If NULL then summary of all environment.
  if (is.null(environment)) {
    if (length(x$GWAResult) != 1) {
      stop(paste("Environment not supplied but multiple environments",
                 "detected in data.\n"))
    } else {
      environment <- 1
    }
  }
  GWAResult <- x$GWAResult[[environment]]
  signSnp <- x$signSnp[[environment]]
  if (plotType != "qtl") {
    if (is.null(trait)) {
      trait <- unique(GWAResult$trait)
      if (length(trait) > 1) {
        if (substr(as.character(x$GWASInfo$call)[1], 1, 9) == "runSingle") {
          stop("Trait not supplied but multiple traits detected in data.\n")
        } else {
          ## For multi trait GWAS p-values are the same for all traits.
          trait <- trait[1]
        }
      }
    } else {
      GWAResult <- GWAResult[GWAResult$trait == trait, ]
      signSnp <- signSnp[signSnp$trait == trait, ]
    }
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
      set.seed(1234)
      mapShw <- map[(!is.na(map$LOD) & map$LOD >= dotArgs$lod) |
                      map$snp %in% dotArgs$effect, ]
      mapRem <- map[(!is.na(map$LOD) & map$LOD < dotArgs$lod) &
                      !map$snp %in% dotArgs$effect, ]
      sampIdx <- sample(seq_len(nrow(mapRem)), floor(nrow(mapRem) / 20),
                        prob = mapRem$LOD ^ 1.5)
      map <- rbind(mapRem[sampIdx, ], mapShw)
    }
    map <- merge(map, addPos, by = "chr")
    map$cumPos <- map$pos + map$add
    ## Extract row numbers for significant SNPs.
    if (!is.null(dotArgs$yThr)) {
      signSnpNr <- which(map$LOD > dotArgs$yThr)
    } else if (!is.null(signSnp)) {
      signSnpNr <- which(map$snp %in%
                           signSnp$snp[as.numeric(signSnp$snpStatus) == 1])
    } else {
      signSnpNr <- integer()
    }
    if (is.null(dotArgs$yThr)) {
      yThr <- x$thr[[environment]][trait]
    } else {
      yThr <- dotArgs$yThr
    }
    xEffects <- which(map$snp %in% dotArgs$effects)
    ## Create manhattan plot.
    do.call(manhattanPlot,
            args = c(list(xValues = map$cumPos, yValues = map$LOD,
                          map = map[, -which(colnames(map) == "LOD")],
                          xSig = signSnpNr, xEffects = xEffects,
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
      signSnp <- GWAResult[!is.na(GWAResult$LOD) &
                             GWAResult$LOD > dotArgs$yThr, ]
    }
    if (is.null(signSnp)) {
      stop("No significant SNPs found. No plot can be made.\n")
    }
    map <- GWAResult[!is.na(GWAResult$pos), c("chr", "pos")]
    if (!is.null(dotArgs$chr)) {
      signSnp <- signSnp[signSnp$chr %in% dotArgs$chr, ]
      map <- map[map$chr %in% dotArgs$chr, ]
      if (nrow(signSnp) == 0) {
        stop("Select at least one valid chromosome for plotting.\n")
      }
    }
    do.call(qtlPlot,
            args = c(list(data = signSnp, map = map, output = output),
                     dotArgs[!(names(dotArgs) %in% c("yThr", "chr"))]
                     ))
  } 
}
