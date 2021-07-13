#' S3 Class gData
#'
#' \code{createGData} creates an object of S3 class gData with genotypic and
#' phenotypic data for usage in further analysis. All input to the function is
#' optional, however at least one input should be provided. It is possible to
#' provide an existing \code{gData} object as additional input in which case
#' data is added to this object. Existing data will be overwritten with a
#' warning.
#'
#' @param gData An optional gData object to be modified. If \code{NULL}, a new
#' gData object is created.
#' @param geno A matrix or data.frame with genotypes in the rows and markers in
#' the columns. A matrix from the \code{matrix} in the base package may be
#' provided as well as as matrix from the Matrix package.\cr
#' If no row names are provided, they are taken from \code{pheno} (if supplied 
#' and dimension matches). If no column names are provided, the row names
#' from \code{map} are used (if supplied and dimension matches).
#' @param map A data.frame with columns \code{chr} for chromosome and
#' \code{pos} for position. Positions can be in base pair (bp) or centimorgan (cM). They
#' should not be cumulative over the chromosomes. Other columns are ignored.
#' Marker names should be in the row names. These should match the marker names
#' in \code{geno} (if supplied).
#' @param kin A kinship matrix or list of kinship matrices with genotype in
#' rows and colums. These matrices can be from the \code{matrix} class, as
#' defined in the base package, or from the \code{dsyMatrix} class, the class
#' of symmetric matrices in the Matrix package.\cr
#' The genotypes should be identical to the genotypes in \code{geno}.\cr
#' If a list of kinship matrices is provided these are supposed to be
#' chromosome specific matrices. In that case their names should match
#' the names of the chromosomes in \code{map}. If no names are
#' provided, the number of matrices should match the number of chromosomes
#' in \code{map} in which case default names are provided.
#' @param pheno A data.frame or a list of data.frames with phenotypic data,
#' with genotypes in the first column \code{genotype} and traits in the
#' following columns. The trait columns should be numerical columns only.
#' A list of data.frames can be used for replications, i.e. different
#' trials.
#' @param covar A data.frame with extra covariates per genotype. Genotypes
#' should be in the rows.
#'
#' @return An object of class \code{gData} with the following components:
#' \item{\code{map}}{a data.frame containing map data. Map is sorted by
#' chromosome and position.}
#' \item{\code{markers}}{a matrix containing marker information.}
#' \item{\code{pheno}}{a list of data.frames containing phenotypic data.}
#' \item{\code{kinship}}{a kinship matrix.}
#' \item{\code{covar}}{a data.frame with extra covariates.}
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{summary.gData}}
#'
#' @examples 
#' set.seed(1234)
#' ## Create genotypic data.
#' geno <- matrix(sample(x = c(0, 1, 2), size = 15, replace = TRUE), nrow = 3)
#' dimnames(geno) <- list(paste0("G", 1:3), paste0("M", 1:5))
#'
#' ## Construct map.
#' map <- data.frame(chr = c(1, 1, 2, 2, 2), pos = 1:5,
#'                   row.names = paste0("M", 1:5))
#'
#' ## Compute kinship matrix.
#' kin <- kinship(X = geno, method = "IBS")
#'
#' ## Create phenotypic data.
#' pheno <- data.frame(paste0("G", 1:3),
#'                     matrix(rnorm(n = 12, mean = 50, sd = 5), nrow = 3),
#'                     stringsAsFactors = FALSE)
#' dimnames(pheno) = list(paste0("G", 1:3), c("genotype", paste0("T", 1:4)))
#'
#' ## Combine all data in gData object.
#' gData <- createGData(geno = geno, map = map, kin = kin, pheno = pheno)
#' summary(gData)
#'
#' ## Construct covariate.
#' covar <- data.frame(C1 = c("a", "a", "b"), row.names = paste0("G", 1:3))
#'
#' ## Compute alternative kinship matrix.
#' kin2 <- kinship(X = geno, method = "astle")
#'
#' ## Add covariates to previously created gData object and overwrite
#' ## current kinship matrix by newly computed one.
#' gData2 <- createGData(gData = gData, kin = kin2, covar = covar)
#'
#' @name gData
NULL

#' @rdname gData
#'
#' @export
createGData <- function(gData = NULL,
                        geno = NULL,
                        map = NULL,
                        kin = NULL,
                        pheno = NULL,
                        covar = NULL) {
  ## Check gData
  if (!is.null(gData) && !inherits(gData, "gData")) {
    stop("Provided gData object should be of class gData.\n")
  }
  ## Check that at least one input argument, other than gData, is provided.
  if (is.null(geno) && is.null(map) && is.null(kin) && is.null(pheno) &&
      is.null(covar)) {
    stop("At least one of geno, map, kin, pheno and covar should be", 
         "provided.\n")
  }
  ## Modify map.
  if (!is.null(map)) {
    if (!is.data.frame(map)) {
      stop("map should be a data.frame.\n")
    }
    if (!all(hasName(x = map, name = c("chr", "pos")))) {
      ## chr and pos are obligatory cols.
      stop("chr and pos should be columns in map.\n")
    }
    if (!is.numeric(map$pos)) {
      stop("pos should be a numeric column in map.\n")
    }
    if (!is.null(gData$map)) {
      ## gData already contained a map object. Overwrite with a warning.
      warning("Existing map will be overwritten.\n", call. = FALSE)
    }
    ## Extract columns and order.
    map <- map[order(map[["chr"]], map[["pos"]]), c("chr", "pos")]
    if (all(rownames(map) %in% as.character(1:nrow(map)))) {
      ## If no marker name in input compute them from chromosome and position.
      ## Names are made unique if necessary by adding a suffix _1, _2, etc.
      replicates <- aggregate(chr ~ chr + pos, data = map, FUN = length)$chr
      suffix <- unlist(sapply(X = replicates, FUN = function(n) {
        if (n == 1) {
          return("")
        } else {
          return(paste0("_", 1:n))
        }
      }))
      rownames(map) <- paste0("chr", map[["chr"]], "_", map[["pos"]], suffix)
      warning("map contains no marker names. Names constructed from ",
              "chromosome and position.\n", call. = FALSE)
    } 
  } else if (!is.null(gData$map)) {
    ## No map input, but available from gData object. Set map to
    ## map from gData for use later on.
    map <- gData$map
  } else {
    map <- NULL
  }
  ## Modify pheno.
  if (!is.null(pheno)) {
    ## Data.frame or list of data.frames allowed.
    if (!is.data.frame(pheno) &&
        !(is.list(pheno) && all(sapply(X = pheno, FUN = is.data.frame)))) {
      stop("pheno should be a data.frame or a list of data.frames.\n")
    }
    if (is.data.frame(pheno)) {
      ## If not a list already put data.frame/matrix in a list.
      pheno <- setNames(list(pheno), deparse(substitute(pheno)))
    } else {
      if (is.null(names(pheno))) {
        ## Pheno is unnamed list.
        ## Add default names for ease of use in functions.
        names(pheno) <- sapply(X = seq_along(pheno), FUN = function(x) {
          paste0("Trial", x)
        })
        warning("pheno contains no trial names. Default names added.\n",
                call. = FALSE)
      } else {
        ## Pheno has at least some names.
        ## Check for possible empty/missing names in list and add
        ## default names for those. Keep other names.
        if (!all(nzchar(names(pheno)))) {
          names(pheno) <- sapply(X = seq_along(pheno), FUN = function(i) {
            if (!nzchar(names(pheno)[i])) {
              paste0("Trial", i)
            } else {
              names(pheno)[i]
            }
          })
          warning("Some data.frames in pheno contain no trial names. ",
                  "Default names added.\n", call. = FALSE)
        }
      }
    }
    ## Check that first column is always named genotype.
    if (!all(sapply(X = pheno, FUN = function(x) {
      colnames(x)[1] == "genotype"
    }))) {
      stop("First column in pheno should be genotype.\n")
    }
    ## Convert genotype to character.
    for (i in seq_along(pheno)) {
      pheno[[i]]$genotype <- as.character(pheno[[i]]$genotype)
    }
    ## Check that all non-genotype columns in pheno are numerical.
    if (!all(sapply(X = pheno, FUN = function(x) {
      all(sapply(X = x[-1], FUN = is.numeric))
    }))) {
      stop("All trait columns in pheno should be numerical.\n")
    }
    if (!is.null(gData$pheno)) {
      ## gData already contained a pheno object. Overwrite with a warning.
      warning("Existing pheno will be overwritten.\n", call. = FALSE)
    }
  } else if (!is.null(gData$pheno)) {
    ## No pheno input, but available from gData object. Set pheno to
    ## pheno from gData for use later on.
    pheno <- gData$pheno
  } else {
    pheno <- NULL
  }
  ## Modify geno.
  if (!is.null(geno)) {
    ## Either a 2D matrix/data.frame or a 3D array of probabilities is allowed.
    if (!is.data.frame(geno) && !inherits(geno, "Matrix") &&
        !is.matrix(geno)) {
      stop("geno should be a matrix or a data.frame.\n")
    }
    if (!is.matrix(geno)) {
      ## Convert geno to matrix.
      ## This should work for matrices with only numerical entries and
      ## matrices with (some) non numeric entries.
      markers <- as.matrix(geno)
    } else {
      markers <- geno
    }
    ## Check for row names in markers. 
    ## If not available:
    ## 1 - take them from pheno.
    ## 2 - use default names.
    if (all(rownames(markers) == as.character(1:nrow(markers)))) {
      if (is.null(pheno)) {
        ## Default names are constructed as g001, g002, etc. with the number
        ## of zeros dependent on the number of rows.
        rownames(markers) <-
          paste0("g", formatC(1:nrow(markers),
                              width = ceiling(log10(nrow(markers))),
                              flag = "0"))
        warning("geno contains no genotype names. Default names used.\n",
                call. = FALSE)
      } else {
        ## Phenotypic data available. Try to copy names of genotypes from
        ## genotypic data. If dimensions don't match throw an error.
        if (nrow(pheno[[1]]) == nrow(markers)) {
          rownames(markers) <- rownames(pheno[[1]])
          warning("geno contains no genotype names. Names taken from pheno.\n",
                  call. = FALSE)
        } else {
          stop("geno contains no genotype names. Dimensions between ",
               "geno and pheno differ.\n")
        }
      }
    } else {
      ## Sort alphabetically by genotypes.
      markers <- markers[order(rownames(markers)), , drop = FALSE]
    }
    if (is.null(colnames(markers))) {
      ## Check for column names in markers. If not available take them from map.
      ## If map not available or dimensions don't match throw an error.
      if (is.null(map)) {
        stop("geno contains no marker names. Map not available.\n")
      }
      if (nrow(map) != ncol(markers)) {
        stop("geno contains no marker names. Dimensions between geno ",
             "and map differ.\n")
      }
      colnames(markers) <- rownames(map)
      warning("geno contains no marker names. Names taken from map.\n",
              call. = FALSE)
    } else if (!is.null(map)) {
      ## Both markers and map available. Make sure markernames in markers match
      ## markernames in map. Remove non-matching markers from markers.
      ## Map may still contain markers that are not in markers.
      if (any(!colnames(markers) %in%
              rownames(map)[rownames(map) %in% colnames(markers)])) {
        warning("Not all markers in geno are in map. Extra markers ",
                "will be removed.\n", call. = FALSE)
      }
      ## This not only removes markers that are not in map but orders them in
      ## the same order as in map as well.
      ## Distinguish between matrix and array because of number of dims.
      markers <- markers[, colnames(markers) %in% rownames(map), drop = FALSE]
    }
    if (!is.null(gData$markers)) {
      ## gData already contained a markers object. Overwrite with a warning.
      warning("Existing geno will be overwritten.\n", call. = FALSE)
    }
  } else if (!is.null(gData$markers)) {
    ## No marker input, but available from gData object. Set markers to
    ## markers from gData for use later on.
    markers <- gData$markers
  } else {
    markers <- NULL
  }
  ## Modify kin.
  if (!is.null(kin)) {
    ## matrix or list of matrices are allowed.
    if (!inherits(kin, "Matrix") && !is.matrix(kin) &&
        !(is.list(kin) && all(sapply(X = kin, FUN = function(k) {
          inherits(k, "Matrix") || is.matrix(k)
        }))
        )
    ) {
      stop("kin should be a matrix or a list of matrices.\n")
    }
    ## If kin is a list of matrices the number of matrices in the list should
    ## match the number of chromosomes in map.
    if (!is.null(map) && is.list(kin) &&
        length(kin) != length(unique(map[["chr"]]))) {
      stop("kin should be the same length as the number of ",
           "chromosomes in map.\n")
    }
    ## If kin is a named list of matrices names of list items should match names
    ## of chromosomes in map.
    if (!is.null(map) && is.list(kin) &&
        !is.null(names(kin)) && all(names(kin) != unique(map[["chr"]]))) {
      stop("Names of kin should correspond to names of chromosomes in map.\n")
    }
    ## If kin is an unnamed list of matrices add default names.
    if (is.list(kin) && is.null(names(kin))) {
      warning("kin contains no names. Default names added.\n", call. = FALSE)
      names(kin) <- unique(map[["chr"]])
    }
    ## Row and colnames should be provided.
    if (!is.list(kin)) {
      if (is.null(rownames(kin)) || is.null(colnames(kin))) {
        stop("Row and column names in kin cannot be NULL.\n")
      }
    } else {
      if (any(sapply(X = kin, FUN = function(k) {
        is.null(rownames(k)) || is.null(colnames(k))
      }))) {
        stop("Row and column names in kin cannot be NULL.\n")
      }
    }
    ## Genotypes in kin should all be in markers.
    if ((!is.null(markers) && is.list(kin) &&
         any(sapply(X = kin, FUN = function(k) {
           !all(rownames(k) %in% rownames(markers)) ||
             !all(colnames(k) %in% rownames(markers))
         }))) ||
        (!is.null(markers) && !is.list(kin) &&
         (!all(rownames(kin) %in% rownames(markers)) ||
          !all(colnames(kin) %in% rownames(markers))))) {
      stop("Row and column names of kin should be in row names of geno.\n")
    }
    ## Order as in geno and convert to matrix.
    ## match is needed since markers may contain more genotypes than kin.
    ## If markers is NULL only the conversion is done.
    if (is.list(kin)) {
      kin <- lapply(X = kin,
                    FUN = function(k) {
                      as.matrix(k[order(match(rownames(k), rownames(markers))),
                                  order(match(colnames(k), rownames(markers)))])
                    })
    } else {
      if (!is.matrix(kin)) {
        kin <- as.matrix(kin)
      }
      kin <- kin[order(match(rownames(kin), rownames(markers))),
                 order(match(colnames(kin), rownames(markers)))]
    }
    if (!is.null(gData$kinship)) {
      ## gData already contained a kinship object. Overwrite with a warning.
      warning("Existing kinship will be overwritten.\n", call. = FALSE)
    }
  } else if (!is.null(gData$kinship)) {
    ## No kin input, but available from gData object. Set kin to
    ## kinship from gData for use later on.
    kin <- gData$kinship
  } else {
    kin <- NULL
  }
  ## Modify covar.
  if (!is.null(covar)) {
    ## Only a data.frame is allowed.
    if (!inherits(covar, "data.frame")) {
      stop("covar should be a data.frame.\n")
    }
    ## All columns should be numerical, character of factors.
    if (!all(sapply(X = covar, FUN = function(x) {
      is.numeric(x) || is.character(x) || is.factor(x)}))) {
      stop("All columns in covar should be numeric, character or ",
           "factor columns.\n")
    }
    ## Convert character columns to factors.
    covar[sapply(covar, is.character)] <-
      lapply(X = covar[sapply(X = covar, FUN = is.character)],
             FUN = as.factor)
    if (!is.null(gData$covar)) {
      ## gData already contained a covar object. Overwrite with a warning.
      warning("Existing covar will be overwritten.\n", call. = FALSE)
    }
  } else if (!is.null(gData$covar)) {
    ## No covar input, but available from gData object. Set covar to
    ## covar from gData for use later on.
    covar <- gData$covar
  } else {
    covar <- NULL
  }
  ## Create gData object.
  gData <- structure(list(map = map,
                          markers = markers,
                          pheno = pheno,
                          kinship = kin,
                          covar = covar),
                     class = "gData")
  return(gData)
}

#' Summary function for the class \code{gData}
#'
#' Gives a summary for an object of S3 class \code{gData}.
#'
#' @param object An object of class \code{gData}.
#' @param ... Not used.
#' @param trials A vector of trials to include in the summary. These can
#' be either numeric indices or character names of list items in \code{pheno}.
#' If \code{NULL}, all trials are included.
#' 
#' @return A list with a most four components:
#' \describe{
#' \item{mapSum}{A list with number of markers and number of chromosomes in 
#' the map.}
#' \item{markerSum}{A list with number of markers, number of genotypes and 
#' the distribution of the values within the markers.}
#' \item{phenoSum}{A list of data.frames, one per trial with a summary of all
#' traits within the trial.}
#' \item{covarSum}{A list of data.frames, one per trial with a summary of all 
#' covariates within the trial.}
#' }
#' All components are only present in the output if the corresponding content is
#' present in the gData object.
#'
#' @export
summary.gData <- function(object, 
                          ...,
                          trials = NULL) {
  chkTrials(trials, object)
  ## If trials is null set trials to all trials in pheno.
  if (is.null(trials)) {
    trials <- seq_along(object$pheno)
  }
  map <- object$map
  markers <- object$markers
  pheno <-  object$pheno
  covar <- object$covar
  totSum <- vector(mode = "list")
  if (!is.null(map)) {
    mapSum <- list(nMarkers = nrow(map), nChr = length(unique(map[["chr"]])))
    totSum$mapSum <- mapSum
  }
  if (!is.null(markers)) {
    markerSum <- list(nMarkers = ncol(markers), nGeno = nrow(markers),
                      markerContent = round(prop.table(table(as.vector(markers), 
                                                             useNA = "always")), 
                                            2))
    totSum$markerSum <- markerSum
  }
  if (!is.null(pheno)) {
    phenoSum <- sapply(X = names(pheno[trials]), FUN = function(trial) {
      trSum <- do.call(cbind, lapply(X = pheno[[trial]][, -1], FUN = summaryNA))
      attr(x = trSum, which = "nGeno") <- 
        length(unique(pheno[[trial]][["genotype"]]))
      return(trSum)
    }, simplify = FALSE)
    totSum$phenoSum <- phenoSum
  }
  if (!is.null(covar)) {
    covarSum <- summary(covar)
    totSum$covarSum <- covarSum
  }
  return(structure(totSum, class = "summary.gData"))
}

#' Printing summazed objects of class gData
#'
#' \code{print} method for object of class summary.gData created by summarizing
#' objects of class gData.
#'
#' @param x An object of class \code{summary.gData}.
#' @param ... Not used.
#'
#' @noRd
#' @export
print.summary.gData <- function(x, 
                                ...) {
  if (!is.null(x$mapSum)) {
    cat("map\n")
    cat("\tNumber of markers:", x$mapSum$nMarkers, "\n")
    cat("\tNumber of chromosomes:", x$mapSum$nChr, "\n\n")
  }
  if (!is.null(x$markerSum)) {
    cat("markers\n")
    cat("\tNumber of markers:", x$markerSum$nMarkers, "\n")
    cat("\tNumber of genotypes:", x$markerSum$nGeno, "\n")
    cat("\tContent:\n")
    cat("\t", names(x$markerSum$markerContent), "\n")
    cat("\t", x$markerSum$markerContent, "\n\n")
  }
  if (!is.null(x$phenoSum)) {
    cat("pheno\n")
    cat("\tNumber of trials:", length(x$phenoSum), "\n\n")
    for (i in seq_along(x$phenoSum)) {
      if (!is.null(names(x$phenoSum)[i])) {
        cat("\t", names(x$phenoSum)[i], ":\n", sep = "")
      } else {
        cat("\tTrial ", i, ":\n", sep = "")
      }
      cat("\t\tNumber of traits:", ncol(x$phenoSum[[i]]), "\n")
      cat("\t\tNumber of genotypes:", 
          attr(x = x$phenoSum[[i]], which = "nGeno"), "\n\n")
      attr(x = x$phenoSum[[i]], which = "nGeno") <- NULL
      print(x$phenoSum[[i]])
      cat("\n")
    }
  }
  if (!is.null(x$covarSum)) {
    cat("covar\n")
    cat("\tNumber of covariates:", ncol(x$covarSum), "\n")
    print(x$covarSum)
  }
}
