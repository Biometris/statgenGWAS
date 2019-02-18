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
#' A three dimensional array of probabilities may be provided as well with
#' genotypes in the first, markers in the second and alleles in the third
#' dimension.\cr
#' If no row names are provided, they are taken from \code{pheno} (if supplied and
#' dimension matches). If no column names are provided, the row names
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
#' environments.
#' @param covar A data.frame with extra covariates per genotype. Genotypes
#' should be in the rows.
#'
#' @return An object of class \code{gData} with the following components:
#' \item{\code{map}}{a data.frame containing map data. Map is sorted by
#' chromosome and position.}
#' \item{\code{markers}}{a sparse matrix from the Matrix package containing
#' marker information in case of numerical genotypic data, a standard matrix
#' otherwise.\cr
#' If \code{geno} is a three dimensional array, \code{markers} is a three dimensional
#' array as well.}
#' \item{\code{pheno}}{a list of data.frames containing phenotypic data.}
#' \item{\code{kinship}}{a kinship matrix of class \code{dsyMatrix} from the
#'  Matrix package.}
#' \item{\code{covar}}{a data.frame with extra covariates.}
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{summary.gData}}
#'
#' @examples set.seed(1234)
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
#' @importFrom methods as
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
    stop(paste("At least one of geno, map, kin, pheno and covar should",
               "be provided.\n"))
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
      warning("existing map will be overwritten.\n", call. = FALSE)
    }
    ## Extract columns and order.
    map <- map[order(map$chr, map$pos), c("chr", "pos")]
    if (all(rownames(map) == as.character(1:nrow(map)))) {
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
      rownames(map) <- paste0("chr", map$chr, "_", map$pos, suffix)
      warning(paste("map contains no marker names. Names constructed from",
                    "chromosome and position.\n"),
              call. = FALSE)
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
      stop("pheno should be a data.frame or a list data.frames.\n")
    }
    if (is.data.frame(pheno)) {
      ## If not a list already put data.frame/matrix in a list.
      pheno <- setNames(list(pheno), deparse(substitute(pheno)))
    } else {
      if (is.null(names(pheno))) {
        ## Pheno is unnamed list.
        ## Add default names for ease of use in functions.
        names(pheno) <- sapply(X = seq_along(pheno), FUN = function(x) {
          paste0("Environment", x)
        })
        warning("pheno contains no environment names. Default names added.\n",
                call. = FALSE)
      } else {
        ## Pheno has at least some names.
        ## Check for possible empty/missing names in list and add
        ## default names for those. Keep other names.
        if (!all(sapply(X = names(pheno), FUN = nzchar))) {
          names(pheno) <- sapply(X = seq_along(pheno), FUN = function(x) {
            if (!nzchar(names(pheno)[x])) {
              paste0("Environment", x)
            } else {
              names(pheno)[x]
            }
          })
          warning(paste("Some data.frames in pheno contain no environment",
                        "names. Default names added.\n", call. = FALSE))
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
    ## Check that all non-genotype columns are numerical.
    if (!all(sapply(X = pheno, FUN = function(x) {
      all(sapply(X = x[-1], FUN = is.numeric))
    }))) {
      stop("all trait columns in pheno should be numerical.\n")
    }
    if (!is.null(gData$pheno)) {
      ## gData already contained a pheno object. Overwrite with a warning.
      warning("existing pheno will be overwritten.\n", call. = FALSE)
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
        !is.matrix(geno) && !(is.array(geno) && length(dim(geno)) == 3)) {
      stop("geno should be a matrix, data.frame or an array.\n")
    }
    isMat <- length(dim(geno)) == 2
    isPMat <- length(dim(geno)) == 3
    if (isMat) {
      if (is.data.frame(geno) || is.matrix(geno)) {
        if (is.numeric(unlist(geno))) {
          ## Convert geno to Matrix of class Matrix.
          markers <- as(geno, "Matrix")
        } else {
          ## Markers contain non numeric entries. Matrix class cannot handle
          ## this but it is needed for use in recodeMarkers.
          markers <- as.matrix(geno)
        }
      } else {
        markers <- geno
      }
    } else if (isPMat) {
      markers <- geno
    }
    ## Check for row names in markers. If not available take them from pheno or
    ## use default names.
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
          stop(paste("geno contains no genotype names. Dimensions between",
                     "geno and pheno differ.\n"))
        }
      }
    } else {
      ## Sort alphabetically by genotypes.
      ## Distinguish between matrix and array because of number of dims.
      if (isMat) {
        markers <- markers[order(rownames(markers)), , drop = FALSE]
      } else if (isPMat) {
        markers <- markers[order(rownames(markers)), , , drop = FALSE]
      }
    }
    if (is.null(colnames(markers))) {
      ## Check for column names in markers. If not available take them from map.
      ## If map not available or dimensions don't match throw an error.
      if (is.null(map)) {
        stop("geno contains no marker names. Map not available.\n")
      }
      if (nrow(map) != ncol(markers)) {
        stop(paste("geno contains no marker names. Dimensions between geno",
                   "and map differ.\n"))
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
        warning(paste("not all markers in geno are in map. Extra markers",
                      "will be removed.\n"), call. = FALSE)
      }
      if (isMat) {
        ## This not only removes markers that are not in map but orders them in
        ## the same order as in map as well.
        ## Distinguish between matrix and array because of number of dims.
        markers <- markers[, colnames(markers) %in% rownames(map), drop = FALSE]
      } else if (isPMat) {
        markers <- markers[, colnames(markers) %in% rownames(map), ,
                           drop = FALSE]
        ## Check that probabilities for each marker sum to one.
        ## Always rescale values.
        ## Throw a warning if difference from one too large.
        genoMrk <- setNames(as.data.frame(matrix(nrow = 0, ncol = 2)),
                            c("geno", "marker"))
        for (mrk in colnames(markers)) {
          mrkProbs <- rowSums(markers[, mrk, ], na.rm = TRUE)
          if (any(abs(mrkProbs - 1) > 1e-2)) {
            genoMrk <-
              rbind(genoMrk,
                    data.frame(geno =
                                 names(mrkProbs[abs(mrkProbs - 1) > 1e-2]),
                               marker = mrk, stringsAsFactors = FALSE))
          }
          markers[, mrk, ] <- markers[, mrk, ] / mrkProbs
        }
        if (nrow(genoMrk) > 0) {
          genoMrk$genoMarker <- paste(genoMrk$geno, "\t", genoMrk$marker)
          warning(paste0("Probabilities differ from 1 for the following ",
                         "combinations of genotype and markers:\n",
                         paste(genoMrk$genoMarker, collapse = "\n")),
                  call. = FALSE)
        }
      }
    }
    if (!is.null(gData$markers)) {
      ## gData already contained a markers object. Overwrite with a warning.
      warning("existing geno will be overwritten.\n", call. = FALSE)
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
        length(kin) != length(unique(map$chr))) {
      stop(paste("kin should be the same length as the number of",
                 "chromosomes in map.\n"))
    }
    ## If kin is a named list of matrices names of list items should match names
    ## of chromosomes in map.
    if (!is.null(map) && is.list(kin) &&
        !is.null(names(kin)) && names(kin) != unique(map$chr)) {
      stop("names of kin should correspond to names of chromosomes in map.\n")
    }
    ## If kin is an unnamed list of matrices add default names.
    if (is.list(kin) && is.null(names(kin))) {
      warning("kin contains no names. Default names added.\n", call. = FALSE)
      names(kin) <- unique(map$chr)
    }
    ## Row and colnames should be provided.
    if (!is.list(kin)) {
      if (is.null(rownames(kin)) || is.null(colnames(kin))) {
        stop("row and column names in kin cannot be NULL.\n")
      }
    } else {
      if (any(sapply(X = kin, FUN = function(k) {
        is.null(rownames(k)) || is.null(colnames(k))
      }))) {
        stop("row and column names in kin cannot be NULL.\n")
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
      stop(paste("row and column names of kin should be in row",
                 "names of geno.\n"))
    }
    ## Order as in geno and convert to symmetric matrix in Matrix class.
    ## match is needed since markers may contain more genotypes than kin.
    ## If markers is NULL only the conversion is done.
    if (is.list(kin)) {
      kin <- lapply(X = kin,
                    FUN = function(k) {
                      as(k[order(match(rownames(k), rownames(markers))),
                           order(match(colnames(k), rownames(markers)))],
                         "dsyMatrix")
                    })
    } else {
      if (is.matrix(kin)) {
        kin <- as(kin, "dsyMatrix")
      }
      kin <- kin[order(match(rownames(kin), rownames(markers))),
                 order(match(colnames(kin), rownames(markers)))]
    }
    if (!is.null(gData$kinship)) {
      ## gData already contained a kinship object. Overwrite with a warning.
      warning("existing kinship will be overwritten.\n", call. = FALSE)
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
    if (!is.data.frame(covar)) {
      stop("covar should be a data.frame.\n")
    }
    ## All columns should be numerical, character of factors.
    if (!all(sapply(X = covar, FUN = function(x) {
      is.numeric(x) || is.character(x) || is.factor(x)}))) {
      stop(paste("all columns in covar should be numeric, character or",
                 "factor columns.\n"))
    }
    ## Convert character columns to factors.
    covar[sapply(covar, is.character)] <-
      lapply(X = covar[sapply(X = covar, FUN = is.character)],
             FUN = as.factor)
    if (!is.null(gData$covar)) {
      ## gData already contained a covar object. Overwrite with a warning.
      warning("existing covar will be overwritten.\n", call. = FALSE)
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

#' Summary function for the class \code{GData}
#'
#' Gives a summary for an object of S3 class \code{GData}.
#'
#' @param object An object of class \code{GData}
#' @param ... Not used
#'
#' @export
summary.gData <- function(object, ...) {
  map <- object$map
  markers <- object$markers
  pheno <-  object$pheno
  kinship <- object$kinship
  covar <- object$covar
  if (!is.null(map)) {
    cat("map\n")
    cat("\tNumber of markers:", nrow(map), "\n")
    cat("\tNumber of chromosomes:", length(unique(map$chr)), "\n\n")
  }
  if (!is.null(markers)) {
    cat("markers\n")
    cat("\tNumber of markers:", ncol(markers), "\n")
    cat("\tNumber of genotypes:", nrow(markers), "\n")
    cat("\tContent:\n")
    tab <- (round(prop.table(table(as.vector(markers), useNA = "always")), 2))
    cat("\t", names(tab), "\n")
    cat("\t", tab, "\n\n")
  }
  if (!is.null(pheno)) {
    cat("pheno\n")
    cat("\tNumber of environments:", length(pheno), "\n\n")
    for (i in 1:length(pheno)) {
      if (!is.null(names(pheno)[i])) {
        cat("\t", names(pheno)[i], ":\n", sep = "")
      } else {
        cat("\tEnvironment ", i, ":\n", sep = "")
      }
      cat("\t\tNumber of traits:", ncol(pheno[[i]]) - 1, "\n")
      cat("\t\tNumber of genotypes:", length(unique(pheno[[i]]$genotype)),
          "\n\n")
      print(summary(pheno[[i]][, -1]))
      cat("\n")
    }
  }
  if (!is.null(covar)) {
    cat("covar\n")
    cat("\tNumber of covariates:", ncol(covar), "\n")
    print(summary(covar))
  }
}







