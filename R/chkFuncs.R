#' @noRd
#' @keywords internal
chkGData <- function(gData = NULL,
                     comps = c("map", "markers", "pheno")) {
  if (is.null(gData) || !inherits(gData, "gData")) {
    stop("gData should be a valid gData object", call. = FALSE)
  }
  for (comp in comps) {
    if (is.null(gData[[comp]]))
      stop("gData should contain ", comp, call. = FALSE)
  }
}

#' @noRd
#' @keywords internal
chkTrials <- function(trials,
                      gData) {
  if (!is.null(trials) && !is.numeric(trials) && !is.character(trials)) {
    stop("trials should be a numeric or character vector.\n",
         call. = FALSE)
  }
  if ((is.character(trials) && !all(trials %in% names(gData$pheno))) ||
      (is.numeric(trials) && any(trials > length(gData$pheno)))) {
    stop("trials should be in pheno.\n", call. = FALSE)
  }
}

#' @noRd
#' @keywords internal
chkTraits <- function(traits,
                      trials,
                      gData,
                      multi) {
  if (!is.null(traits) && !is.numeric(traits) && !is.character(traits)) {
    stop("traits should be a numeric or character vector.\n", call. = FALSE)
  } 
  if (!multi && length(traits) > 1) {
    stop("traits should be a single numeric or character value.\n", 
         call. = FALSE)
  }
  for (trial in trials) {
    if ((is.character(traits) &&
         !all(hasName(x = gData$pheno[[trial]], traits))) ||
        (is.numeric(traits) &&
         (any(traits == 1) || any(traits > ncol(gData$pheno[[trial]]))))) {
      stop("For ", trial, " not all traits are columns in pheno.\n",
           call. = FALSE)
    }
  }
}

#' @noRd
#' @keywords internal
chkNum <- function(x,
                   min = NULL,
                   max = NULL) {
  if (missing(x) || length(x) > 1 || !is.numeric(x) || isTRUE(x < min) ||
      isTRUE(x > max)) {
    if (!is.null(min) && !is.null(max)) {
      txt <- paste(" between", min, "and", max)
    } else if (!is.null(min)) {
      txt <- paste(" greater than", min)
    } else if (!is.null(max)) {
      txt <- paste(" smaller than", max)
    } else {
      txt <- ""
    }
    stop(match.call()$x, " should be a single numerical value", txt, ".\n", 
         call. = FALSE)
  }
}

#' @noRd
#' @keywords internal
chkMarkers <- function(markers,
                       dim = 2) {
  if (dim == 2) {
    if (!is.numeric(markers)) {
      stop("markers in gData should be a numerical matrix. Use ",
           "codeMarkers first for recoding.\n", call. = FALSE)
    }
  } else if (dim == 3) {
    if (!inherits(markers, "array") || length(dim(markers)) != 3) {
      stop("markers should be a three-dimensional array.\n", call. = FALSE)
    }
  }
  if (anyNA(markers)) {
    stop("markers contains missing values. Impute or remove these first.\n",
         call. = FALSE)
  }
}

#' @noRd
#' @keywords internal
chkCovar <- function(covar,
                     gData) {
  if (!is.null(covar) && !is.numeric(covar) && !is.character(covar)) {
    stop("covar should be a numeric or character vector.\n", call. = FALSE)
  }
  if ((is.character(covar) && !all(hasName(x = gData$covar, name = covar))) ||
      (is.numeric(covar) && any(covar > ncol(gData$covar)))) {
    stop("covar should be columns in covar in gData.\n", call. = FALSE)
  }
}

#' @noRd
#' @keywords internal
chkSnpCov <- function(snpCov,
                      gData) {
  if (!is.null(snpCov) && !all(snpCov %in% colnames(gData$markers))) {
    stop("All snpCov should be in markers.\n", call. = FALSE)
  }
}

#' @noRd
#' @keywords internal
chkKin <- function(kin,
                   gData,
                   GLSMethod) {
  if (GLSMethod == "single" && !is.null(kin) &&
      !(inherits(kin, "Matrix") || is.matrix(kin))) {
    stop("kin should be a matrix.\n", call. = FALSE)
  }
  if (GLSMethod == "multi" && !is.null(kin) &&
      (!is.list(kin) || !all(sapply(kin, FUN = function(k) {
        is.matrix(k) || inherits(k, "Matrix")})) ||
       length(kin) != length(unique(gData$map$chr)))) {
    stop("kin should be a list of matrices of length equal to the ",
         "number of chromosomes in the map.\n", call. = FALSE)
  }
}




