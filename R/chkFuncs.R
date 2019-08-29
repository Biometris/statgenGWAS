#' @keywords internal
chkGData <- function(gData = NULL,
                     comps = c("map", "markers", "pheno")) {
  if (is.null(gData) || !inherits(gData, "gData")) {
    stop("gData should be a valid gData object", call. = FALSE)
  }
  for (comp in comps) {
    if (is.null(gData[[comp]]))
      stop(paste("gData should contain", comp), call. = FALSE)
  }
}

chkMarkers <- function(markers,
                       dim = 2) {
  if (dim == 2) {
    if (!is.numeric(markers)) {
      stop(paste("markers in gData should be a numerical matrix. Use",
                 "recodeMarkers first for recoding.\n"), call. = FALSE)
    }
  } else if (dim == 3) {
    if (!inherits(markers, "array")) {
      stop("markers should be a three-dimensional array.\n", call. = FALSE)
    }
  }
  if (anyNA(markers)) {
    stop("markers contains missing values. Impute or remove these first.\n",
         call. = FALSE)
  }
}

chkEnvs <- function(envs,
                    gData) {
  if (!is.null(envs) && !is.numeric(envs) && !is.character(envs)) {
    stop("environments should be a numeric or character vector.\n",
         call. = FALSE)
  }
  if ((is.character(envs) && !all(envs %in% names(gData$pheno))) ||
      (is.numeric(envs) && any(envs > length(gData$pheno)))) {
    stop("environments should be in pheno.\n", call. = FALSE)
  }
}

chkTraits <- function(traits,
                      envs,
                      gData,
                      multi) {
  if (!is.null(traits) && !is.numeric(traits) && !is.character(traits)) {
    if (multi) {
      stop("traits should be a numeric or character vector.\n", call. = FALSE)
    } else if (length(traits) > 1) {
      stop("trait should be a single numeric or character.\n", call. = FALSE)
    }
  }
  for (env in envs) {
    if ((is.character(traits) &&
         !all(hasName(x = gData$pheno[[env]], traits))) ||
        (is.numeric(traits) &&
         (any(traits == 1) || any(traits > ncol(gData$pheno[[env]]))))) {
      stop(paste("For", env, "not all traits are columns in pheno.\n"),
           call. = FALSE)
    }
  }
}

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

chkSnpCov <- function(snpCov,
                      gData) {
  if (!is.null(snpCov) && !all(snpCov %in% colnames(gData$markers))) {
    stop("All snpCov should be in markers.\n", call. = FALSE)
  }
}

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
    stop(paste("kin should be a list of matrices of length equal to the",
               "number of chromosomes in the map.\n", call. = FALSE))
  }
}

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
    stop(paste0(match.call()$x, " should be a single numerical value",
                txt, ".\n"), call. = FALSE)
  }
}


