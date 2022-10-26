#' Add covariates and snpCovariates to phenotypic data and convert covariate
#' factors to dummy varables.
#' 
#' @noRd
#' @keywords internal
expandPheno <- function(gData,
                        trial,
                        covar,
                        snpCov,
                        ref = NULL) {
  phTr <- gData$pheno[[trial]]
  ## Add covariates to pheno data.
  if (is.null(covar)) {
    covTr <- NULL
  } else {
    ## Remove columns in covar from pheno.
    phTr <- phTr[, !colnames(phTr) %in% covar, drop = FALSE]
    ## Append covariates to pheno data. Merge to remove values from pheno that
    ## are missing in covar.
    phTr <- merge(phTr, as.data.frame(gData$covar)[covar],
                  by.x = "genotype", by.y = "row.names")
    ## Remove rows from phTr with missing covar check if there are
    ## missing values.
    phTr <- phTr[complete.cases(phTr[covar]), ]
    ## Expand covariates that are a factor (i.e. dummy variables are created)
    ## using model.matrix. The new dummies are attached to phTr, and covar
    ## is changed accordingly
    factorCovs <- which(sapply(X = gData$covar[covar], FUN = is.factor))
    if (length(factorCovs) > 0) {
      ## Create dummy variables without intercept.
      covFormula <- as.formula(paste("genotype ~ ",
                                     paste(covar[factorCovs], collapse = "+")))
      extraCov <- as.data.frame(suppressWarnings(
        model.matrix.lm(object = covFormula, data = droplevels(phTr))))[, -1]
      ## Add dummy variables to pheno data.
      phTr <- cbind(phTr[, -which(colnames(phTr) %in% names(factorCovs))],
                    extraCov)
      ## Modify covar to suit newly defined columns
      covTr <- c(covar[-factorCovs], colnames(extraCov))
    } else {
      covTr <- covar
    }
  }
  if (!is.null(snpCov)) {
    ## Add snp covariates to covar.
    covTr <- c(covTr, snpCov)
    ## Add snp covariates to pheno data.
    phTr <- merge(phTr, gData$markers[, snpCov, drop = FALSE],
                  by.x = "genotype", by.y = "row.names")
    colnames(phTr)[(ncol(phTr) - length(snpCov) + 1):ncol(phTr)] <- snpCov
  }
  return(list(phTr = phTr, covTr = covTr))
}

#' Helper function for computing (or extracting kinship matrices)
#' 1 - If kin is supplied use kin
#' 2 - Get kin from gData object
#' 3 - Compute kin from markers (and map for GLSMethod multi)
#' 
#' @noRd
#' @keywords internal
computeKin <- function(GLSMethod,
                       kin = NULL,
                       gData = NULL,
                       markers = NULL,
                       map = NULL,
                       kinshipMethod = NULL,
                       MAF = NULL) {
  if (GLSMethod == "single") {
    if (!is.null(kin)) {
      ## kin is supplied as input. Convert to matrix.
      K <- as.matrix(kin)
      kinshipMethod <- NULL
    } else if (!is.null(gData$kinship) && !inherits(gData$kinship, "list")) {
      ## Get kin from gData object.
      K <- gData$kinship
      kinshipMethod <- NULL
    } else {
      ## Compute K from markers.
      K <- kinship(X = markers, method = kinshipMethod, MAF = MAF)
    }
    K <- K[order(match(rownames(K), rownames(markers))),
           order(match(colnames(K), rownames(markers)))]
  } else if (GLSMethod == "multi") {
    if (!is.null(kin)) {
      ## kin is supplied as input. Convert to matrices.
      K <- lapply(X = kin, FUN = as.matrix)
      kinshipMethod <- NULL
    } else if (!is.null(gData$kinship) && inherits(gData$kinship, "list")) {
      ## Get kin from gData object.
      K <- gData$kinship
      kinshipMethod <- NULL
    } else {
      ## Compute chromosome specific kinship matrices.
      K <- chrSpecKin(markers = markers, map = map, 
                      kinshipMethod = kinshipMethod, MAF = MAF)
    }
    K <- lapply(X = K, FUN = function(k) {
      k[order(match(rownames(k), rownames(markers))),
        order(match(colnames(k), rownames(markers)))]
    })
  }
  attr(K, which = "method") <- kinshipMethod
  return(K)
}

## Compute chromosome specific kinship matrices.
chrSpecKin <- function(markers, 
                       map,
                       kinshipMethod,
                       MAF) {
  chrs <- unique(map[rownames(map) %in% colnames(markers), "chr"])
  if (length(chrs) == 1) {
    stop("Chromosome specific kinship calculation not possible since ",
         "map contains only 1 chromosome.\n")
  }
  ## Create list of zero matrices.
  KChr <- setNames(replicate(n = length(chrs), 
                             matrix(data = 0, nrow = nrow(markers),
                                    ncol = nrow(markers),
                                    dimnames = list(rownames(markers),
                                                    rownames(markers))),
                             simplify = FALSE), chrs)
  ## Create vector of marker numbers per chromosome.
  denom <- setNames(rep(x = 0, times = length(chrs)), chrs)
  for (chr in chrs) {
    ## Extract markers for current chromosome.
    chrMrk <- which(colnames(markers) %in% rownames(map[map[["chr"]] == chr, ]))
    ## Compute kinship for current chromosome only. Denominator = 1, division
    ## is done later.
    K <- kinship(X = markers[, chrMrk, drop = FALSE], method = kinshipMethod, 
                 MAF = MAF, denominator = 1)
    ## Compute number of markers for other chromosomes.
    denom[which(chrs == chr)] <- ncol(markers[, -chrMrk, drop = FALSE])
    for (i in setdiff(seq_along(chrs), which(chr == chrs))) {
      ## Add computed kinship to all other matrices in KChr.
      KChr[[i]] <- KChr[[i]] + K
    }
  }
  ## Divide matrix for current chromosome by number of markers in other
  ## chromosomes.
  for (i in seq_along(KChr)) {
    KChr[[i]] <- KChr[[i]] / denom[i]
  }
  return(KChr)
}

#' Helper function for creating summaries that always display NA.
#' 
#' @noRd
#' @keywords internal
summaryNA <- function(dat) {
  if (!any(is.na(dat))) {
    return(c(summary(dat), "NA's" = 0))
  } else{
    return(summary(dat))
  }
}


