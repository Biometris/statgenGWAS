## Add covariates and snpCovariates to phenotypic data and convert covariate
## factors to dummy varables.
#' @keywords internal
expandPheno <- function(gData,
                        trial,
                        covar,
                        snpCov,
                        ref = NULL) {
  ## Add covariates to pheno data.
  if (is.null(covar)) {
    phTr <- gData$pheno[[trial]]
    covTr <- NULL
  } else {
    ## Append covariates to pheno data. Merge to remove values from pheno that
    ## are missing in covar.
    phTr <- merge(gData$pheno[[trial]], gData$covar[covar],
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
        model.matrix(object = covFormula, data = droplevels(phTr))))[, -1]
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
    ## Distinguish between 2- and 3-dimensional marker data.
    if (length(dim(gData$markers)) == 2) {
      ## Add snp covariates to covar.
      covTr <- c(covTr, snpCov)
      ## Add snp covariates to pheno data.
      phTr <- merge(phTr, gData$markers[, snpCov, drop = FALSE],
                    by.x = "genotype", by.y = "row.names")
      colnames(phTr)[(ncol(phTr) - length(snpCov) + 1):ncol(phTr)] <- snpCov
    } else if (length(dim(gData$markers)) == 3) {
      for (snpCovar in snpCov) {
        ## Get alleles for current covariate and remove reference allele.
        allCov <- gData$markers[, snpCovar , -ref]
        ## Remove alleles with only zeros.
        allCov <- allCov[, apply(X = allCov, MARGIN = 2, FUN = function(a) {
          any(a > 0)
        })]
        ## Rename columns to combination of allele and marker.
        colnames(allCov) <- paste0(snpCovar, "_", colnames(allCov))
        ## Add snp covariates to covar.
        covTr <- c(covTr, colnames(allCov))
        ## Add snp covariates to pheno data.
        phTr <- merge(phTr, allCov, by.x = "genotype", by.y = "row.names")
      }
    }
  }
  return(list(phTr = phTr, covTr = covTr))
}

## Helper function for computing (or extracting kinship matrices)
## 1 - If kin is supplied use kin
## 2 - Get kin from gData object
## 3 - Compute kin from markers (and map for GLSMethod multi)
#' @keywords internal
computeKin <- function(GLSMethod,
                       kin = NULL,
                       gData = NULL,
                       markers = NULL,
                       map = NULL,
                       kinshipMethod = NULL) {
  if (GLSMethod == "single") {
    if (!is.null(kin)) {
      ## kin is supplied as input. Convert to matrix.
      K <- as.matrix(kin)
    } else if (!is.null(gData$kinship) && !inherits(gData$kinship, "list")) {
      ## Get kin from gData object.
      K <- gData$kinship
    } else {
      ## Compute K from markers.
      K <- kinship(X = markers, method = kinshipMethod)
    }
    K <- K[order(match(rownames(K), rownames(markers))),
           order(match(colnames(K), rownames(markers)))]
  } else if (GLSMethod == "multi") {
    if (!is.null(kin)) {
      ## kin is supplied as input. Convert to matrices.
      K <- lapply(X = kin, FUN = as.matrix)
    } else if (!is.null(gData$kinship) && inherits(gData$kinship, "list")) {
      ## Get kin from gData object.
      K <- gData$kinship
    } else {
      ## Compute chromosome specific kinship matrices.
      K <- chrSpecKin(gData = createGData(geno = markers, map = map),
                      kinshipMethod = kinshipMethod)
    }
    K <- lapply(X = K, FUN = function(k) {
      k[order(match(rownames(k), rownames(markers))),
        order(match(colnames(k), rownames(markers)))]
    })
  }
  return(K)
}

## Compute chromosome specific kinship matrices.
chrSpecKin <- function(gData,
                       kinshipMethod) {
  chrs <- unique(gData$map$chr[rownames(gData$map) %in%
                                 colnames(gData$markers)])
  if (length(chrs) == 1) {
    stop(paste("Chromosome specific kinship calculation not possible since",
               "map contains only 1 chromosome.\n"))
  }
  ## Create list of zero matrices.
  KChr <- setNames(
    replicate(n = length(chrs),
              matrix(data = 0, nrow = nrow(gData$markers),
                     ncol = nrow(gData$markers),
                     dimnames = list(rownames(gData$markers),
                                     rownames(gData$markers))),
              simplify = FALSE),
    paste0("KChr", chrs))
  ## Create vector of marker numbers per chromosome.
  denom <- setNames(rep(x = 0, times = length(chrs)), chrs)
  for (chr in chrs) {
    ## Extract markers for current chromosome.
    chrMrk <- which(colnames(gData$markers) %in%
                      rownames(gData$map[gData$map$chr == chr, ]))
    ## Compute kinship for current chromosome only. Denominator = 1, division
    ## is done later.
    if (length(dim(gData$markers)) == 2) {
      K <- kinship(X = gData$markers[, chrMrk, drop = FALSE],
                   method = kinshipMethod, denominator = 1)
      ## Compute number of markers for other chromosomes.
      denom[which(chrs == chr)] <-
        ncol(gData$markers[, -chrMrk, drop = FALSE])
    } else if (length(dim(gData$markers)) == 3) {
      K <- kinship(X = gData$markers[, chrMrk, , drop = FALSE],
                   method = kinshipMethod, denominator = 1)
      ## Compute chromosome length.
      ## Add extra bits for first and last marker as in kinship calculation.
      pos <- gData$map[gData$map$chr == chr, "pos"]
      chrLen <- max(pos) - min(pos) +
        (pos[2] - pos[1] + rev(pos)[1] - rev(pos)[2]) / 2
    }
    for (i in setdiff(1:length(chrs), which(chr == chrs))) {
      ## Add computed kinship to all other matrices in KChr.
      KChr[[i]] <- KChr[[i]] + K
      if (length(dim(gData$markers)) == 3) {
        ## Add chromosome length to all other denominators in denom.
        denom[i] <- denom[i] + chrLen
      }
    }
  }
  ## Divide matrix for current chromosome by number of markers in other
  ## chromosomes.
  for (i in 1:length(KChr)) {
    KChr[[i]] <- KChr[[i]] / denom[i]
  }
  return(KChr)
}

#' Row bind data.frames
#'
#' Helper function for row binding data.frames with diffent columns.
#'
#' @param dfList A list of data.frames.
#'
#' @keywords internal
dfBind <- function(dfList) {
  ## Filter empty data.frames from dfList
  dfList <- Filter(f = function(x) nrow(x) > 0, x = dfList)
  if (length(dfList) == 0) {
    return(data.frame())
  }
  ## Get variable names from all data.frames.
  allNms <- unique(unlist(lapply(dfList, names)))
  ## rbind all data.frames setting values for missing columns to NA.
  do.call(rbind,
          c(lapply(X = dfList, FUN = function(x) {
            nwDat <- sapply(X = setdiff(allNms, names(x)), FUN = function(y) {
              NA
            })
            data.frame(c(x, nwDat), stringsAsFactors = FALSE)
          }), make.row.names = FALSE)
  )
}

