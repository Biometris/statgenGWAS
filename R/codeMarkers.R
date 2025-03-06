#' Code and impute markers
#'
#' \code{codeMarkers} codes markers in a \code{gData} object and optionally
#' performs imputation of missing values as well.\cr
#' The function performs the following steps:\cr
#' \enumerate{
#' \item{replace strings in \code{naStrings} by \code{NA}.}
#' \item{remove genotypes with a fraction of missing values higher than
#' \code{nMissGeno}.}
#' \item{remove SNPs with a fraction of missing values higher than
#' \code{nMiss}.}
#' \item{recode SNPs to numerical values.}
#' \item{remove SNPs with a minor allele frequency lower than \code{MAF}.}
#' \item{optionally remove duplicate SNPs.}
#' \item{optionally impute missing values.}
#' \item{repeat steps 5. and 6. if missing values are imputed.}
#' }
#'
#' @param gData An object of class \code{gData} containing at least
#' \code{markers}.
#' @param refAll A character string indicating the reference allele used when
#' recoding markers.\cr
#' If "minor", then the recoding is done using the minor allele as reference
#' allele. Alternatively a single character can be supplied as a reference
#' allele for the whole set of SNPs, or a character vector with a reference
#' allele per SNP.
#' @param nMissGeno A numerical value between 0 and 1. Genotypes with a
#' fraction of missing values higher than \code{nMissGeno} will be removed. 
#' Genotypes with only missing values will always be removed.
#' @param nMiss A numerical value between 0 and 1. SNPs with a fraction of
#' missing values higher than \code{nMiss} will be removed. SNPs with only 
#' missing values will always be removed.
#' @param MAF A numerical value between 0 and 1. SNPs with a Minor Allele
#' Frequency (MAF) below this value will be removed. Only one of \code{MAF} 
#' and \code{MAC} may be specified.
#' @param MAC A numerical value. SNPs with Minor Allele Count (MAC) below this 
#' value will be removed. Only one of \code{MAF} and \code{MAC} may be 
#' specified.
#' @param removeDuplicates Should duplicate SNPs be removed?
#' @param keep A vector of SNPs that should never be removed in the whole
#' process.
#' @param impute Should imputation of missing values be done?
#' @param imputeType A character string indicating what kind of imputation of
#' values should be done.\cr
#' \itemize{
#' \item{fixed - missing values will be replaced by a given fixed value.}
#' \item{random - missing values will be replaced by a random value calculated
#' using allele frequencies per SNP.}
#' \item{beagle - missing values will be imputed using beagle software, version
#' 5.2. Beagle only accepts integers as map positions. If you use this option,
#' please cite the original papers in your publication (see references).}
#' }
#' @param fixedValue A numerical value used for replacing missing values in
#' case \code{inputType} is fixed.
#' @param naStrings A character vector of strings to be treated as NA.
#' @param verbose Should a summary of the performed steps be printed?
#'
#' @returns A copy of the input \code{gData} object with markers replaced by
#' coded and imputed markers.
#'
#' @references S R Browning and B L Browning (2007) Rapid and accurate haplotype
#' phasing and missing data inference for whole genome association studies by 
#' use of localized haplotype clustering. Am J Hum Genet 81:1084-1097. 
#' \doi{10.1086/521987}
#'
#' @examples ## Create markers
#' markers <- matrix(c(
#' "AA",   "AB",   "AA",   "BB",   "BA",   "AB",   "AA",   "AA",   NA,  "AA",
#' "AA",   "AA",   "BB",   "BB",   "AA",   "AA",   "BB",   "AA",   NA,  "AA",
#' "AA",   "BA",   "AB",   "BB",   "AB",   "AB",   "AA",   "BB",   NA,  "AA",
#' "AA",   "AA",   "BB",   "BB",   "AA",   "AA",   "AA",   "AA",   NA,  "AA",
#' "AA",   "AA",   "BB",   "BB",   "AA",   "BB",   "BB",   "BB",  "AB", "AA",
#' "AA",   "AA",   "BB",   "BB",   "AA",    NA,    "BB",   "AA",   NA,  "AA",
#' "AB",   "AB",   "BB",   "BB",   "BB",   "AA",   "BB",   "BB",   NA,  "AB",
#' "AA",   "AA",    NA,    "BB",    NA,    "AA",   "AA",   "AA",  "AA", "AA",
#' "AA",    NA,     NA,    "BB",   "BB",   "BB",   "BB",   "BB",  "AA", "AA",
#' "AA",    NA,    "AA",   "BB",   "BB",   "BB",   "AA",   "AA",   NA,  "AA"),
#' ncol = 10, byrow = TRUE, dimnames = list(paste0("IND", 1:10),
#' paste0("SNP", 1:10)))
#'
#' ## create object of class 'gData'.
#' gData <- createGData(geno = markers)
#'
#' ## Code markers by minor allele, no imputation.
#' gDataCoded1 <- codeMarkers(gData = gData, impute = FALSE)
#'
#' ## Code markers by reference alleles, impute missings by fixed value.
#' gDataCoded2 <- codeMarkers(gData = gData,
#'                            refAll = rep(x = c("A", "B"), times =  5),
#'                            impute = TRUE, imputeType = "fixed",
#'                            fixedValue = 1)
#'
#' ## Code markers by minor allele, impute by random value.
#' gDataCoded3 <- codeMarkers(gData = gData, impute = TRUE,
#'                            imputeType = "random")
#'
#' @importFrom utils read.table write.table
#' @export
codeMarkers <- function(gData,
                        refAll = "minor",
                        nMissGeno = 1,
                        nMiss = 1,
                        MAF = NULL,
                        MAC = NULL,
                        removeDuplicates = TRUE,
                        keep = NULL,
                        impute = TRUE,
                        imputeType =  c("random", "fixed", "beagle"),
                        fixedValue = NULL,
                        naStrings = NA,
                        verbose = FALSE) {
  ## Checks.
  chkGData(gData, comps = "markers")
  if (ncol(gData$markers) == 0) {
    stop("At least one marker should be present.\n")
  }
  if (length(refAll) > 1 && !length(refAll) == ncol(gData$markers)) {
    stop("Number of reference alleles should either be 1 or equal to",
         "the amount of SNPs in markers.\n")
  }
  chkNum(nMissGeno, min = 0, max = 1)
  chkNum(nMiss, min = 0, max = 1)
  if (!is.null(MAF) && !is.null(MAC)) {
    stop("Only one of MAF and MAC can be specified.\n")
  }
  if (!is.null(MAF)) {
    chkNum(MAF, min = 0, max = 1)
  }
  if (!is.null(MAC)) {
    chkNum(MAC, min = 0, max = 2 * nrow(gData$markers))
  }
  if (!is.null(keep) && (!is.character(keep) ||
                         !all(keep %in% colnames(gData$markers)))) {
    stop("All items in keep should be SNPs in markers.\n")
  }
  if (impute) {
    imputeType <- match.arg(imputeType)
    if (imputeType == "fixed") {
      if (is.null(fixedValue)) {
        stop("When using imputeType = fixed, fixedValue cannot be NULL.\n")
      }
      chkNum(fixedValue, min = 0, max = 2)
    }
    if (imputeType == "beagle" && !nzchar(Sys.which("java"))) {
      stop("When using beagle imputation an installation of java is ", 
           "required.\n")
    }
    if (imputeType == "beagle" &&
        (is.null(gData$map) || any(gData$map$pos != round(gData$map$pos)))) {
      stop("When using beagle imputation gData should contain a map with only ",
           "integer positions.\n")
    }
    if (imputeType == "beagle" && 
        ## Beagle needs at least two different values for position for each
        ## chromosome.
        min(by(data = gData$map$pos, INDICES = gData$map$chr, 
               FUN = function(x) {
                 length(unique(x))
               })) == 1) {
      stop("When using beagle imputation gData should contain a map with at ", 
           "least two different positions for each chromosome.\n")
    }
  }
  markersOrig <- gData$markers
  if (verbose) {
    codeSum <- paste0("Input contains ", ncol(markersOrig), " SNPs for ", 
                      nrow(markersOrig), " genotypes.\n")
  }
  ## Replace naStrings by NA.
  if (any(!is.na(naStrings))) {
    if (verbose) {
      codeSum <- c(codeSum, paste(sum(markersOrig %in% naStrings), 
                                  "values replaced by NA.\n"))
    }
    markersOrig[markersOrig %in% naStrings] <- NA
  }
  ## Remove genotypes with too many missings.
  if (!is.null(nMissGeno)) {
    genoMiss <- rowMeans(is.na(markersOrig)) < nMissGeno
    markersClean <- markersOrig[genoMiss, ]
    if (verbose) {
      codeSum <- c(codeSum, 
                   paste0(sum(!genoMiss), " genotypes removed because ", 
                          "proportion of missing values larger than or equal ",
                          "to ", nMissGeno, ".\n"))
    }
  }
  snpKeep <- colnames(markersClean) %in% keep
  ## Remove markers with too many missings.
  if (!is.null(nMiss)) {
    snpMiss <- colMeans(is.na(markersClean)) < nMiss
    markersClean <- markersClean[, snpMiss | snpKeep, drop = FALSE]
    if (verbose) {
      codeSum <- c(codeSum, 
                   paste0(sum(!(snpMiss | snpKeep)), " SNPs removed because ", 
                          "proportion of missing values larger than or equal ", 
                          "to ", nMiss, ".\n"))
    }
    snpKeep <- snpKeep[snpMiss | snpKeep]
    if (length(refAll) > 1) {
      refAll <- refAll[snpMiss | snpKeep]
    }
  }
  ## Recode markers.
  if (!is.numeric(markersClean)) {
    ## . is the default NA value in vcf format.
    ## Perform an extra check that it is being removed.
    if (any(markersClean == ".", na.rm = TRUE) && !"." %in% naStrings) {
      stop("SNPs contain '.', but these are not set to missing values.\n
           Specify '.' in naStrings to set them to missing values.") 
    }
    if (refAll[1] == "minor") {
      refAlls <- character()
    } else if (length(refAll) > 1) {
      refAlls <- refAll
    } else {
      refAlls <- rep(x = refAll, times = ncol(markersClean))
    }
    markersCharRecoded <- codeCharMarkers(markersClean, refAlls, 
                                          refAll[1] == "minor", NULL)
    markersRecoded <- markersCharRecoded$markersRecoded
    maxAll <- markersCharRecoded$maxAll
  } else { # Numeric input.
    markersRecoded <- markersClean
    maxAll <- max(markersRecoded, na.rm = TRUE)
  }
  ## Remove markers with low MAF.
  ## If MAC is specified convert it to MAF.
  if (!is.null(MAC)) {
    MAF <- MAC / (maxAll * nrow(markersRecoded)) - 1e-5
  }
  if (!is.null(MAF)) {
    snpMAFs <- colMeans(markersRecoded, na.rm = TRUE)
    snpMAF <- snpMAFs >= maxAll * MAF & snpMAFs <= maxAll * (1 - MAF)
    markersRecoded <- markersRecoded[, snpMAF | snpKeep, drop = FALSE]
    if (verbose) {
      codeSum <- c(codeSum, 
                   paste0(sum(!(snpMAF | snpKeep)), " SNPs removed because ", 
                          "MAF smaller than ", MAF, ".\n"))
    }
    snpKeep <- snpKeep[snpMAF | snpKeep]
  }
  ## Remove duplicated markers.
  ## Before duplicated was used for this but this only removes one occurrence
  ## per duplicate. Unique does not.
  if (removeDuplicates) {
    ## Only using unique would always remove the first occurrence.
    ## Using sample to make it random, always putting keep SNPs first.
    randOrder <- c(c(1:ncol(markersRecoded))[snpKeep],
                   sample(x = c(1:ncol(markersRecoded))[!snpKeep]))
    markersDedup <- unique(markersRecoded[, randOrder], MARGIN = 2)
    markersRecoded <- cbind(markersRecoded[, colnames(markersRecoded[, snpKeep])
                                           [!colnames(markersRecoded[, snpKeep])
                                             %in% colnames(markersDedup)],
                                           drop = FALSE],
                            markersDedup)
    snpKeep <- c(rep(TRUE, sum(snpKeep)),
                 rep(FALSE, ncol(markersRecoded) - sum(snpKeep)))
    if (verbose) {
      codeSum <- c(codeSum, 
                   paste(length(randOrder) - ncol(markersRecoded), 
                         "duplicate SNPs removed.\n"))
    }
  }
  ## Impute missing values.
  if (impute) {
    ## Count number of missings in markers to print after imputing.
    noMiss <- sum(is.na(markersRecoded))
    if (imputeType == "fixed") {
      ## Replace missing values by fixed value.
      markersRecoded[is.na(markersRecoded)] <- fixedValue
    } else if (imputeType == "random") {
      ## Replace missing values by random value based on probabilities per SNP.
      snpNA <- apply(X = markersRecoded, MARGIN = 2, FUN = anyNA)
      het <- maxAll != 2 || 
        (maxAll == 2 && any(markersRecoded == 1, na.rm = TRUE))
      markersRecoded[, snpNA] <-
        apply(X = markersRecoded[, snpNA, drop = FALSE], MARGIN = 2, 
              FUN = function(x) {
                p <- mean(x, na.rm = TRUE) / maxAll
                if (het) {
                  ## At least one heterozygous marker. Imputation may add more.
                  x[is.na(x)] <- sample(x = 0:maxAll, size = sum(is.na(x)),
                                        replace = TRUE,
                                        prob = choose(maxAll, 0:maxAll) *
                                          (1 - p) ^ (maxAll:0) * p ^ (0:maxAll))
                } else {
                  ## only homozygous markers.
                  x[is.na(x)] <- sample(x = c(0, maxAll), size = sum(is.na(x)),
                                        replace = TRUE, prob = c(1 - p, p))
                }
                return(x)
              })
    } else if (imputeType == "beagle") {
      markersRecoded <- imputeBeagle(markersRecoded = markersRecoded, 
                                     map = gData$map)
    }
    if (verbose) {
      ## Print after imputing.
      ## Mainly relevant for beagle since it may take very long to run.
      codeSum <- c(codeSum, paste(noMiss, "missing values imputed.\n"))
    }
    if (refAll[1] == "minor") {
      ## Correct for position of minor allele after imputation.
      markersRecoded[, colMeans(markersRecoded, na.rm = TRUE) > maxAll / 2] <-
        maxAll -
        markersRecoded[, colMeans(markersRecoded, na.rm = TRUE) > maxAll / 2]
    }
    ## Remove markers with low MAF after imputation.
    if (!is.null(MAF)) {
      snpMAFs <- colMeans(markersRecoded, na.rm = TRUE)
      snpMAF <- snpMAFs >= maxAll * MAF & snpMAFs <= maxAll * (1 - MAF)
      markersRecoded <- markersRecoded[, snpMAF | snpKeep, drop = FALSE]
      if (verbose) {
        codeSum <- c(codeSum, 
                     paste(sum(!(snpMAF | snpKeep)), "SNPs removed because", 
                           "MAF smaller than", MAF, "after imputation.\n"))
      }
      snpKeep <- snpKeep[snpMAF | snpKeep]
    }
    ## Remove duplicated markers after imputation.
    ## Before duplicated was used for this but this only removes one occurrence
    ## per duplicate. Unique does not.
    if (removeDuplicates) {
      ## Only using unique would always remove the first occurrence.
      ## Using sample to make it random, always putting keep SNPs first.
      randOrder <- c(c(1:ncol(markersRecoded))[snpKeep],
                     sample(x = c(1:ncol(markersRecoded))[!snpKeep]))
      markersDedup <- unique(markersRecoded[, randOrder], MARGIN = 2)
      markersRecoded <-
        cbind(markersRecoded[, colnames(markersRecoded[, snpKeep])
                             [!colnames(markersRecoded[, snpKeep]) %in%
                                 colnames(markersDedup)], drop = FALSE],
              markersDedup)
      if (verbose) {
        codeSum <- c(codeSum, 
                     paste(length(randOrder) - ncol(markersRecoded),
                           "duplicate SNPs removed after imputation.\n"))
      }
    }
  }
  if (removeDuplicates) {
    ## Columns are put in random order when removing duplicates.
    ## Reorder these columns.
    markersRecoded <-
      markersRecoded[, colnames(markersOrig)[colnames(markersOrig) %in%
                                               colnames(markersRecoded)],
                     drop = FALSE]
  }
  ## Remove removed SNPs and genotypes from map, pheno, kinship and covar.
  ## Replace corresponding slots in original gData object since recreating it 
  ## is both needlessly time consuming and produces warnings for replacing
  ## existing outputs.
  genoNew <- rownames(markersRecoded)
  if (!is.null(gData$map)) {
    mapNew <- gData$map[rownames(gData$map) %in% colnames(markersRecoded), ]
    gData$map <- mapNew
  }
  if (!is.null(gData$pheno)) {
    gData$pheno <- lapply(X = gData$pheno, FUN = function(trial) {
      trial[trial$genotype %in% genoNew, ]
    })
  }
  if (!is.null(gData$kinship)) {
    if (!is.list(gData$kinship)) {
      gData$kinship <-
        gData$kinship[rownames(gData$kinship) %in% genoNew,
                      colnames(gData$kinship) %in% genoNew]
    } else {
      gData$kinship <- lapply(X = gData$kinship, FUN = function(k) {
        k[rownames(k) %in% genoNew, colnames(k) %in% genoNew]
      })
    }
  }
  if (!is.null(gData$covar)) {
    gData$covar <- gData$covar[rownames(gData$covar) %in% genoNew, , 
                               drop = FALSE]
  }
  if (verbose) {
    codeSum <- c(codeSum, 
                 paste("Output contains", ncol(markersRecoded), "SNPs for",
                       nrow(markersRecoded), "genotypes.\n"))
    message(codeSum)
  }
  ## Return gData object with recoded and imputed markers.
  ## Replace corresponding slots in original gData object since recreating it 
  ## is both needlessly time consuming and produces warnings for replacing
  ## existing outputs.
  gData$markers <- markersRecoded
  return(gData)
}

#' Helper function for imputation using beagle.
#' 
#' @noRd
#' @keywords internal
imputeBeagle <- function(markersRecoded,
                         map) {
  ## Imputation of missing values using beagle software.
  ## Create temporary files for writing beagle input and output.
  outdir <- tempdir()
  tmpMap <- tempfile(tmpdir = outdir, fileext = ".map")
  tmpVcf <- tempfile(tmpdir = outdir, fileext = "input.vcf")
  tmpVcfOut <- tempfile(tmpdir = outdir, fileext = "out")
  ## Convert map to format suitable for beagle input.
  map <- map[colnames(markersRecoded), ]
  mapBeagle <- data.frame(map$chr, rownames(map), map$pos, map$pos)
  mapBeagle <- mapBeagle[order(mapBeagle$map.chr, mapBeagle$map.pos), ]
  ## Beagle doesn't accept duplicate map entries. To get around this
  ## duplicate entries are made distinct by adding 1 to the position.
  while (anyDuplicated(mapBeagle[, c(1, 4)])) {
    mapBeagle[duplicated(mapBeagle[, c(1, 4)]), 4] <-
      mapBeagle[duplicated(mapBeagle[, c(1, 4)]), 4] + 1
  }
  ## Convert chromosome column to a format suitable for beagle:
  ## Always 'chr' + chr number.
  if (!is.integer(mapBeagle[, 1])) {
    mapBeagle[, 1] <- as.integer(as.factor(mapBeagle[, 1]))
  }
  mapBeagle[, 1] <- paste0("chr", mapBeagle[, 1])
  ## Convert position columns to character format to avoid coversion
  ## to scientific format when writing to file.
  ## Beagle cannot handle scientific input format.
  mapBeagle[, 3] <- format(mapBeagle[, 3], trim = TRUE, scientific = FALSE)
  mapBeagle[, 4] <- format(mapBeagle[, 4], trim = TRUE, scientific = FALSE)
  ## Write map to .map file
  write.table(mapBeagle, file = tmpMap, col.names = FALSE, 
              row.names = FALSE, quote = FALSE, na = ".", sep = "\t")
  ## Convert markers to format suitable for beagle input.
  all00 <- "0/0"
  all01 <- "0/1"
  all11 <- "1/1"
  all10 <- "1/0"
  markersBeagle <- as.data.frame(t(markersRecoded),
                                 stringsAsFactors = FALSE)
  markersBeagle[markersBeagle == 0] <- all00
  markersBeagle[markersBeagle == 1] <- all01
  markersBeagle[markersBeagle == 2] <- all11
  markersBeagle[markersBeagle == -1] <- all10
  ## Write markers to vcf file.
  vcfBeagle <- cbind(data.frame(CHROM = mapBeagle[, 1], 
                                POS = mapBeagle[, 4], ID = rownames(map), 
                                REF = "A", ALT = "G", QUAL = ".", 
                                FILTER = "PASS", INFO = ".", FORMAT = "GT",
                                stringsAsFactors = FALSE),
                     markersBeagle)
  cat(file = tmpVcf,
      "##fileformat=VCFv4.1",
      "\n##filedate=", Sys.Date(),
      "\n##source=\"codeMarkers of statgenGWAS\"",
      "\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "\n#")
  cat(file = tmpVcf, paste(colnames(vcfBeagle), collapse = "\t"), "\n",
      append = TRUE)
  write.table(vcfBeagle, file = tmpVcf, quote = FALSE, col.names = FALSE,
              row.names = FALSE, append = TRUE, sep = "\t",
              na = paste(".", ".", sep = "/"))
  ## Run beagle with default settings.
  system(paste0("java -Xmx3000m -jar ",
                shQuote(paste0(sort(path.package()[grep("statgenGWAS",
                                                        path.package())])[1],
                               "/java/beagle.jar")), " gtgl=", tmpVcf, 
                " out=", tmpVcfOut, " gprobs=true seed=1234 nthreads=", 1,
                " map=", tmpMap), intern = TRUE)
  ## Read beagle output.
  beagleOut <- read.table(gzfile(paste0(tmpVcfOut, ".vcf.gz")),
                          stringsAsFactors = FALSE)
  ## Convert beagle output to format suitable for gData.
  markersRecoded <- t(beagleOut[, 10:ncol(beagleOut)])
  markersRecoded[substr(markersRecoded, 1, 3) == all00] <- 0
  markersRecoded[substr(markersRecoded, 1, 3) == all01] <- 1
  markersRecoded[substr(markersRecoded, 1, 3) == all11] <- 2
  markersRecoded[substr(markersRecoded, 1, 3) == all10] <- NA
  markersRecoded <- apply(X = markersRecoded, MARGIN = 2, FUN = as.numeric)
  markersRecoded[markersRecoded == -1] <- 0
  rownames(markersRecoded) <- colnames(markersBeagle)
  colnames(markersRecoded) <- beagleOut[, 3]
  return(markersRecoded)
}



