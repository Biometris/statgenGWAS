#' Read PLINK binary data 
#' 
#' Read PLINK binary data and save in gData format. This is a wrapper around 
#' \code{\link[snpStats]{read.plink}} in the Bioconductor package 
#' \code{snpStats}. This package needs to be installed for the function to 
#' work. 
#' 
#' @param bed The name of the file containing the packed binary SNP genotype 
#' data. It should have the extension .bed; If it doesn't, then this extension
#' will be appended.
#' @param bim The file containing the SNP descriptions. If not specified 
#' \code{bed} is used with its file extension replaced by bim.
#' @param fam The file containing subject (and, possibly, family) identifiers. 
#' This is basically a tab-delimited "pedfile". If not specified 
#' \code{bed} is used with its file extension replaced by fam.
#' @param ... Further arguments passed to \code{\link[snpStats]{read.plink}}.
#' 
#' @return An object of class \code{gData}.
#' 
#' @export
readPLINK <- function(bed,
                      bim,
                      fam,
                      ...) {
  if (!requireNamespace("snpStats", quietly = TRUE)) {
    stop("Package snpStats needed for reading PLINK files.")
  }
  genoPLINK <- snpStats::read.plink(bed = bed, bim = bim, fam = fam, ...)
  
  ## Get makers.  
  markers <- as(genoPLINK$genotypes, "numeric")
  ## Get map.
  map <- genoPLINK$map
  names(map) <- c("chr", "snp.name", "cM", "pos", "allele.1", "allele.2")
  # Create a gData object containing map and marker information.
  res <- createGData(geno = markers, map = map)
  return(res)
}
