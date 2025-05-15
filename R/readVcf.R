#' Read variant call format data 
#' 
#' Read variant call format (VCF) data and save in gData format. This is a 
#' wrapper around [vcfR::read.vcfR] in the package \code{vcfR}. This package 
#' needs to be installed for the function to work. 
#' 
#' @param vcfFile The name of the vcf file. This can either be a plain text file 
#' with extension (.vcf), or a gzipped file with extencion (.vcf.gz).
#' @param ... Further arguments passed to [vcfR::read.vcfR].
#' 
#' @returns An object of class \code{gData}.
#' 
#' @references Knaus BJ, Grünwald NJ (2017). “VCFR: a package to manipulate and 
#' visualizevariant call format data in R.” _Molecular Ecology Resources_, 
#' *17*(1), 44-53. ISSN 757, \doi{10.1111/1755-0998.12549}.
#' 
#' @export
readVcf <- function(vcfFile,
                    ...) {
  if (!requireNamespace("vcfR", quietly = TRUE)) {
    stop("Package vcfR needed for reading VCF files.")
  }
  genoVcf <- vcfR::read.vcfR(vcfFile, ...)
  ## Get makers.  
  markers <- t(vcfR::extract.gt(genoVcf, element = "GT", as.numeric = TRUE))
  ## Get map.
  map <- vcfR::getFIX(genoVcf)
  map <- map[!is.na(map[, "ID"]), ]
  rownames(map) <- map[, "ID"]
  map <- as.data.frame(map[, c("CHROM", "POS")])
  colnames(map) <- c("chr", "pos")
  map[["pos"]] <- as.numeric(map[["pos"]])
  # Create a gData object containing map and marker information.
  res <- createGData(geno = markers, map = map)
  return(res)
}
