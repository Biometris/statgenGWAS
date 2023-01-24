#' Annotate significant SNPs
#' 
#' Annotation of significant SNPs with gene info.
#' 
#' @export
annotateSNP <- function(GWAS,
                        gff) {
  signSnp <- do.call(rbind, GWAS$signSnp)
  ## Read gff (https://en.wikipedia.org/wiki/General_feature_format).
  gffDat <- data.table::fread(gff, skip = 1, na.strings = c("###", "."), 
                              fill = TRUE, sep = "\t")
  ## Set colnames as described on wiki.
  colnames(gffDat) <- c("sequence", "source", "feature", "start", "end", "score",
                        "strand", "phase", "attributes")
  ## Remove empty sequences and select genes.
  gffDat <- gffDat[!is.na(gffDat[["sequence"]]), ]
  gffDat <- gffDat[gffDat[["feature"]] == "gene"]
  ## Replace chromosomes by their numeric value.
  ## This might give warnings for non-numeric chromosomes.
  ## These are converted to NA and have to be removed.
  gffDat[["chr"]] <- suppressWarnings(as.numeric(gsub(pattern = "Chr", 
                                                      replacement = "",
                                                      x = gffDat[["sequence"]])))
  gffDat <- gffDat[!is.na(gffDat[["chr"]]), ]
  ## Not allocating memory causes an error in the merge step.
  ## This is a known problem in data.table.
  data.table::setalloccol(gffDat, 2048)
  annoDat <- gffDat[signSnp, on = c("chr", "start <= pos", "end >= pos"), 
                    nomatch = NA]
  ## Process attributes to extract gene IDs.
  ## Attributes are assumed to be a ; separated list with items of the form a=b. 
  subAttribs <- strsplit(annoDat[["attributes"]], split = ";")
  annoDat[["ID"]] <- sapply(X = subAttribs, FUN = function(subAttrib) {
    if (!is.na(subAttrib[1])) {
      items <- strsplit(subAttrib, split = "=")
      items <- setNames(sapply(items, `[[`, 2), sapply(items, `[[`, 1))
      items[["ID"]]
    } else {
      NA
    }
  })
  ## Convert to list with genes per SNP.
  annoDat2 <- split(annoDat[["ID"]], annoDat[["snp"]])
  return(annoDat2)
}


