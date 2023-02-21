#' Annotate significant SNPs
#' 
#' Annotation of significant SNPs with gene info.
#' 
#' @export
annotateSNP <- function(GWAS,
                        gff) {
  signSnp <- do.call(rbind, GWAS$signSnp)
  ## Read gff (https://en.wikipedia.org/wiki/General_feature_format).
  # gffDat <- data.table::fread(gff, skip = 1, na.strings = c("###", "."), 
  #                             fill = TRUE, sep = "\t")
  # ## Set colnames as described on wiki.
  # colnames(gffDat) <- c("seqid", "source", "type", "start", "end", "score",
  #                       "strand", "phase", "attributes")
  
  gffDat <- ape::read.gff(gff)
  gffDat <- data.table::as.data.table(gffDat)
  
  ## Remove empty sequences and select genes.
  gffDat <- gffDat[!is.na(gffDat[["seqid"]]), ]
  gffDat <- gffDat[gffDat[["type"]] == "gene", ]
  ## Replace chromosomes by their numeric value.
  ## This might give warnings for non-numeric chromosomes.
  ## These are converted to NA and have to be removed.
  gffDat[["chr"]] <- gsub(pattern = "Chr", replacement = "",
                          x = gffDat[["seqid"]])
  gffDat <- gffDat[!is.na(gffDat[["chr"]]), ]
  
  signSnp[["chr"]] <- gsub(pattern = "chr", replacement = "",
                           x = signSnp[["chr"]])
  
  ## Not allocating memory causes an error in the merge step.
  ## This is a known problem in data.table.
  data.table::setalloccol(gffDat, 2048)
  annoDat <- gffDat[signSnp, on = c("chr", "start <= pos", "end >= pos"), 
                    nomatch = NULL]
  ## Process attributes to extract gene IDs.
  ## Attributes are assumed to be a ; separated list with items of the form a=b. 
  subAttribs <- strsplit(annoDat[["attributes"]], split = ";")
  attrDat <- sapply(X = subAttribs, FUN = function(subAttrib) {
    if (!is.na(subAttrib[1])) {
      items <- strsplit(subAttrib, split = "=")
      items <- setNames(sapply(items, `[[`, 2), sapply(items, `[[`, 1))
      as.data.frame(t(items))
    } else {
      data.frame()
    }
  })
  attrDat <- do.call(dfBind, args = list(dfList = attrDat))
  annoDat <- cbind(annoDat, attrDat)
  ## Convert to list with genes per SNP.
  return(annoDat)
}


