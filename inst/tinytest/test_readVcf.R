### Test reading VCF files.

if (requireNamespace("vcfR")) {
  ## Use data from vcfR for testing.
  data(vcfR_test, package = "vcfR")
  
  ## Write data to temp directory.
  filePath <- file.path(tempdir(), "vcfRTest.vcf.gz")
  vcfR::write.vcf(vcfR_test, file = filePath)
  
  ## Read VCF file.
  expect_warning(gdVcf <- readVcf(filePath, verbose = FALSE),
                 "Not all markers in geno are in map. Extra markers")
  
  expect_inherits(gdVcf, "gData")
  
  expect_equal(dim(gdVcf[["map"]]), c(3, 2))
  expect_equal(dim(gdVcf[["markers"]]), c(3, 3))
}
