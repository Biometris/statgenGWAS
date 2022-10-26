### Test codeMarkers

set.seed(1234)

# Create object containing character markers.
markers <- matrix(c("AA", "AB", "AA", "AB", NA, "AA", NA, "AA",
                    "AA", "BB", "AB", "BB"), ncol = 4, byrow = TRUE,
                  dimnames = list(paste0("IND", 1:3), paste0("SNP", 1:4)))
map <- data.frame(chr = c(1, 1, 2, 2), pos = 1:4, 
                  row.names = paste0("SNP", 1:4))
pheno <- data.frame(genotype = paste0("IND", 1:3), t1 = 1) 
kin <- matrix(1, nrow = 3, ncol = 3, 
              dimnames = list(paste0("IND", 1:3), paste0("IND", 1:3)))
covar <- data.frame(C1 = 1:3, row.names = paste0("IND", 1:3))
gData <- createGData(geno = markers, map = map, kin = kin, 
                     pheno = pheno, covar = covar)
                     
# Create object containing numeric markers.
markers2 <- matrix(c(0, 1, 0, 1, NA, 0, NA, 0, 0, 2, 1, 2), ncol = 4,
                   byrow = TRUE,
                   dimnames = list(paste0("IND", 1:3), paste0("SNP", 1:4)))
gData2 <- createGData(geno = markers2)

# Check with no imputation, remove duplicates.

expect_equivalent(codeMarkers(gData = gData, impute = FALSE)$markers,
                  matrix(c(0, NA, 0, 0, NA, 1, 1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, nMiss = 0.1, 
                              impute = FALSE)$markers,
                  matrix(c(1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, MAF = 0.2, impute = FALSE)$markers,
                  matrix(c(1, 2, 0, 0, NA, 1), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, nMiss = 0.1, MAF = 0.2, 
                              impute = FALSE)$markers,
                  matrix(c(1, 2, 0), nrow = 3))

# Check with no imputation, no remove duplicates.

expect_equivalent(codeMarkers(gData = gData, removeDuplicates = FALSE,
                              impute = FALSE)$markers,
                  matrix(c(0, NA, 0, 1, 2, 0, 0, NA, 1, 1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, nMiss = 0.1,
                              removeDuplicates = FALSE, impute = FALSE)$markers,
                  matrix(c(1, 2, 0, 1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, MAF = 0.2, 
                              removeDuplicates = FALSE, impute = FALSE)$markers,
                  matrix(c(1, 2, 0, 0, NA, 1, 1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, nMiss = 0.1, MAF = 0.2,
                              removeDuplicates = FALSE, impute = FALSE)$markers,
                  matrix(c(1, 2, 0, 1, 2, 0), nrow = 3))

# Check with imputation.

expect_error(codeMarkers(gData = gData, impute = TRUE, imputeType = "fixed"),
             "When using imputeType = fixed, fixedValue cannot be NULL")
expect_error(codeMarkers(gData = gData, impute = TRUE, imputeType = "fixed",
                         fixedValue = 3),
             "fixedValue should be a single numerical value between 0 and 2")

expect_equivalent(codeMarkers(gData = gData, impute = TRUE, 
                              imputeType = "fixed", fixedValue = 1)$markers,
                  matrix(c(0, 1, 0, 1, 2, 0, 0, 1, 1), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, impute = TRUE,
                              imputeType = "random")$markers,
                  matrix(c(0, 0, 0, 0, 0, 1, 1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, removeDuplicates = FALSE,
                              impute = TRUE, imputeType = "fixed",
                              fixedValue = 1)$markers,
                  matrix(c(0, 1, 0, 1, 2, 0, 0, 1, 1, 1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, removeDuplicates = FALSE,
                              impute = TRUE,  imputeType = "random")$markers,
                  matrix(c(0, 0, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0), nrow = 3))

# Beagle imputation cannot be checked on CRAN since it needs 
# java to be installed. 

if (at_home()) {
  # map has to be present for beagle.
  gData3 <- createGData(geno = markers)
  expect_error(codeMarkers(gData = gData3, impute = TRUE, 
                           imputeType = "beagle"),
               "When using beagle imputation gData should contain a map")
  
  # map can only contain integer positions for beagle.
  map2 <- data.frame(chr = c(1, 1, 2, 2), pos = c(1, 1.1, 2, 2))
  gData4 <- suppressWarnings(createGData(gData = gData, map = map2))
  expect_error(codeMarkers(gData = gData4, impute = TRUE, 
                           imputeType = "beagle"),
               "gData should contain a map with only integer positions")
  
  # map needs at least two positions per chromosome for beagle.
  map3 <- data.frame(chr = c(1, 1, 1, 2), pos = 1:4)
  gData5 <- suppressWarnings(createGData(gData = gData, map = map3))
  expect_error(codeMarkers(gData = gData5, impute = TRUE, 
                           imputeType = "beagle"),
               "at least two different positions for each chromosome")
  
  expect_equivalent(codeMarkers(gData = gData, impute = TRUE, 
                                removeDuplicates = FALSE, 
                                imputeType = "beagle")$markers, 
                    matrix(c(0, 0, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0), nrow = 3))
}

# Check option keep.

expect_error(codeMarkers(gData = gData, keep = "SNP5"),
             "All items in keep should be SNPs in markers")

expect_equivalent(codeMarkers(gData = gData, nMiss = 0.1, keep = "SNP1",
                              impute = FALSE)$markers,
                  matrix(c(0, NA, 0, 1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, MAF = 0.2, keep = "SNP1",
                              impute = FALSE)$markers,
                  matrix(c(0, NA, 0, 0, NA, 1, 1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, keep = "SNP2",
                              impute = FALSE)$markers,
                  matrix(c(0, NA, 0, 1, 2, 0, 0, NA, 1), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, keep = c("SNP2", "SNP4"),
                              impute = TRUE, imputeType = "fixed",
                              fixedValue = 1)$markers,
                  matrix(c(0, 1, 0, 1, 2, 0, 0, 1, 1, 1, 2, 0), nrow = 3))

# Check option naStrings.
gData6 <- gData
gData6$markers[is.na(gData6$markers)] <- "*"
expect_equal(codeMarkers(gData6, naStrings = "*", impute = FALSE, 
                         removeDuplicates = FALSE), 
             codeMarkers(gData, impute = FALSE, removeDuplicates = FALSE))

# Check option refAll.

expect_error(codeMarkers(gData = gData, refAll = c("A", "B"), impute = FALSE),
             "Number of reference alleles should either be")
expect_equivalent(codeMarkers(gData = gData, refAll = "A",
                              impute = FALSE)$markers,
                  matrix(c(2, NA, 2, 2, NA, 1, 1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, refAll = c("A", "B", "B", "A"),
                              impute = FALSE)$markers,
                  matrix(c(2, NA, 2, 1, 0, 2, 0, NA, 1, 1, 2, 0), nrow = 3))

# Check option MAC

expect_error(codeMarkers(gData = gData, MAF = 0.5, MAC = 3),
             "Only one of MAF and MAC can be specified")
expect_error(codeMarkers(gData = gData, MAC = 10),
             "MAC should be a single numerical value between 0 and 6")

## MAC of 4 should correspond to MAF of 2/3.
expect_equal(codeMarkers(gData = gData, MAF = 2/3, impute = FALSE),
             codeMarkers(gData = gData, MAC = 4, impute = FALSE))

# Check for using numerical input.

expect_equivalent(codeMarkers(gData = gData2, impute = FALSE)$markers,
                  matrix(c(0, NA, 0, 0, NA, 1, 1, 0, 2), nrow = 3))
expect_equivalent(codeMarkers(gData = gData2, removeDuplicates = FALSE,
                              impute = TRUE, imputeType = "fixed",
                              fixedValue = 1)$markers,
                  matrix(c(0, 1, 0, 1, 0, 2, 0, 1, 1, 1, 0, 2), nrow = 3))

# Check option verbose.

expect_equal(capture.output(gd <- codeMarkers(gData = gData, MAF = 0.2,
                                              impute = FALSE)),
             character())
expect_message(codeMarkers(gData = gData, MAF = 0.2,
                           impute = FALSE, verbose = TRUE), 
               "4 SNPs for 3 genotypes")
expect_message(codeMarkers(gData = gData, MAF = 0.2,
                           impute = FALSE, verbose = TRUE), 
               "1 SNPs removed because MAF")
expect_message(codeMarkers(gData = gData, MAF = 0.2,
                           impute = FALSE, verbose = TRUE), 
               "1 duplicate SNPs removed")
expect_message(codeMarkers(gData = gData, MAF = 0.2,
                           impute = FALSE, verbose = TRUE), 
               "2 SNPs for 3 genotypes")

expect_message(codeMarkers(gData = gData6, MAF = 0.2, imputeType = "fixed",
                           fixedValue = 1, naStrings = "*", verbose = TRUE),
               "2 values replaced by NA")
expect_message(codeMarkers(gData = gData6, MAF = 0.2, imputeType = "fixed",
                           fixedValue = 1, naStrings = "*", verbose = TRUE),
               "1 missing values imputed")
expect_message(codeMarkers(gData = gData6, MAF = 0.2, imputeType = "fixed",
                           fixedValue = 1, naStrings = "*", verbose = TRUE),
               "0 SNPs removed because MAF")
expect_message(codeMarkers(gData = gData6, MAF = 0.2, imputeType = "fixed",
                           fixedValue = 1, naStrings = "*", verbose = TRUE),
               "0 duplicate SNPs removed after imputation")

# Check that markers with only NA are always removed.

gDataNA <- gData
gDataNA$markers[, 1] <- NA
expect_equivalent(codeMarkers(gData = gDataNA, impute = FALSE)$markers,
                  matrix(c(0, NA, 1, 1, 2, 0), nrow = 3))

# Check that genotypes with only NA are always removed.

gDataNA <- gData
gDataNA$markers[1, ] <- NA
expect_equivalent(codeMarkers(gData = gDataNA, impute = FALSE)$markers,
                  matrix(c(NA, 0, NA, 1, 2, 0), nrow = 2))

# Check that genotypes and individuals are removed from 
# map, pheno, kinship and covar as well.
gd <- codeMarkers(gData = gData, nMiss = 0.1, nMissGeno = 0.4, keep = "SNP1",
                  impute = FALSE)
expect_equal(rownames(gd$map), c("SNP1", "SNP2", "SNP3"))
expect_equal(gd$pheno$pheno$genotype, c("IND1", "IND3"))
expect_equal(rownames(gd$kinship), c("IND1", "IND3"))
expect_equal(colnames(gd$kinship), c("IND1", "IND3"))
expect_equal(rownames(gd$covar), c("IND1", "IND3"))


