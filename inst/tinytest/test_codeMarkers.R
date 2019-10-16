### Test codeMarkers

# Create object containing character markers.
markers <- matrix(c("AA", "AB", "AA", "AB", NA, "AA", NA, "AA",
                    "AA", "BB", "AB", "BB"), ncol = 4, byrow = TRUE,
                  dimnames = list(paste0("IND", 1:3), paste0("SNP", 1:4)))
gData <- createGData(geno = markers)
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
                  matrix(c(0, NA, 1, 1, 2, 0), nrow = 3))
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
expect_equivalent(codeMarkers(gData = gData, impute = TRUE, 
                              imputeType = "fixed", fixedValue = 1)$markers,
                  matrix(c(0, 1, 0, 0, 1, 1, 1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, impute = TRUE,
                              imputeType = "random")$markers,
                  matrix(c(0, 0, 0, 0, 1, 1, 1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, removeDuplicates = FALSE,
                              impute = TRUE, imputeType = "fixed",
                              fixedValue = 1)$markers,
                  matrix(c(0, 1, 0, 1, 2, 0, 0, 1, 1, 1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, removeDuplicates = FALSE,
                              impute = TRUE,  imputeType = "random")$markers,
                  matrix(c(0, 0, 0, 1, 2, 0, 0, 1, 1, 1, 2, 0), nrow = 3))

# Check option keep.
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

# Check option refAll.
expect_error(codeMarkers(gData = gData, refAll = c("A", "B"), impute = FALSE),
             "number of reference alleles should either be")
expect_equivalent(codeMarkers(gData = gData, refAll = "A",
                              impute = FALSE)$markers,
                  matrix(c(2, NA, 2, 2, NA, 1, 1, 2, 0), nrow = 3))
expect_equivalent(codeMarkers(gData = gData, refAll = c("A", "B", "B", "A"),
                              impute = FALSE)$markers,
                  matrix(c(2, NA, 2, 1, 0, 2, 0, NA, 1, 1, 2, 0), nrow = 3))

# Check for using numrical input.
expect_equivalent(codeMarkers(gData = gData2, impute = FALSE)$markers,
                  matrix(c(0, NA, 0, 0, NA, 1, 1, 0, 2), nrow = 3))
expect_equivalent(codeMarkers(gData = gData2, removeDuplicates = FALSE,
                              impute = TRUE, imputeType = "fixed",
                              fixedValue = 1)$markers,
                  matrix(c(0, 1, 0, 1, 0, 2, 0, 1, 1, 1, 0, 2), nrow = 3))

# Check option verbose.
expect_equal(capture.output(gDataRec <- codeMarkers(gData = gData, MAF = 0.2,
                                                    impute = FALSE)),
             character())
consOut <- capture.output(gDataRec <- codeMarkers(gData = gData, MAF = 0.2,
                                                  impute = FALSE,
                                                  verbose = TRUE))
expect_true(any(grepl(pattern = "4 markers for 3 genotypes", x = consOut)))
expect_true(any(grepl(pattern = "1 markers removed because MAF", x = consOut)))
expect_true(any(grepl(pattern = "1 duplicate markers removed", x = consOut)))
expect_true(any(grepl(pattern = "2 markers for 3 genotypes", x = consOut)))

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
