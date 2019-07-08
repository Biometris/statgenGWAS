context("codeMarkers")

## Create object containing character markers.
markers <- matrix(c("AA", "AB", "AA", "AB", NA, "AA", NA, "AA",
                    "AA", "BB", "AB", "BB"), ncol = 4, byrow = TRUE,
                  dimnames = list(paste0("IND", 1:3), paste0("SNP", 1:4)))
gData <- createGData(geno = markers)
## Create object containing numeric markers.
markers2 <- matrix(c(0, 1, 0, 1, NA, 0, NA, 0, 0, 2, 1, 2), ncol = 4,
                   byrow = TRUE,
                   dimnames = list(paste0("IND", 1:3), paste0("SNP", 1:4)))
gData2 <- createGData(geno = markers2)

test_that(paste("codeMarkers produces correct results when no imputation,",
                "remove duplicates"), {
                  expect_equal(as.numeric(codeMarkers(gData = gData,
                                                      impute = FALSE)$markers),
                               c(0, NA, 0, 0, NA, 1, 1, 2, 0))
                  expect_equal(as.numeric(codeMarkers(gData = gData, nMiss = 0.1,
                                                      impute = FALSE)$markers),
                               c(1, 2, 0))
                  expect_equal(as.numeric(codeMarkers(gData = gData, MAF = 0.2,
                                                      impute = FALSE)$markers),
                               c(0, NA, 1, 1, 2, 0))
                  expect_equal(as.numeric(codeMarkers(gData = gData, nMiss = 0.1,
                                                      MAF = 0.2,
                                                      impute = FALSE)$markers),
                               c(1, 2, 0))
                })

test_that(paste("codeMarkers produces correct results when no imputation,",
                "no remove duplicates"), {
                  expect_equal(as.numeric(codeMarkers(gData = gData,
                                                      removeDuplicates = FALSE,
                                                      impute = FALSE)$markers),
                               c(0, NA, 0, 1, 2, 0, 0, NA, 1, 1, 2, 0))
                  expect_equal(as.numeric(codeMarkers(gData = gData, nMiss = 0.1,
                                                      removeDuplicates = FALSE,
                                                      impute = FALSE)$markers),
                               c(1, 2, 0, 1, 2, 0))
                  expect_equal(as.numeric(codeMarkers(gData = gData, MAF = 0.2,
                                                      removeDuplicates = FALSE,
                                                      impute = FALSE)$markers),
                               c(1, 2, 0, 0, NA, 1, 1, 2, 0))
                  expect_equal(as.numeric(codeMarkers(gData = gData, nMiss = 0.1,
                                                      MAF = 0.2,
                                                      removeDuplicates = FALSE,
                                                      impute = FALSE)$markers),
                               c(1, 2, 0, 1, 2, 0))
                })

test_that("codeMarkers produces correct results with imputation", {
  expect_equal(as.numeric(codeMarkers(gData = gData, impute = TRUE,
                                      imputeType = "fixed",
                                      fixedValue = 1)$markers),
               c(0, 1, 0, 0, 1, 1, 1, 2, 0))
  expect_equal(as.numeric(codeMarkers(gData = gData, impute = TRUE,
                                      imputeType = "random")$markers),
               c(0, 0, 0, 0, 1, 1, 1, 2, 0))
  expect_equal(as.numeric(codeMarkers(gData = gData, removeDuplicates = FALSE,
                                      impute = TRUE, imputeType = "fixed",
                                      fixedValue = 1)$markers),
               c(0, 1, 0, 1, 2, 0, 0, 1, 1, 1, 2, 0))
  expect_equal(as.numeric(codeMarkers(gData = gData, removeDuplicates = FALSE,
                                      impute = TRUE,
                                      imputeType = "random")$markers),
               c(0, 0, 0, 1, 2, 0, 0, 1, 1, 1, 2, 0))
})

test_that("option keep in codeMarkers works properly", {
  expect_equal(as.numeric(codeMarkers(gData = gData, nMiss = 0.1, keep = "SNP1",
                                      impute = FALSE)$markers),
               c(0, NA, 0, 1, 2, 0))
  expect_equal(as.numeric(codeMarkers(gData = gData, MAF = 0.2, keep = "SNP1",
                                      impute = FALSE)$markers),
               c(0, NA, 0, 0, NA, 1, 1, 2, 0))
  expect_equal(as.numeric(codeMarkers(gData = gData, keep = "SNP2",
                                      impute = FALSE)$markers),
               c(0, NA, 0, 1, 2, 0, 0, NA, 1))
  expect_equal(as.numeric(codeMarkers(gData = gData, keep = c("SNP2", "SNP4"),
                                      impute = TRUE, imputeType = "fixed",
                                      fixedValue = 1)$markers),
               c(0, 1, 0, 1, 2, 0, 0, 1, 1, 1, 2, 0))
})

test_that("option refAll in codeMarkers works properly", {
  expect_equal(as.numeric(codeMarkers(gData = gData, refAll = "A",
                                      impute = FALSE)$markers),
               c(2, NA, 2, 2, NA, 1, 1, 2, 0))
  expect_equal(as.numeric(codeMarkers(gData = gData,
                                      refAll = c("A", "B", "B", "A"),
                                      impute = FALSE)$markers),
               c(2, NA, 2, 1, 0, 2, 0, NA, 1, 1, 2, 0))
  expect_error(codeMarkers(gData = gData, refAll = c("A", "B"), impute = FALSE))
})

test_that("codeMarkers functions properly when using numeric input", {
  expect_equal(as.numeric(codeMarkers(gData = gData2, impute = FALSE)$markers),
               c(0, NA, 0, 0, NA, 1, 1, 0, 2))
  expect_equal(as.numeric(codeMarkers(gData = gData2, removeDuplicates = FALSE,
                                      impute = TRUE, imputeType = "fixed",
                                      fixedValue = 1)$markers),
               c(0, 1, 0, 1, 0, 2, 0, 1, 1, 1, 0, 2))
})

test_that("option verbose in codeMarkers works properly", {
  expect_equal(capture.output(gDataRec <- codeMarkers(gData = gData, MAF = 0.2,
                                                      impute = FALSE)),
               character())
  consOut <- capture.output(gDataRec <- codeMarkers(gData = gData, MAF = 0.2,
                                                    impute = FALSE,
                                                    verbose = TRUE))
  expect_true(any(grepl(pattern = "4 markers for 3 genotypes", x = consOut)))
  expect_true(any(grepl(pattern = "1 markers removed because MAF",
                        x = consOut)))
  expect_true(any(grepl(pattern = "1 duplicate markers removed", x = consOut)))
  expect_true(any(grepl(pattern = "2 markers for 3 genotypes", x = consOut)))
})

test_that("codeMarkers always removes markers with only NA", {
  gData$markers[, 1] <- NA
  expect_equal(as.numeric(codeMarkers(gData = gData, impute = FALSE)$markers),
               c(0, NA, 1, 1, 2, 0))
})

test_that("codeMarkers always removes genotypes with only NA", {
  gData$markers[1, ] <- NA
  expect_equal(as.numeric(codeMarkers(gData = gData, impute = FALSE)$markers),
               c(NA, 0, NA, 1, 2, 0))
})
