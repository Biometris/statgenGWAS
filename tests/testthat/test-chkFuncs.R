context("check Functions")

set.seed(1234)
X <- matrix(sample(x = c(0, 1), size = 30, replace = TRUE), nrow = 10)
Sigma <- matrix(runif(n = 100), nrow = 10)
Sigma <- tcrossprod(Sigma)
covs <- matrix(runif(n = 20, max = 100), nrow = 10)
pheno <- data.frame(genotype = paste0("G", 1:10),
                    matrix(rnorm(50, mean = 10, sd = 2), nrow = 10))
map <- data.frame(chr = c(1, 1, 2), pos = 1:3)
rownames(X) <- rownames(Sigma) <- colnames(Sigma) <- rownames(covs) <-
  paste0("G", 1:10)
colnames(X) <- rownames(map) <- paste0("M", 1:3)
gDataTest <- createGData(map = map, geno = X, kin = Sigma,
                         pheno = list(ph1 = pheno, ph2 = pheno),
                         covar = as.data.frame(covs))

test_that("chkGData functions properly", {
  expect_error(chkGData(), "valid gData object")
  expect_error(chkGData(X), "valid gData object")
  expect_silent(chkGData(gDataTest))
  gDataTest2 <- gDataTest
  gDataTest2$markers <- NULL
  expect_error(chkGData(gDataTest2), "should contain markers")
  expect_silent(chkGData(gDataTest2, comps = "map"))
})

test_that("chkEnvs functions properly", {
  expect_silent(chkEnvs("ph1", gDataTest))
  expect_silent(chkEnvs(c("ph1", "ph2"), gDataTest))
  expect_error(chkEnvs("ph3", gDataTest), "should be in pheno")
  expect_silent(chkEnvs(1, gDataTest))
  expect_silent(chkEnvs(c(1, 2), gDataTest))
  expect_error(chkEnvs(3, gDataTest), "should be in pheno")
})

test_that("chkTraits functions properly", {
  expect_silent(chkTraits("X1", "ph1", gDataTest))
  expect_silent(chkTraits(c("X1", "X2"), "ph1", gDataTest))
  expect_silent(chkTraits(c("X1", "X2"), c("ph1", "ph2"), gDataTest))
  expect_error(chkTraits("X6", "ph1", gDataTest), "For ph1 not all traits")
  expect_error(chkTraits(c("X1", "X6"), "ph1", gDataTest),
               "For ph1 not all traits")
  expect_silent(chkTraits(2, 1, gDataTest))
  expect_silent(chkTraits(c(2, 3), 1, gDataTest))
  expect_silent(chkTraits(c(2, 3), c(1, 2), gDataTest))
  expect_error(chkTraits(1, 1, gDataTest), "For 1 not all traits")
  expect_error(chkTraits(c(2, 7), 1, gDataTest), "For 1 not all traits")
})

test_that("chkNum functions properly", {
  expect_silent(chkNum(1, min = 0, max = 2))
  expect_error(chkNum("a"), "single numerical value")
  expect_error(chkNum(c(1, 2)), "single numerical value")
  expect_error(chkNum(1, min = NULL, max = 0), "smaller than 0")
  expect_error(chkNum(1, min = 2, max = NULL), "greater than 2")
  expect_error(chkNum(1, min = 2, max = 3), "between 2 and 3")
})
