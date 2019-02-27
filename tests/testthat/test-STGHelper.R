context("STG Helper functions")

set.seed(1234)
y <- 1:10
X <- matrix(sample(x = c(0, 1), size = 30, replace = TRUE), nrow = 10)
Sigma <- matrix(runif(n = 100), nrow = 10)
Sigma <- Sigma %*% t(Sigma)
covs <- matrix(runif(n = 20, max = 100), nrow = 10)
pheno <- data.frame(genotype = paste0("G", 1:10),
                    matrix(rnorm(50, mean = 10, sd = 2), nrow = 10))
map <- data.frame(chr = c(1, 1, 2), pos = 1:3)
rownames(X) <- rownames(Sigma) <- colnames(Sigma) <- rownames(covs) <-
  paste0("G", 1:10)
colnames(X) <- rownames(map) <- paste0("M", 1:3)
gDataTest <- createGData(map = map, geno = X, kin = Sigma, pheno = pheno,
                         covar = as.data.frame(covs))

test_that("EMMA produces correct results with default settings", {
  expect_equivalent(EMMA(gData = gDataTest, trait = 2, environment = 1)[[1]],
                    c(0.020597492367456, 1.85412717490278))
  expect_equivalent(EMMA(gData = gDataTest, trait = "X1",
                         environment = 1)[[1]],
                    c(0.020597492367456, 1.85412717490278))
  expect_equivalent(EMMA(gData = gDataTest, trait = "X1",
                         environment = "pheno")[[1]],
                    c(0.020597492367456, 1.85412717490278))
  expect_equivalent(EMMA(gData = gDataTest, trait = 2,
                         environment = 1, K = Sigma)[[1]],
                    c(0.020597492367456, 1.85412717490278))
})

test_that("EMMA produces correct results with covariates", {
  expect_equivalent(EMMA(gData = gDataTest, trait = 2,
                         environment = 1, covar = 1)[[1]],
                    c(8.76962021955844e-05, 1.93163739799549))
  expect_equivalent(EMMA(gData = gDataTest, trait = 2,
                         environment = 1, covar = "V1")[[1]],
                    c(8.76962021955844e-05, 1.93163739799549))
  expect_equivalent(EMMA(gData = gDataTest, trait = 2,
                         environment = 1, K = Sigma, covar = 1)[[1]],
                    c(8.76962021955844e-05, 1.93163739799549))
  expect_equivalent(EMMA(gData = gDataTest, trait = 2,
                         environment = 1, snpName = "M1")[[1]],
                    c(0.184076137051078, 1.83600897652493))
  expect_equivalent(EMMA(gData = gDataTest, trait = 2,
                         environment = 1, K = Sigma, snpName = "M1")[[1]],
                    c(0.184076137051078, 1.83600897652493))
  expect_equivalent(EMMA(gData = gDataTest, trait = 2,
                         environment = 1, covar = 1, snpName = "M1")[[1]],
                    c(0.140370911511506, 1.71006778174736))
  expect_equivalent(EMMA(gData = gDataTest, trait = 2, environment = 1,
                         K = Sigma, covar = 1, snpName = "M1")[[1]],
                    c(0.140370911511506, 1.71006778174736))
})

test_that("extra options in EMMA don't significantly change results", {
  ## Compute base value
  EMMA0 <- EMMA(gData = gDataTest, trait = 2, environment = 1)[[1]]
  expect_equal(EMMA(gData = gDataTest, trait = 2, environment = 1,
                    nGrids = 50)[[1]], EMMA0, tolerance = 1e-6)
  expect_equal(EMMA(gData = gDataTest, trait = 2, environment = 1,
                    nGrids = 500)[[1]], EMMA0, tolerance = 1e-6)
  expect_equal(EMMA(gData = gDataTest, trait = 2, environment = 1,
                    lLim = -100)[[1]], EMMA0, tolerance = 1e-6)
  expect_equal(EMMA(gData = gDataTest, trait = 2, environment = 1,
                    lLim = -20)[[1]], EMMA0, tolerance = 1e-6)
  expect_equal(EMMA(gData = gDataTest, trait = 2, environment = 1,
                    uLim = 20)[[1]], EMMA0, tolerance = 1e-6)
  expect_equal(EMMA(gData = gDataTest, trait = 2, environment = 1,
                    uLim = 100)[[1]], EMMA0, tolerance = 1e-6)
  expect_equal(EMMA(gData = gDataTest, trait = 2, environment = 1,
                    eps = 1e-5)[[1]], EMMA0, tolerance = 1e-6)
  expect_equal(EMMA(gData = gDataTest, trait = 2, environment = 1,
                    nGrids = 500, lLim = -100, uLim = 100)[[1]],
               EMMA0, tolerance = 1e-6)
})

test_that("EMMAREMLLL produces correct output", {
  expect_equivalent(emmaREMLLL(logDelta = -1, lambda = 1, etas1 = 2, n = 0,
                               t = 0, etas2 = 0), -2.11208571376462)
  expect_equivalent(emmaREMLLL(logDelta = -1, lambda = 1:10, etas1 = 1:10,
                               n = 0, t = 0, etas2 = 0), -30.4461091155684)
  expect_equivalent(emmaREMLLL(logDelta = -1, lambda = 1:10, etas1 = 1:10,
                               n = 3, t = 2, etas2 = 1:2),
                    c(-31.9439236241564, -32.2122231996107))
})

GLS0 <- fastGLS(y = y, X = X, Sigma = Sigma)
test_that("fastGLS produces correct output structure", {
  expect_is(GLS0, "data.frame")
  expect_equal(dim(GLS0), c(3, 4))
  expect_equal(rownames(GLS0), paste0("M", 1:3))
  expect_equal(colnames(GLS0), c("pValue", "beta", "betaSe", "RLR2"))
})

test_that("fastGLS without covariates produces correct output", {
  expect_equal(GLS0[, 1],
               c(0.191990244479038, 0.0346367487131218, 0.297099155797429))
  expect_equal(GLS0[, 2],
               c(-0.765513606856416, -2.9009783711453, 0.760762479686889))
  expect_equal(GLS0[, 3],
               c(0.125433669816756, 0.319994381086895, 0.152892637284377))
  expect_equal(GLS0[, 4],
               c(0.975876824782575, 0.999730440588029, 0.91590885222019))
})

test_that("fastGLS with covariates produces correct output", {
  GLS1 <- fastGLS(y = y, X = X, Sigma = Sigma, covs = covs)
  expect_equal(GLS1[, 1],
               c(0.729670715779269, 0.0632229836715346, 0.489762145590089))
  expect_equal(GLS1[, 2],
               c(0.779388868965725, -3.50930856822299, 1.51419592119873))
  expect_equal(GLS1[, 3],
               c(0.483570573233199, 0.467906979316023, 0.47775066615577))
  expect_equal(GLS1[, 4],
               c(0.228770901452724, 0.996393508816207, 0.633782123371107))
})

test_that("fastGLS is independent of dimensions", {
  expect_equal(fastGLS(y = y, X = X[, 1:2], Sigma = Sigma),
               fastGLS(y = y, X = X, Sigma = Sigma)[1:2, ])
})

test_that("function exclMarkers functions properly", {
  markers <- matrix(c(0, 1, 0, 1, 2, 1, 0, 1, 0, 2, 1, 2), ncol = 4,
                    dimnames = list(paste0("IND", 1:3), paste0("SNP", 1:4)))
  allFreq <- colMeans(markers, na.rm = TRUE)
  expect_length(exclMarkers(snpCov = NULL, markers = markers,
                            allFreq = allFreq), 0)
  expect_equal(exclMarkers(snpCov = "SNP2", markers = markers,
                           allFreq = allFreq), 2)
  expect_equal(exclMarkers(snpCov = "SNP1", markers = markers,
                           allFreq = allFreq), c(1, 3))
  expect_equal(exclMarkers(snpCov = c("SNP1", "SNP3"), markers = markers,
                           allFreq = allFreq), c(1, 3))
  expect_equal(exclMarkers(snpCov = c("SNP1", "SNP2", "SNP3"),
                           markers = markers, allFreq = allFreq), c(1, 3, 2))
})

test_that("genCtrlPVals produces correct new p-values", {
  expect_equal(genCtrlPVals(pVals = .5, nObs = 10)[[1]], 0.5)
  expect_equal(genCtrlPVals(pVals = c(0.25, 0.5), nObs = 10)[[1]],
               c(0.410638105484779, 0.63432328826532))
  expect_equal(genCtrlPVals(pVals = c(0.25, 0.5), nObs = 1e6)[[1]],
               c(0.410588927034021, 0.629472060364479))
})

test_that("genCtrlPVals produces correct inflation factor", {
  expect_equal(genCtrlPVals(pVals = .5, nObs = 10)[[2]], 1)
  expect_equal(genCtrlPVals(pVals = c(0.25, 0.5), nObs = 10)[[2]],
               2.04152779518634)
  expect_equal(genCtrlPVals(pVals = c(0.25, 0.5), nObs = 1e6)[[2]],
               1.95438310682772)
})

