context("utils")

test_that("function dfBind copies columns properly", {
  df1 <- data.frame(a = 1:2, b = 1:2)
  df2 <- data.frame(a = 1:2, c = 1:2)
  df3 <- data.frame(c = 1:2, d = 1:2)
  expect_equal(colnames(dfBind(list(df1, df1))), c("a", "b"))
  expect_equal(colnames(dfBind(list(df1, df2))), c("a", "b", "c"))
  expect_equal(colnames(dfBind(list(df1, df3))), c("a", "b", "c", "d"))
  expect_equal(colnames(dfBind(list(df1, df2, df3))), c("a", "b", "c", "d"))
})

test_that("function dfBind inserts NAs for missing columns", {
  df1 <- data.frame(a = 1:2, b = 1:2)
  df2 <- data.frame(a = 1:2, c = 1:2)
  expect_equivalent(unlist(dfBind(list(df1, df2))),
                    c(1, 2, 1, 2, 1, 2, NA, NA, NA, NA, 1, 2))
  expect_equivalent(unlist(dfBind(list(df1, df2, df1))),
                    c(1, 2, 1, 2, 1,2, 1, 2, NA, NA, 1, 2, NA, NA, 1, 2, NA, NA))
})

test_that(paste("function dfBind removes empty data.frames lists from",
                "input before binding"), {
                  df1 <- data.frame(a = 1:2, b = 1:2)
                  expect_equal(dfBind(list(data.frame(), df1)), df1)
                  expect_equal(dfBind(list(df1, data.frame())), df1)
                  expect_equal(dfBind(list(data.frame())), data.frame())
                })

test_that("function matrixRoot functions properly", {
  M1 <- matrix(1:4, nrow = 2)
  M2 <- matrix(c(1:2, 2:1), nrow = 2)
  expect_error(matrixRoot(M1), "should be a symmetric positive definite matrix")
  expect_error(matrixRoot(M2), "should be a symmetric positive definite matrix")
  expect_equal(as.numeric(matrixRoot(crossprod(M2))), c(2, 1, 1, 2))
})

test_that("function reduceKinship functions properly", {
  K1 <- reduceKinship(K = K, nPca = 2)
  expect_is(K1, "matrix")
  expect_equal(dim(K1), dim(K))
  expect_known_output(K1, file = test_path("test-redKin1.txt"), print = TRUE)
  K2 <- reduceKinship(K = K, nPca = 5)
  expect_known_output(K2, file = test_path("test-redKin2.txt"), print = TRUE)
})

set.seed(1234)
M <- matrix(data = runif(9), nrow = 3)
## Assure symmetry.
M <- M + t(M)
test_that("function nearestPD functions properly", {
  M1 <- nearestPD(M)
  expect_is(M1, "matrix")
  expect_true(isSymmetric(M1))
  expect_equivalent(M1,
                    c(0.604727709553899, 1.04559967845123, 0.586415044934283,
                      1.04559967845123, 1.82792523675371, 0.890017983604708,
                      0.586415044934283, 0.890017983604709, 1.33494200875219))
})

test_that("parameters in nearestPD function properly", {
  M1 <- nearestPD(M, corr = TRUE)
  expect_equivalent(M1, c(1, 0.984775432717047, 0.680679114763,
                          0.984775432717047, 1, 0.7976616570699, 0.680679114763,
                          0.7976616570699, 1))
  expect_warning(nearestPD(M, keepDiag = TRUE),
                 "did not converge in 100 iterations")
  M2 <- nearestPD(M, keepDiag = TRUE, maxIter = 120)
  expect_equivalent(M2,
                    c(0.227406822610646, 0.614704685973532, 0.414723676348398,
                      0.614704685973532, 1.72183076711372, 0.934827342556341,
                      0.414723676348398, 0.934827342556341, 1.3321675164625))
  M3 <- nearestPD(M, do2eigen = FALSE)
  expect_equivalent(M3,
                    c(0.604727709553899, 1.04559971368566, 0.586415058601556,
                      1.04559971368566, 1.82792523675371, 0.890017984212184,
                      0.586415058601556, 0.890017984212184, 1.33494200875219))
  M4 <- nearestPD(M, doDykstra = FALSE)
  expect_equivalent(M4,
                    c(0.604727709553899, 1.04559967845123, 0.586415044934283,
                      1.04559967845123, 1.82792523675371, 0.890017983604708,
                      0.586415044934283, 0.890017983604708, 1.33494200875219))
  M5 <- nearestPD(M, doSym = TRUE)
  expect_equivalent(M5,
                    c(0.604727709553899, 1.04559967845123, 0.586415044934283,
                      1.04559967845123, 1.82792523675371, 0.890017983604708,
                      0.586415044934283, 0.890017983604708, 1.33494200875219))
})
