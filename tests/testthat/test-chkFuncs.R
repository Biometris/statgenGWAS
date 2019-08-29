context("check Functions")

test_that("chkGData functions properly", {
  expect_error(chkGData(), "valid gData object")
  expect_error(chkGData(X), "valid gData object")
  expect_silent(chkGData(gDataTest))
  gDataTest2 <- gDataTest
  gDataTest2$markers <- NULL
  expect_error(chkGData(gDataTest2), "should contain markers")
  expect_silent(chkGData(gDataTest2, comps = "map"))
})

test_that("chkTrials functions properly", {
  expect_silent(chkTrials("ph1", gDataTest))
  expect_silent(chkTrials(c("ph1", "ph2"), gDataTest))
  expect_error(chkTrials("ph3", gDataTest), "should be in pheno")
  expect_silent(chkTrials(1, gDataTest))
  expect_silent(chkTrials(c(1, 2), gDataTest))
  expect_error(chkTrials(3, gDataTest), "should be in pheno")
})

test_that("chkTraits functions properly", {
  expect_silent(chkTraits("X1", "ph1", gDataTest, FALSE))
  expect_silent(chkTraits(c("X1", "X2"), "ph1", gDataTest, TRUE))
  expect_silent(chkTraits(c("X1", "X2"), c("ph1", "ph2"), gDataTest, TRUE))
  expect_error(chkTraits("X6", "ph1", gDataTest, FALSE), 
               "For ph1 not all traits")
  expect_error(chkTraits(c("X1", "X6"), "ph1", gDataTest, TRUE),
               "For ph1 not all traits")
  expect_silent(chkTraits(2, 1, gDataTest, FALSE))
  expect_silent(chkTraits(c(2, 3), 1, gDataTest, TRUE))
  expect_silent(chkTraits(c(2, 3), c(1, 2), gDataTest, TRUE))
  expect_error(chkTraits(1, 1, gDataTest, FALSE), "For 1 not all traits")
  expect_error(chkTraits(c(2, 7), 1, gDataTest, TRUE), "For 1 not all traits")
})

test_that("chkNum functions properly", {
  expect_silent(chkNum(1, min = 0, max = 2))
  expect_error(chkNum("a"), "single numerical value")
  expect_error(chkNum(c(1, 2)), "single numerical value")
  expect_error(chkNum(1, min = NULL, max = 0), "smaller than 0")
  expect_error(chkNum(1, min = 2, max = NULL), "greater than 2")
  expect_error(chkNum(1, min = 2, max = 3), "between 2 and 3")
})
