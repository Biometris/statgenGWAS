context("kinship functions")

test_that("kinship functions give correct output", {
  X <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_equivalent(astleCPP(X),c(2/3, -2/3, -2/3, 2/3))
  expect_equivalent(IBSCPP(X), c(1, 0, 0, 1))
  expect_equivalent(vanRadenCPP(X), c(2/3, -2/3, -2/3, 2/3))
})

test_that(paste("kinship functions give correct output with user",
                "definded denominator"), {
  X <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_equivalent(astleCPP(X, denom = 4), c(1/3, -1/3, -1/3, 1/3))
  expect_equivalent(IBSCPP(X, denom = 4), c(0.5, 0, 0, 0.5))
  expect_equivalent(vanRadenCPP(X, denom = 4), c(1/8, -1/8, -1/8, 1/8))
})

test_that("function kinship functions properly for 2d markers", {
  X <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_equal(kinship(X = X, method = "astle"), astleCPP(X))
  expect_equal(kinship(X = X, method = "IBS", denominator = 2),
               IBSCPP(X, denom = 2))
  expect_equal(kinship(X = X, method = "vanRaden", denominator = 3),
               vanRadenCPP(X, denom = 3))
})


