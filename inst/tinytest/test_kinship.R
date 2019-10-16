### Test kinship functions.

## Test with default setting.
X <- matrix(c(1, 0, 0, 1), nrow = 2)
expect_equal(statgenGWAS:::astleCPP(X), 
             matrix(c(2/3, -2/3, -2/3, 2/3), nrow = 2))
expect_equal(statgenGWAS:::IBSCPP(X), X)
expect_equal(statgenGWAS:::vanRadenCPP(X), 
             matrix(c(2/3, -2/3, -2/3, 2/3), nrow = 2))

## Test with user defined denominator.
expect_equal(statgenGWAS:::astleCPP(X, denom = 4), 
             matrix(c(1/3, -1/3, -1/3, 1/3), nrow = 2))
expect_equal(statgenGWAS:::IBSCPP(X, denom = 4), 
             matrix(c(0.5, 0, 0, 0.5), nrow = 2))
expect_equal(statgenGWAS:::vanRadenCPP(X, denom = 4), 
             matrix(c(1/8, -1/8, -1/8, 1/8), nrow = 2))

## Test actual kinship functions.
expect_equal(kinship(X = X, method = "astle"), statgenGWAS:::astleCPP(X))
expect_equal(kinship(X = X, method = "IBS", denominator = 2),
             statgenGWAS:::IBSCPP(X, denom = 2))
expect_equal(kinship(X = X, method = "vanRaden", denominator = 3),
             statgenGWAS:::vanRadenCPP(X, denom = 3))



