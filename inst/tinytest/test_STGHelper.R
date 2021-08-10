load(file = "testdata.rda")

### Test GWAS helper functions.

## Test EMMA function

pheno <- gDataTest$pheno$ph1
covar <- gDataTest$covar

# Check with default settings.

expect_equivalent(statgenGWAS:::EMMA(dat = pheno, trait = 2, K = Sigma)[[1]],
                  c(9.80733999117857e-05, 2.16021038853735))
expect_equivalent(statgenGWAS:::EMMA(dat = pheno, trait = "X1", K = Sigma)[[1]],
                  c(9.80733999117857e-05, 2.16021038853735))

# Check with covariates.
# Use combinations of covariates and SNP covariates.

expect_equivalent(statgenGWAS:::EMMA(dat = pheno, trait = 2, K = Sigma,
                                     covar = covar[, 1, drop = FALSE])[[1]],
                  c(0.000107734030018285, 2.37299992713443))

# Check that extra options don't significantly change results.
# Changing the number of grid points or setting different upper or lower
# limits should give output very similar to the default settings.

## Compute base value
EMMA0 <- statgenGWAS:::EMMA(dat = pheno, trait = 2, K = Sigma)[[1]]
expect_equal(statgenGWAS:::EMMA(dat = pheno, trait = 2, K = Sigma,
                                nGrids = 50)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(dat = pheno, trait = 2, K = Sigma, 
                                nGrids = 500)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(dat = pheno, trait = 2, K = Sigma, 
                                lLim = -100)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(dat = pheno, trait = 2, K = Sigma, 
                                lLim = -20)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(dat = pheno, trait = 2, K = Sigma, 
                                uLim = 20)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(dat = pheno, trait = 2, K = Sigma, 
                                uLim = 100)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(dat = pheno, trait = 2, K = Sigma, 
                                eps = 1e-5)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(dat = pheno, trait = 2, K = Sigma, 
                                nGrids = 500, lLim = -100, uLim = 100)[[1]],
             EMMA0, tolerance = 1e-6)

## Test emmaREMLLL

# All single input values.
expect_equivalent(statgenGWAS:::emmaREMLLL(logDelta = -1, lambda = 1, etas1 = 2,
                                           n = 0, t = 0, etas2 = 0), 
                  matrix(-2.11208571376462))

# Vector for etas1 should still give single result.
expect_equivalent(statgenGWAS:::emmaREMLLL(logDelta = -1, lambda = 1:10, 
                                           etas1 = 1:10, n = 0, t = 0, 
                                           etas2 = 0), 
                  matrix(-30.4461091155684))

# Vector for etas2 should give vector of the same length as result.
expect_equivalent(statgenGWAS:::emmaREMLLL(logDelta = -1, lambda = 1:10, 
                                           etas1 = 1:10, n = 3, t = 2, 
                                           etas2 = 1:2),
                  matrix(c(-31.9439236241564, -32.2122231996107)))

## Test fastGLS

# Check output format.

GLS0 <- statgenGWAS:::fastGLS(y = y, X = X, Sigma = Sigma)
expect_true(inherits(GLS0, "data.table"))
expect_equal(dim(GLS0), c(4, 5))
expect_equal(colnames(GLS0), c("rn", "pValue", "effect", "effectSe", "RLR2"))

# Check without covariates.

expect_equivalent(
  as.matrix(GLS0[, 2:5]), 
  matrix(c(0.00902917885035805, 0.348617547569426, 0.0943676219534976, 
           0.0128620260615326, -16.7349901664247, -9.90568316516165, 
           -13.2646467301295, 8.80467920925613, 0.249876050144977, 
           0.343464035166123, 0.274126541008359, 0.135523209809584, 
           1, 1, 1, 1), nrow = 4))

# Check with covariates.

GLS1 <- statgenGWAS:::fastGLS(y = y, X = X, Sigma = Sigma, covs = covs)
expect_equivalent(
  as.matrix(GLS1[, 2:5]), 
  matrix(c(0.762937476428709, 0.983140796459948, 0.109726700344376, 
           0.274137725349961, 1.69578638166878, 0.167226553742237, 
           6.47091225567247, -3.98457205198297, 0.61037456571308, 
           0.855589996753001, 0.489578911347709, 0.415729973147865, 
           0.537856284290155, 0.00381285082817118, 0.999999974117709, 
           0.99989756698506), nrow = 4))

# Check that fastGLS is independent of dimensions.

expect_equal(statgenGWAS:::fastGLS(y = y, X = X[, 1:2], Sigma = Sigma),
             statgenGWAS:::fastGLS(y = y, X = X, Sigma = Sigma)[1:2, ])

## Test exclMarkers

markers <- matrix(c(0, 1, 0, 1, 2, 1, 0, 1, 0, 2, 1, 2), ncol = 4,
                  dimnames = list(paste0("IND", 1:3), paste0("SNP", 1:4)))
allFreq <- colMeans(markers, na.rm = TRUE)
expect_equal(statgenGWAS:::exclMarkers(snpCov = NULL, markers = markers,
                                       allFreq = allFreq), integer(0))
expect_equal(statgenGWAS:::exclMarkers(snpCov = "SNP2", markers = markers,
                                       allFreq = allFreq), 2)
expect_equal(statgenGWAS:::exclMarkers(snpCov = "SNP1", markers = markers,
                                       allFreq = allFreq), c(1, 3))
expect_equal(statgenGWAS:::exclMarkers(snpCov = c("SNP1", "SNP3"), 
                                       markers = markers, allFreq = allFreq), 
             c(1, 3))
expect_equal(statgenGWAS:::exclMarkers(snpCov = c("SNP1", "SNP2", "SNP3"),
                                       markers = markers, allFreq = allFreq), 
             c(1, 3, 2))

## Test genCtrlPVals 

# Check p-Values.

expect_equal(statgenGWAS:::genCtrlPVals(pVals = .5, nObs = 10)[[1]], 0.5)
expect_equal(statgenGWAS:::genCtrlPVals(pVals = c(0.25, 0.5), nObs = 10)[[1]],
             c(0.410638105484779, 0.63432328826532))
expect_equal(statgenGWAS:::genCtrlPVals(pVals = c(0.25, 0.5), nObs = 1e6)[[1]],
             c(0.410588927034021, 0.629472060364479))

# Check inflation factor.

expect_equal(statgenGWAS:::genCtrlPVals(pVals = .5, nObs = 10)[[2]], 1)
expect_equal(statgenGWAS:::genCtrlPVals(pVals = c(0.25, 0.5), nObs = 10)[[2]],
             2.04152779518634)
expect_equal(statgenGWAS:::genCtrlPVals(pVals = c(0.25, 0.5), nObs = 1e6)[[2]],
             1.95438310682772)

