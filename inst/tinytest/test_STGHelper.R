load(file = "testdata.rda")

### Test GWAS helper functions.

## Test EMMA function

# Check with default settings.

expect_equivalent(statgenGWAS:::EMMA(gData = gDataTest, trait = 2, 
                                     trial = 1)[[1]],
                  c(0.020597492367456, 1.85412717490278))
expect_equivalent(statgenGWAS:::EMMA(gData = gDataTest, trait = "X1",
                                     trial = 1)[[1]],
                  c(0.020597492367456, 1.85412717490278))
expect_equivalent(statgenGWAS:::EMMA(gData = gDataTest, trait = "X1",
                                     trial = "ph1")[[1]],
                  c(0.020597492367456, 1.85412717490278))
expect_equivalent(statgenGWAS:::EMMA(gData = gDataTest, trait = 2,
                                     trial = 1, K = Sigma)[[1]],
                  c(0.020597492367456, 1.85412717490278))

# Check with covariates.
# Use combinations of covariates and SNP covariates.

expect_equivalent(statgenGWAS:::EMMA(gData = gDataTest, trait = 2,
                                     trial = 1, covar = 1)[[1]],
                  c(8.76962021955837e-05, 1.93163739799548))
expect_equivalent(statgenGWAS:::EMMA(gData = gDataTest, trait = 2,
                                     trial = 1, covar = "V1")[[1]],
                  c(8.76962021955837e-05, 1.93163739799548))
expect_equivalent(statgenGWAS:::EMMA(gData = gDataTest, trait = 2,
                                     trial = 1, K = Sigma, covar = 1)[[1]],
                  c(8.76962021955837e-05, 1.93163739799548))
expect_equivalent(statgenGWAS:::EMMA(gData = gDataTest, trait = 2,
                                     trial = 1, snpName = "M1")[[1]],
                  c(9.11116376851807e-05, 2.00686737098146))
expect_equivalent(statgenGWAS:::EMMA(gData = gDataTest, trait = 2,
                                     trial = 1, K = Sigma, snpName = "M1")[[1]],
                  c(9.11116376851807e-05, 2.00686737098146))
expect_equivalent(statgenGWAS:::EMMA(gData = gDataTest, trait = 2,
                                     trial = 1, covar = 1, snpName = "M1")[[1]],
                  c(8.77313241082566e-05, 1.93241100960362))
expect_equivalent(statgenGWAS:::EMMA(gData = gDataTest, trait = 2, trial = 1,
                                     K = Sigma, covar = 1, snpName = "M1")[[1]],
                  c(8.77313241082566e-05, 1.93241100960362))

# Check that extra options don't significantly change results.
# Changing the number of grid points or setting different upper or lower
# limits should give output very similar to the default settings.

## Compute base value
EMMA0 <- statgenGWAS:::EMMA(gData = gDataTest, trait = 2, trial = 1)[[1]]
expect_equal(statgenGWAS:::EMMA(gData = gDataTest, trait = 2, trial = 1,
                                nGrids = 50)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(gData = gDataTest, trait = 2, trial = 1,
                                nGrids = 500)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(gData = gDataTest, trait = 2, trial = 1,
                                lLim = -100)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(gData = gDataTest, trait = 2, trial = 1,
                                lLim = -20)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(gData = gDataTest, trait = 2, trial = 1,
                                uLim = 20)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(gData = gDataTest, trait = 2, trial = 1,
                                uLim = 100)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(gData = gDataTest, trait = 2, trial = 1,
                                eps = 1e-5)[[1]], EMMA0, tolerance = 1e-6)
expect_equal(statgenGWAS:::EMMA(gData = gDataTest, trait = 2, trial = 1,
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
expect_true(inherits(GLS0, "matrix"))
expect_equal(dim(GLS0), c(3, 4))
expect_null(rownames(GLS0))
expect_null(colnames(GLS0))

# Check without covariates.

expect_equal(GLS0, 
             matrix(c(0.0271120646525595, 0.556661238666598, 0.203405705505351, 
                      -2.90711358246727, -1.01229725829335, -1.73072077802875, 
                      0.310442837982955, 0.352135174778732, 0.290108955585292,
                      0.999844554418642, 0.562383932125061, 0.971533598465952),
                    nrow = 3))

# Check with covariates.

GLS1 <- statgenGWAS:::fastGLS(y = y, X = X, Sigma = Sigma, covs = covs)
expect_equal(GLS1, 
             matrix(c(0.158456432903659, 0.592448089085206, 0.861012555940774, 
                      -4.0813171867358, -1.23573112888003, 0.547427202536095, 
                      0.67415615801093, 0.498746367260275, 0.667641616295594, 
                      0.974397061281524, 0.458757080015615, 0.0650202741814823),
                    nrow = 3))

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

