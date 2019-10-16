load(file = "testdata.rda")

### Test check functions.

## Test chkGData.
expect_error(statgenGWAS:::chkGData(), "valid gData object")
expect_error(statgenGWAS:::chkGData(X), "valid gData object")
expect_silent(statgenGWAS:::chkGData(gDataTest))

# Copy and remove markers to create an invalid gData object.
gDataTest2 <- gDataTest
gDataTest2$markers <- NULL
expect_error(statgenGWAS:::chkGData(gDataTest2), "should contain markers")

# Map is still present so if only map is needed check should succeed.
expect_silent(statgenGWAS:::chkGData(gDataTest2, comps = "map"))

## Test chkTrials.

# Input should be character or numeric.
expect_error(statgenGWAS:::chkTrials(trials = TRUE, gData = gDataTest), 
             "should be a numeric or character vector")

# Check character input.
expect_silent(statgenGWAS:::chkTrials(trials = "ph1", gData = gDataTest))
expect_silent(statgenGWAS:::chkTrials(trials = c("ph1", "ph2"), 
                                      gData = gDataTest))
expect_error(statgenGWAS:::chkTrials(trials = "ph3", gData = gDataTest),
             "should be in pheno")

# Check numeric input.
expect_silent(statgenGWAS:::chkTrials(trials = 1, gData = gDataTest))
expect_silent(statgenGWAS:::chkTrials(trials = c(1, 2), gData = gDataTest))
expect_error(statgenGWAS:::chkTrials(trials = 3, gData = gDataTest), 
             "should be in pheno")

## Test chkTraits.

# Input should be character or numeric.
expect_error(statgenGWAS:::chkTraits(traits = TRUE, trials = "ph1", 
                                     gData = gDataTest, multi = TRUE), 
             "should be a numeric or character vector")

# When multi = FALSE only single numeric/character values are accepted.
expect_error(statgenGWAS:::chkTraits(traits = 1:2, trials = "ph1", 
                                     gData = gDataTest, multi = FALSE), 
             "should be a single numeric or character value")

# Check character input.
expect_silent(statgenGWAS:::chkTraits(traits = "X1", trials = "ph1", 
                                      gData = gDataTest, multi = FALSE))
expect_silent(statgenGWAS:::chkTraits(traits = c("X1", "X2"), trials = "ph1", 
                                      gData = gDataTest, multi = TRUE))
expect_silent(statgenGWAS:::chkTraits(traits = c("X1", "X2"), 
                                      trials = c("ph1", "ph2"), 
                                      gData = gDataTest, multi = TRUE))
expect_error(statgenGWAS:::chkTraits(traits = "X6", trials = "ph1", 
                                     gData = gDataTest, multi = FALSE), 
             "For ph1 not all traits")
expect_error(statgenGWAS:::chkTraits(traits = c("X1", "X6"), trials = "ph1", 
                                     gData = gDataTest, multi = TRUE),
             "For ph1 not all traits")

# Check numeric input.
expect_silent(statgenGWAS:::chkTraits(traits = 2, trials = 1, gData = gDataTest, 
                                      multi = FALSE))
expect_silent(statgenGWAS:::chkTraits(traits = c(2, 3), trials = 1, 
                                      gData = gDataTest, multi = TRUE))
expect_silent(statgenGWAS:::chkTraits(traits = c(2, 3), trials = c(1, 2), 
                                      gData = gDataTest, multi = TRUE))
expect_error(statgenGWAS:::chkTraits(traits = 1, trials = 1, gData = gDataTest, 
                                     multi = FALSE), 
             "For 1 not all traits")
expect_error(statgenGWAS:::chkTraits(traits = c(2, 7), trials = 1, 
                                     gData = gDataTest, multi = TRUE), 
             "For 1 not all traits")

## Test chkNum

expect_silent(statgenGWAS:::chkNum(1, min = 0, max = 2))
expect_error(statgenGWAS:::chkNum("a"), "single numerical value")
expect_error(statgenGWAS:::chkNum(c(1, 2)), "single numerical value")
expect_error(statgenGWAS:::chkNum(1, min = NULL, max = 0), "smaller than 0")
expect_error(statgenGWAS:::chkNum(1, min = 2, max = NULL), "greater than 2")
expect_error(statgenGWAS:::chkNum(1, min = 2, max = 3), "between 2 and 3")

## Test chkMarkers

# Input should be numerical.
testMrk <- matrix(letters[1:4], nrow = 2)
expect_error(statgenGWAS:::chkMarkers(markers = testMrk), 
             "should be a numerical matrix")

# No missing values allowed.
testMrk2 <- matrix(c(1:3, NA), nrow = 2)
expect_error(statgenGWAS:::chkMarkers(markers = testMrk2), 
             "markers contains missing values")

expect_error(statgenGWAS:::chkMarkers(markers = gDataTest$markers, dim = 3),
             "markers should be a three-dimensional array")

expect_silent(statgenGWAS:::chkMarkers(markers = gDataTest$markers))

## Test chkCovar

# Input should be character or numeric.
expect_error(statgenGWAS:::chkCovar(covar = TRUE, gData = gDataTest), 
             "should be a numeric or character vector")

# Covariates should be in covar in gData.
expect_error(statgenGWAS:::chkCovar(covar = "V3", gData = gDataTest),
             "covar should be columns in covar in gData")

expect_silent(statgenGWAS:::chkCovar(covar = "V1", gData = gDataTest))

## Test chkSnpCov

# SNP Covariates should be in markers in gData.
expect_error(statgenGWAS:::chkSnpCov(snpCov = "SNP1", gData = gDataTest),
             "All snpCov should be in markers")

expect_silent(statgenGWAS:::chkSnpCov(snpCov = "M1", gData = gDataTest))

## Test chkKin

# When GLSMethod = "single" kin should be a matrix.
expect_error(statgenGWAS:::chkKin(kin = 1, gData = gDataTest,
                                  GLSMethod = "single"),
             "kin should be a matrix")
expect_silent(statgenGWAS:::chkKin(kin = gDataTest$kinship, gData = gDataTest,
                                   GLSMethod = "single"))

# When GLSMethod = "multi" kin should be a list of matrices.
expect_error(statgenGWAS:::chkKin(kin = gDataTest$kinship, gData = gDataTest,
                                   GLSMethod = "multi"),
             "kin should be a list of matrices of length equal to the number")
expect_silent(statgenGWAS:::chkKin(kin = list(gDataTest$kinship,
                                              gDataTest$kinship),
                                   gData = gDataTest, GLSMethod = "multi"))

