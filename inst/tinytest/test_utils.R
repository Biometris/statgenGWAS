load(file = "testdata.rda")

### Test utils

## Test dfBind

# Check that columns are copied correctly.
df1 <- data.frame(a = 1:2, b = 1:2)
df2 <- data.frame(a = 1:2, c = 1:2)
df3 <- data.frame(c = 1:2, d = 1:2)
expect_equal(colnames(statgenGWAS:::dfBind(list(df1, df1))), c("a", "b"))
expect_equal(colnames(statgenGWAS:::dfBind(list(df1, df2))), c("a", "b", "c"))
expect_equal(colnames(statgenGWAS:::dfBind(list(df1, df3))), 
             c("a", "b", "c", "d"))
expect_equal(colnames(statgenGWAS:::dfBind(list(df1, df2, df3))), 
             c("a", "b", "c", "d"))

# Check that NAs are inserted for missing columns.
expect_equivalent(unlist(statgenGWAS:::dfBind(list(df1, df2))),
                  c(1, 2, 1, 2, 1, 2, NA, NA, NA, NA, 1, 2))
expect_equivalent(unlist(statgenGWAS:::dfBind(list(df1, df2, df1))),
                  c(1, 2, 1, 2, 1,2, 1, 2, NA, NA, 1, 2, NA, NA, 1, 2, NA, NA))

# Check that empty data.frames are removed before binding.
expect_equal(statgenGWAS:::dfBind(list(data.frame(), df1)), df1)
expect_equal(statgenGWAS:::dfBind(list(df1, data.frame())), df1)
expect_equal(statgenGWAS:::dfBind(list(data.frame())), data.frame())

## Test matrixRoot

M1 <- matrix(1:4, nrow = 2)
M2 <- matrix(c(1:2, 2:1), nrow = 2)
expect_error(statgenGWAS:::matrixRoot(M1), 
             "should be a symmetric positive definite matrix")
expect_error(statgenGWAS:::matrixRoot(M2), 
             "should be a symmetric positive definite matrix")
expect_equal(statgenGWAS:::matrixRoot(crossprod(M2)),
             matrix(c(2, 1, 1, 2), nrow = 2))

## Test reduceKinship

M3 <- statgenGWAS:::reduceKinship(M2, nPca = 1)
expect_true(inherits(M3, "matrix"))
expect_equal(dim(M3), dim(M2))
expect_equal(M3, matrix(rep(1.5, times = 4), nrow = 2))

## Test nearestPD

M4 <- statgenGWAS:::nearestPD(M1)
expect_true(inherits(M4, "matrix"))
expect_true(isSymmetric(M4))
expect_equal(M4, matrix(c(1.31461827942692, 2.32186609398592, 2.32186609398592, 
                          4.10085766799573), nrow = 2))

# Test parameters in nearestPD

M5 <- statgenGWAS:::nearestPD(M1, corr = TRUE)
expect_equal(M5, matrix(c(1, 0.99999998, 0.99999998, 1), nrow = 2))

M6 <- statgenGWAS:::nearestPD(M1, keepDiag = TRUE)
expect_equal(M6, matrix(c(1, 1.99999993750003, 1.99999993750003, 4), nrow = 2))

M7 <- statgenGWAS:::nearestPD(M1, do2eigen = FALSE)
expect_equal(M7, matrix(c(1.31461827942692, 2.32186609398592, 2.32186609398592, 
                          4.10085766799573), nrow = 2))

M8 <- statgenGWAS:::nearestPD(M1, doDykstra = FALSE)
expect_equal(M8, matrix(c(1.31461827942692, 2.32186609398592, 2.32186609398592, 
                          4.10085766799573), nrow = 2))

M9 <- statgenGWAS:::nearestPD(M1, doSym = TRUE)
expect_equal(M9, matrix(c(1.31461827942692, 2.32186609398592, 2.32186609398592, 
                          4.10085766799573), nrow = 2))
 
## Test computeKin

# Test for GLSMethod single.

K0 <- Sigma + 0.1
gDataTestK0 <- createGData(kin = K0)

# Only kin provided -> return directly.
K1 <- statgenGWAS:::computeKin(GLSMethod = "single", kin = K0)
expect_true(inherits(K1, "matrix"))
expect_equal(K1, K0)

# Only gData provided -> return directly.
K2 <- statgenGWAS:::computeKin(GLSMethod = "single", gData = gDataTestK0)
expect_true(inherits(K2, "matrix"))
expect_equal(K2, K0)

# Both kin and gData provided -> Return kin.
expect_equal(statgenGWAS:::computeKin(GLSMethod = "single", kin = K0, 
                                      gData = gDataTestK0), K0)

# Test for GLSMethod multi.

K0 = Sigma + 0.1
gDataTestK0M <- createGData(kin = list("chr1" = K0, "chr2" = K0))

# Only kin provided -> return directly.
KLst1 <- statgenGWAS:::computeKin(GLSMethod = "multi",
                                  kin = list("chr1" = K0, "chr2" = K0))
expect_true(inherits(KLst1, "list"))
expect_true(inherits(KLst1[[1]], "matrix"))
expect_equal(KLst1[[1]], K0)

# Only gData provided -> return directly.
KLst2 <- statgenGWAS:::computeKin(GLSMethod = "multi", gData = gDataTestK0M)
expect_true(inherits(KLst2, "list"))
expect_true(inherits(KLst2[[1]], "matrix"))
expect_equal(KLst2[[1]], K0)

# Both kin and gData provided -> return directly.
expect_equal(statgenGWAS:::computeKin(GLSMethod = "multi", 
                                      kin = list("chr1" = K0, "chr2" = K0), 
                                      gData = gDataTestK0M), KLst1)

# Test for correct output when kinship matrix is actually computed.

# Test for GLSMethod single.
K3 <- statgenGWAS:::computeKin(GLSMethod = "single", 
                               markers = gDataTest$markers)

expect_true(inherits(K3, "matrix"))
expect_equal(dim(K3), c(10, 10))
expect_equivalent(K3[1:2, 1:2], 
                  matrix(rep(0.161495911495911, times = 4), nrow = 2))

# Test for GLSMethod multi.

# multi is only possible if there are multiple chromosomes in the map.
expect_error(statgenGWAS:::computeKin(GLSMethod = "multi", 
                                      markers = gDataTest$markers[, 1:2], 
                                      map = gDataTest$map[1:2, ]),
             "Chromosome specific kinship calculation not possible")

KLst3 <- statgenGWAS:::computeKin(GLSMethod = "multi", 
                                  markers = gDataTest$markers, 
                                  map = gDataTest$map)

expect_true(inherits(KLst3, "list"))
expect_true(inherits(KLst3[[1]], "matrix"))
expect_equal(dim(KLst3[[1]]), c(10, 10))
expect_equivalent(KLst3[[1]][1:2, 1:2], 
                  matrix(rep(0.0202020202020202, times = 4), nrow = 2))

## Test expand pheno.

# Check for covar and snpCov NULL.
expPh1 <- statgenGWAS:::expandPheno(gData = gDataTest, trial = "ph1", 
                                    covar = NULL, snpCov = NULL)
expect_equal(expPh1$phTr, gDataTest$pheno$ph1)
expect_null(expPh1$covTr)

# Check for covar not NULL and snpCov NULL.

# Add factor covariate.
gDataTest$covar$V3 <- factor(rep(1:5, times = 2))

# Single covariate.
expPh2 <- statgenGWAS:::expandPheno(gData = gDataTest, trial = "ph1", 
                                    covar = "V1", snpCov = NULL)
expect_true(inherits(expPh2$phTr, "data.frame"))
expect_equal(colnames(expPh2$phTr), c(colnames(gDataTest$pheno$ph1), "V1"))
expect_equal(expPh2$covTr, "V1")

# Multiple covariates.
expPh3 <- statgenGWAS:::expandPheno(gData = gDataTest, trial = "ph1", 
                                    covar = c("V1", "V2"), snpCov = NULL)
expect_equal(colnames(expPh3$phTr), 
             c(colnames(gDataTest$pheno$ph1), c("V1", "V2")))
expect_equal(expPh3$covTr, c("V1", "V2"))

# Single factor covariate.
expPh4 <- statgenGWAS:::expandPheno(gData = gDataTest, trial = "ph1", 
                                    covar = "V3", snpCov = NULL)
expect_equal(colnames(expPh4$phTr), 
             c(colnames(gDataTest$pheno$ph1), "V32", "V33", "V34", "V35"))
expect_equal(expPh4$covTr, c("V32", "V33", "V34", "V35"))

# Check for covar NULL and snpCov not NULL.

# Single SNP covariate.

expPh5 <- statgenGWAS:::expandPheno(gData = gDataTest, trial = "ph1", 
                                    covar = NULL, snpCov = "M1")
expect_true(inherits(expPh5$phTr, "data.frame"))
expect_equal(colnames(expPh5$phTr), c(colnames(gDataTest$pheno$ph1), "M1"))
expect_equal(expPh5$covTr, "M1")

# Multiple SNP covariates.

expPh6 <- statgenGWAS:::expandPheno(gData = gDataTest, trial = "ph1", 
                                    covar = NULL, snpCov = c("M1", "M2"))
expect_equal(colnames(expPh6$phTr), c(colnames(gDataTest$pheno$ph1), 
                                      "M1", "M2"))
expect_equal(expPh6$covTr, c("M1", "M2"))

# Check for both covar and snpCov not NULL.

# Multiple covariates and SNP covariates.

expPh7 <- statgenGWAS:::expandPheno(gData = gDataTest, trial = "ph1", 
                                    covar = c("V1", "V3"), 
                                    snpCov = c("M1", "M3"))
expect_true(inherits(expPh5$phTr, "data.frame"))
expect_equal(colnames(expPh7$phTr), 
             c(colnames(gDataTest$pheno$ph1), 
               "V1", "V32", "V33", "V34", "V35", "M1", "M3"))
expect_equal(expPh7$covTr, c("V1", "V32", "V33", "V34", "V35", "M1", "M3"))


