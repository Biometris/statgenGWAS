load(file = "testdata.rda")

### Test runSingleTraitGwas with single kinship.

## Test option remlAlgo.

expect_error(runSingleTraitGwas(gData = gDataTest, traits = "X1", trials = 1, 
                                remlAlgo = "a"),
             'should be one of "EMMA", "NR"')

stg1 <- runSingleTraitGwas(gData = gDataTest, traits = "X1", trials = 1, 
                          remlAlgo = "NR")
stgM1 <- runSingleTraitGwas(gData = gDataTest, traits = "X1", trials = 1, 
                          remlAlgo = "NR", GLSMethod = "multi")

expect_equal(stg1$GWASInfo$remlAlgo, "NR")
expect_equal(stg1$GWASInfo$varComp$ph1$X1,
             c(Vg = 0.01816883008257, Ve = 1.8560701334347))

expect_equal(stgM1$GWASInfo$remlAlgo, "NR")
expect_equivalent(unlist(stgM1$GWASInfo$varComp$ph1$X1),
                  c(0, 1.87085447659007, 0, 1.87085447659007))

## Test option sizeInclRegion.

expect_error(runSingleTraitGwas(gData = gDataTest, traits = "X1", trials = 1, 
                                sizeInclRegion = -1),
             "should be a single numerical value greater than 0")

# Only relevant when there are significant SNPs
# Lower threshold to get a significant SNP.
# No difference here for GLSMethod single and multi so only check single.

stg2 <- runSingleTraitGwas(gData = gDataTest, traits = "X1", trials = 1, 
                           thrType = "fixed", LODThr = 0.2,
                           sizeInclRegion = 1, minR2 = 0.01)

sig2 <- stg2$signSnp$ph1

expect_true(inherits(sig2, "data.table"))
expect_equal(nrow(sig2), 2)
expect_equal(as.character(sig2[["snpStatus"]]), 
             c("significant SNP", "within 1 of a significant SNP"))
expect_equal(sig2[["LOD"]] > 0.2, c(TRUE, FALSE))

## Test option useMAF.

expect_error(runSingleTraitGwas(gData = gDataTest, traits = "X1", trials = 1, 
                                useMAF = FALSE, MAC = -1),
             "MAC should be a single numerical value greater than 0")

# No difference here for GLSMethod single and multi so only check single.

stg3 <- runSingleTraitGwas(gData = gDataTest, traits = "X1", trials = 1, 
                           useMAF = FALSE, MAC = 2)
# MAC = 2 causes M3 to be non-segregating.
expect_equal(stg3$GWASInfo$MAF, 0.19999)
expect_equal(stg3$GWAResult$ph1$pValue[3], NA_real_)

## Test option thrType.

# Fixed and bonf are tested implicitely in other tests. 
# Only test small here.

expect_error(runSingleTraitGwas(gData = gDataTest, traits = "X1", trials = 1, 
                                thrType = "small", nSnpLOD = -1),
             "nSnpLOD should be a single numerical value greater than 0")

stg4 <- runSingleTraitGwas(gData = gDataTest, traits = "X1", trials = 1, 
                           thrType = "small", nSnpLOD = 1)
expect_equivalent(stg4$thr$ph1, 0.254604498474883)
expect_equal(nrow(stg4$signSnp$ph1), 1)

## Test option genomicControl.

stg5 <- runSingleTraitGwas(gData = gDataTest, traits = "X1", trials = 1, 
                           genomicControl = TRUE)

# Should only affect pValue and LOD
expect_equal(stg5$GWAResult$ph1$pValue, 
             c(0.392467215642598, 0.5, 0.730980175370005))
expect_equal(stg5$GWAResult$ph1$LOD, 
             c(0.406196615759945, 0.301029995663981, 0.136094401214728))


