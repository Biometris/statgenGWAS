load(file = "testdata.rda")

### Test runSingleTraitGwas with single kinship.

stg0 <- runSingleTraitGwas(gData = gDataTest, trials = 1)
stg1 <- runSingleTraitGwas(gData = gDataTest)
result1 <- runSingleTraitGwas(gData = gDataTest, trials = 1,
                              covar = "V1")[["GWAResult"]]
result1a <- runSingleTraitGwas(gData = gDataTest, traits = 2:6, trials = 1, 
                               covar = "V1")[["GWAResult"]]
result1b <- runSingleTraitGwas(gData = gDataTest, trials = 1,
                               covar = 1)[["GWAResult"]]
result2 <- runSingleTraitGwas(gData = gDataTest, trials = 1,
                              snpCov = "M2")[["GWAResult"]]
result3 <- runSingleTraitGwas(gData = gDataTest, trials = 1, covar = "V1",
                              snpCov = "M2")[["GWAResult"]]

## Check output structure.

expect_true(inherits(stg0, "GWAS"))
expect_equal(length(stg0), 5)
expect_equal(names(stg0), 
             c("GWAResult", "signSnp", "kinship", "thr", "GWASInfo"))
expect_equal(names(stg0[["GWASInfo"]]), 
             c("call", "remlAlgo", "thrType", "MAF", "GLSMethod",
               "kinshipMethod", "varComp", "genomicControl", "inflationFactor"))
expect_equal(length(stg0[["GWASInfo"]][["varComp"]][["ph1"]]), 5)
expect_true(inherits(stg0[["GWAResult"]], "list"))
expect_equal(length(stg0[["GWAResult"]]), 1)
expect_equal(names(stg0[["GWAResult"]]), "ph1")
expect_equal(length(stg1[["GWAResult"]]), 2)
expect_equal(names(stg1[["GWAResult"]]), c("ph1", "ph2"))

## Check results.

expect_equal_to_reference(stg0[["GWAResult"]][["ph1"]],
                          "stg0_GWARes_ph1")
expect_equal_to_reference(result1[["GWAResult"]][["ph1"]],
                          "result1_GWARes_ph1")
expect_equal_to_reference(result2[["GWAResult"]][["ph1"]],
                          "result2_GWARes_ph1")
expect_equal_to_reference(result3[["GWAResult"]][["ph1"]],
                          "result3_GWARes_ph1")

## Check results for traits containing NAs.

expect_equal_to_reference(stg1[["GWAResult"]][["ph2"]],
                          "stg1_GWARes_ph2")

## Check that specifying traits by number functions correctly.
expect_equal(result1, result1a)

## Check that specifying covar by number functions correctly.
expect_equal(result1, result1b)

### Test runSingleTraitGwas with chromosome specific kinship.

stgM0 <- runSingleTraitGwas(gData = gDataTest, trials = 1, GLSMethod = "multi")
stgM1 <- runSingleTraitGwas(gData = gDataTest, GLSMethod = "multi")
resultM1 <- runSingleTraitGwas(gData = gDataTest, trials = 1, 
                               GLSMethod = "multi", covar = "V1")[["GWAResult"]]
resultM2 <- runSingleTraitGwas(gData = gDataTest, trials = 1, 
                               GLSMethod = "multi", snpCov = "M2")[["GWAResult"]]
resultM3 <- runSingleTraitGwas(gData = gDataTest, trials = 1, covar = "V1",
                               GLSMethod = "multi", snpCov = "M2")[["GWAResult"]]

## Check output structure.

expect_true(inherits(stgM0, "GWAS"))
expect_equal(length(stgM0), 5)
expect_equal(names(stgM0), 
             c("GWAResult", "signSnp", "kinship", "thr", "GWASInfo"))
expect_equal(names(stgM0[["GWASInfo"]]), 
             c("call", "remlAlgo", "thrType", "MAF", "GLSMethod", 
               "kinshipMethod", "varComp", "genomicControl", "inflationFactor"))
expect_true(inherits(stgM0[["kinship"]], "list"))
expect_equal(length(stgM0[["kinship"]]), 2)
expect_true(inherits(stgM0[["kinship"]][[1]], "matrix"))
expect_equal(length(stgM0[["GWASInfo"]][["varComp"]][["ph1"]]), 5)
expect_equal(length(stgM0[["GWASInfo"]][["varComp"]][["ph1"]][[1]]), 2)
expect_true(inherits(stgM0[["GWAResult"]], "list"))
expect_equal(length(stgM0[["GWAResult"]]), 1)
expect_equal(names(stgM0[["GWAResult"]]), "ph1")
expect_equal(length(stgM1[["GWAResult"]]), 2)
expect_equal(names(stgM1[["GWAResult"]]), c("ph1", "ph2"))
expect_equal(stgM0[["GWASInfo"]][["GLSMethod"]] , "multi")

## Check results.

expect_equal_to_reference(stgM0[["GWAResult"]][["ph1"]],
                          "stgM0_GWARes_ph1")
expect_equal_to_reference(resultM1[["GWAResult"]][["ph1"]],
                          "resultM1_GWARes_ph1")
expect_equal_to_reference(resultM2[["GWAResult"]][["ph1"]],
                          "resultM2_GWARes_ph1")
expect_equal_to_reference(resultM3[["GWAResult"]][["ph1"]],
                          "resultM3_GWARes_ph1")

## Check results for traits containing NAs.

expect_equal_to_reference(stgM1[["GWAResult"]][["ph2"]],
                          "stgM1_GWARes_ph2")

