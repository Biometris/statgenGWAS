load(file = "testdata.rda")

### Test summary.GWAS

## Run some GWAS analsyis for testing.

stg0 <- runSingleTraitGwas(gData = gDataTest, trials = 1)
stg1 <- runSingleTraitGwas(gData = gDataTest, traits = 2, trials = "ph1", 
                           thrType = "small", nSnpLOD = 1)
stgM0 <- runSingleTraitGwas(gData = gDataTest, trials = 1, GLSMethod = "multi")

expect_error(summary(stg0, trials = TRUE),
             "trials should be a character or numeric vector")
expect_error(summary(stg0, trials = 7),
             "All trials should be in object")

sumStg0 <- capture.output(summary(stg0))
expect_true(any(grepl(pattern = "Traits analysed: X1, X2, X3, X4, X5", 
                      x = sumStg0)))
expect_true(any(grepl(pattern = "Data are available for 4 SNPs", x = sumStg0)))
expect_true(any(grepl(pattern = "0 of them were not analyzed", x = sumStg0)))
expect_true(any(grepl(pattern = "Genetic variance: 9.80734e-05", x = sumStg0)))
expect_true(any(grepl(pattern = "Residual variance: 2.16021", x = sumStg0)))
expect_true(any(grepl(pattern = "LOD-threshold: 1.90309", x = sumStg0)))
expect_true(any(grepl(pattern = "No significant SNPs found", x = sumStg0)))
expect_true(any(grepl(pattern = "No genomic control correction", x = sumStg0)))
expect_true(any(grepl(pattern = "inflation-factor: 1.179", x = sumStg0)))

sumStg1 <- capture.output(summary(stg1))
expect_true(any(grepl(pattern = "Number of significant SNPs: 1", x = sumStg1)))
expect_true(any(grepl(pattern = "Smallest p-value among the significant SNPs: 0.2249729",
                      x = sumStg1)))
expect_true(any(grepl(pattern = "Largest p-value among the significant SNPs: 0.2249729",
                      x = sumStg1)))
expect_true(any(grepl(pattern = "LOD-score: 0.6478698", x = sumStg1)))

sumStgM0 <- capture.output(summary(stgM0))
expect_false(any(grepl(pattern = "Genetic variance", x = sumStgM0)))
expect_false(any(grepl(pattern = "Residual variance", x = sumStgM0)))

