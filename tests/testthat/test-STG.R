context("single trait GWAS")

set.seed(1234)
y <- 1:10
X <- matrix(sample(x = c(0, 1), size = 30, replace = TRUE), nrow = 10)
Sigma <- matrix(runif(n = 100), nrow = 10)
Sigma <- tcrossprod(Sigma)
covs <- matrix(runif(n = 20, max = 100), nrow = 10)
pheno <- data.frame(genotype = paste0("G", 1:10),
                    matrix(rnorm(50, mean = 10, sd = 2), nrow = 10))
map <- data.frame(chr = c(1, 1, 2), pos = 1:3)
rownames(X) <- rownames(Sigma) <- colnames(Sigma) <- rownames(covs) <-
  paste0("G", 1:10)
colnames(X) <- rownames(map) <- paste0("M", 1:3)
gDataTest <- createGData(map = map, geno = X, kin = Sigma,
                         pheno = list(ph1 = pheno, ph2 = pheno),
                         covar = as.data.frame(covs))

stg0 <- runSingleTraitGwas(gData = gDataTest, environments = 1)
stg01 <- runSingleTraitGwas(gData = gDataTest)
result1 <- runSingleTraitGwas(gData = gDataTest, environments = 1,
                              covar = "V1")$GWAResult
result2 <- runSingleTraitGwas(gData = gDataTest, environments = 1,
                              snpCov = "M2")$GWAResult
result3 <- runSingleTraitGwas(gData = gDataTest, environments = 1, covar = "V1",
                              snpCov = "M2")$GWAResult

test_that("runSingleTraitGwas produces correct output structure", {
  expect_is(stg0, "GWAS")
  expect_length(stg0, 5)
  expect_named(stg0, c("GWAResult", "signSnp", "kinship", "thr", "GWASInfo"))
  expect_is(stg0$GWAResult, "list")
  expect_length(stg0$GWAResult, 1)
  expect_named(stg0$GWAResult, "ph1")
  expect_length(stg01$GWAResult, 2)
  expect_named(stg01$GWAResult, c("ph1", "ph2"))
})

test_that("runSingleTraitGWas produces correct p-values", {
  expect_equal(stg0$GWAResult[[1]]$pValue,
               c(0.517279205274389, 0.913713471945949, 0.628864085064797,
                 0.0807864803940613, 0.857734879152358, 0.0951298087141793,
                 0.616384737625595, 0.994787446465712, 0.421203350051828,
                 0.183886590676029, 0.973491209528074, 0.57064757354885,
                 0.588456561785548, 0.367143146285209, 0.905504194974234))
  expect_equal(result1[[1]]$pValue,
               c(0.269933551571042, 0.984965446392588, 0.648441699544188,
                 0.0626857428326674, 0.866711176390765, 0.120055932703493,
                 0.770861361333917, 0.940364514318174, 0.451352288953829,
                 0.22823190214397, 0.946236384390807, 0.592108077080844,
                 0.619302653579761, 0.405860998685857, 0.911116358438535))
  expect_equal(result2[[1]]$pValue,
               c(0.52489597499325, 0.894152021294656, 0.602218414048816,
                 0.0970424392791122, 0.857734879152358, 0.0928126890021481,
                 0.626210599310249, 0.99304196419329, 0.437394594066136,
                 0.212511818760432, 0.973491209528084, 0.579458353874296,
                 0.648912926295298, 0.367143146285209, 0.918479093445766))
  expect_equal(result3[[1]]$pValue,
               c(0.308728292832841, 0.984965446392575, 0.661002824234325,
                 0.0830816331311562, 0.866711176390759, 0.121327025768191,
                 0.802520086738191, 0.928635114511299, 0.464942717455927,
                 0.26495153781054, 0.946236384390817, 0.598229839636105,
                 0.657680037715817, 0.405860998685858, 0.926119151352394))
})
