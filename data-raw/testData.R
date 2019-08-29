## Create data for testing functions.
## Restricted and anonymized version of DROPs data.
load(system.file("extdata", "testData.RData", package = "statgenGWAS"))

## Create a dataset for unit testing.
set.seed(1234)
y <- 1:10
X <- matrix(sample(x = c(0, 1), size = 30, replace = TRUE), nrow = 10)
Sigma <- matrix(runif(n = 100), nrow = 10)
## Assure sigma is symmetric.
Sigma <- tcrossprod(Sigma)
covs <- matrix(runif(n = 20, max = 100), nrow = 10)
pheno <- data.frame(genotype = paste0("G", 1:10),
                    matrix(rnorm(50, mean = 10, sd = 2), nrow = 10))
## Add random missing values to pheno2.
pheno2 <- pheno
for (i in 2:6) {
  pheno2[sample(x = 1:10, size = 2), i] <- NA
}
map <- data.frame(chr = c(1, 1, 2), pos = 1:3)
rownames(X) <- rownames(Sigma) <- colnames(Sigma) <- rownames(covs) <-
  paste0("G", 1:10)
colnames(X) <- rownames(map) <- paste0("M", 1:3)
## Create gData object.
gDataTest <- createGData(map = map, geno = X, kin = Sigma,
                         pheno = list(ph1 = pheno, ph2 = pheno2),
                         covar = as.data.frame(covs))
## Export to package
usethis::use_data(gDataTest, X, y, Sigma, map, covs,
                  internal = TRUE, overwrite = TRUE)
