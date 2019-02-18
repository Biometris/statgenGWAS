## Create data for testing functions.
## Restricted and anonymized version of DROPs data.
load(system.file("extdata", "testData.RData", package = "gwas"))

# Export to package
usethis::use_data(Y, K, X, internal = TRUE, overwrite = TRUE)
