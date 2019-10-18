### Test createGData.

## Create some data for testing.

# Create genotypic data.
geno <- matrix(sample(x = c(0, 1, 2), size = 15, replace = TRUE), nrow = 3)
dimnames(geno) <- list(paste0("G", 1:3), paste0("M", 1:5))

# Construct map.
map <- data.frame(chr = c(1, 1, 2, 2, 2), pos = 1:5,
                  row.names = paste0("M", 1:5))

# Compute kinship matrices.
kin <- kinship(X = geno, method = "IBS")

# Create phenotypic data.
pheno <- data.frame(paste0("G", 1:3),
                    matrix(rnorm(n = 12, mean = 50, sd = 5), nrow = 3),
                    stringsAsFactors = FALSE)
dimnames(pheno) = list(paste0("G", 1:3), c("genotype", paste0("T", 1:4)))

# Construct covariate.
covar <- data.frame(C1 = c("a", "a", "b"), row.names = paste0("G", 1:3))

# Create gData object with all inputs used for checking that 
# overwriting works properly.
gd0 <- createGData(geno = geno, map = map, kin = kin, pheno = pheno, 
                   covar = covar)

## Check inputs.

expect_error(createGData(),
             "At least one of geno, map, kin, pheno and covar should")
expect_error(createGData(gData = 1),
             "gData object should be of class gData")

# Checks for map.
expect_error(createGData(map = 1), "map should be a data.frame")
expect_error(createGData(map = map[, 1, drop = FALSE]), 
             "chr and pos should be columns in map")
# Create map with character pos column.
map2 <- map
map2$pos <- as.character(map2$pos)
expect_error(createGData(map = map2), "pos should be a numeric column in map")
# If gData object already contains a map a warning should be returned.
expect_warning(createGData(gData = gd0, map = map), 
               "Existing map will be overwritten")

# Create map without rownames.
# createGData should generate rownames based on chr and pos.
# Should be made unique in case of duplicate chr x pos.
map3 <- rbind(map, map[1, ])
rownames(map3) <- NULL
expect_warning(gd <- createGData(map = map3),
               "Names constructed from chromosome and position")
expect_equal(rownames(gd$map), 
             c("chr1_1_1", "chr1_1_2", "chr1_2", "chr2_3", "chr2_4", "chr2_5"))

# Checks for pheno.
expect_error(createGData(pheno = 1),
             "pheno should be a data.frame or a list of data.frames")
gd <- createGData(pheno = pheno)
# A single input data.frame should be converted to a list
expect_true(inherits(gd$pheno, "list"))
expect_equal(names(gd$pheno), "pheno")
expect_true(inherits(gd$pheno$pheno, "data.frame"))

# An unnamed list of data.frames should get default names.
expect_warning(gd <- createGData(pheno = list(pheno, pheno)),
               "pheno contains no trial names. Default names added")
expect_equal(names(gd$pheno), c("Trial1", "Trial2"))

# A partially unnamed list should retain its names and have defaults added.
expect_warning(gd <- createGData(pheno = list(ph1 = pheno, pheno)),
               "Some data.frames in pheno contain no trial names.")
expect_equal(names(gd$pheno), c("ph1", "Trial2"))

# A fully named list should be returned directly.
gd <- createGData(pheno = list(ph1 = pheno, ph2 = pheno))
expect_equal(names(gd$pheno), c("ph1", "ph2"))

# First column of all data.frames should be "genotype".
# Create data.frame with missing genotype column.
pheno2 <- pheno
colnames(pheno2)[1] <- "geno"
expect_error(createGData(pheno = pheno2), 
             "First column in pheno should be genotype")
expect_error(createGData(pheno = list(ph1 = pheno, ph2 = pheno2)), 
             "First column in pheno should be genotype")

# All other columns in pheno should be numerical.
pheno3 <- pheno
pheno3$T3 <- as.character(pheno3$T3)
expect_error(createGData(pheno = pheno3), 
             "All trait columns in pheno should be numerical")
expect_error(createGData(pheno = list(ph1 = pheno, ph2 = pheno3)), 
             "All trait columns in pheno should be numerical")
                  
# If gData object already contains a pheno object a warning should be given.
expect_warning(createGData(gData = gd0, pheno = pheno), 
               "Existing pheno will be overwritten")

# Checks for geno.
expect_error(createGData(geno = 1),
             "geno should be a matrix, data.frame or an array")

# Check that input is converted to matrix.
gd <- createGData(geno = as.data.frame(geno))
expect_true(inherits(gd$markers, "matrix"))
expect_equal(gd$markers, geno)

# When geno has no rownames default names should be added.
# If pheno is used as argument names should be taken from pheno.

# Create geno data without rownames.
geno2 <- geno
rownames(geno2) <- NULL

expect_warning(gd <- createGData(geno = geno2, pheno = pheno),
               "geno contains no genotype names. Names taken from pheno")
expect_equal(rownames(gd$markers), c("G1", "G2", "G3"))

# If dimensions differ between geno and pheno function should stop.
expect_error(createGData(geno = geno2, pheno = pheno[1:2, ]),
             "Dimensions between geno and pheno differ")

# If pheno is not used default names are constructed.
expect_warning(gd <- createGData(geno = geno2),
               "geno contains no genotype names. Default names used")
expect_equal(rownames(gd$markers), c("g1", "g2", "g3"))

# When geno has no rownames marker names are taken from map.
# No default names are constructed.

# Create geno data without colnames.
geno3 <- geno
colnames(geno3) <- NULL

expect_error(createGData(geno = geno3),
             "geno contains no marker names. Map not available")
expect_error(createGData(geno = geno3, map = map[1:2, ]),
             "Dimensions between geno and map differ")
expect_warning(gd <- createGData(geno = geno3, map = map),
               "geno contains no marker names. Names taken from map")
expect_equal(colnames(gd$markers), c("M1", "M2", "M3", "M4", "M5"))

# Markers that are in geno but not in map should be removed.
# Map information is always needed for functions involving markers.
expect_warning(gd <- createGData(geno = geno, map = map[1:2, ]),
               "Extra markers will be removed")
expect_equal(colnames(gd$markers), c("M1", "M2"))

# If gData object already contains a marker object a warning should be given.
expect_warning(createGData(gData = gd0, geno = geno), 
               "Existing geno will be overwritten")

# Checks for kin.
expect_error(createGData(kin = 1), 
             "kin should be a matrix or a list of matrices")

# If kin is a list of matrices its length should match the number of chr in map.
expect_error(createGData(kin = list(kin), map = map), 
             "kin should be the same length as the number of chromosomes")

# kin list is named, names should match names of chromosomes.
expect_error(createGData(kin = list(a = kin, b = kin), map = map), 
             "Names of kin should correspond to names of chromosomes in map")
expect_silent(createGData(kin = list("1" = kin, "2" = kin), map = map))

# An unnamed kin list should have default names added.
expect_warning(gd <- createGData(kin = list(kin, kin), map = map),
               "kin contains no names. Default names added.")
expect_equal(names(gd$kinship), c("1", "2"))

# kin should have row and column names

# Create kin without rownames.
kin2 <- kin
rownames(kin2) <- NULL
expect_error(createGData(kin = kin2), 
             "Row and column names in kin cannot be NULL")
expect_error(createGData(kin = list("1" = kin, "2" = kin2)), 
             "Row and column names in kin cannot be NULL")

# Row and column names of kin should be in genotypes of markers.
# No need for an exact match.
expect_error(createGData(kin = kin, geno = geno[1:2, ]), 
             "Row and column names of kin should be in row names of geno")
expect_error(createGData(kin = list("1" = kin, "2" = kin), geno = geno[1:2, ]), 
             "Row and column names of kin should be in row names of geno")

# kin should be converted to a matrix.
gd <- createGData(kin = Matrix::Matrix(kin))
expect_true(inherits(gd$kinship, "matrix"))
expect_equal(gd$kinship, kin)

# If gData object already contains a kinship matrix a warning should be given.
expect_warning(createGData(gData = gd0, kin = kin), 
               "Existing kinship will be overwritten")

# Checks for covar.

expect_error(createGData(covar = 1), "covar should be a data.frame")

# covar can only contain factors, characters and numeric columns.

# Create covar with boolean column.
covar2 <- covar
covar2$C2 <- c(TRUE, FALSE, TRUE)
expect_error(createGData(covar = covar2),
             "columns in covar should be numeric, character or factor columns")

# If gData object already contains a kinship matrix a warning should be given.
expect_warning(createGData(gData = gd0, covar = covar), 
               "Existing covar will be overwritten")

## Test summary.gData
expect_silent(sumGd <- capture.output(summary(gd0)))
expect_true(any(grepl(pattern = "Number of markers: 5", x = sumGd)))
expect_true(any(grepl(pattern = "Number of chromosomes: 2", x = sumGd)))
expect_true(any(grepl(pattern = "Number of genotypes: 3", x = sumGd)))
expect_true(any(grepl(pattern = "Number of trials: 1", x = sumGd)))
expect_true(any(grepl(pattern = "Number of traits: 4", x = sumGd)))
expect_true(any(grepl(pattern = "Number of covariates: 1", x = sumGd)))
