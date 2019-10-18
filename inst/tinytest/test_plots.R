load(file = "testdata.rda")

### Test plotting funtions.

## General input tests.

# Run a simple GWAS used for most plotting tests.
stg <- runSingleTraitGwas(gDataTest)

# These tests are identical for all plotTypes. 
# Only need to be checked once.

expect_error(plot(stg, trial = TRUE),
             "trial should be a character or numerical value")
expect_error(plot(stg, trial = 3), "trial should be in x")

# Plotting is always done for a single trial.
expect_error(plot(stg, "multiple trials detected"))

## Test qq plot

# Check on random p-Values.

expect_error(statgenGWAS:::qqPlot(pValues = 0:2, output = FALSE),
             "pValues should be an numeric vector with values between 0 and 1")

pVals <- runif(n = 50)
p <- statgenGWAS:::qqPlot(pValues = pVals, output = FALSE)
expect_true(inherits(p, "ggplot"))
p1 <- statgenGWAS:::qqPlot(pValues = pVals, title = "Test title", 
                           output = FALSE)
expect_equal(p1$labels$title, "Test title")

# Check for result of GWAS.

expect_error(plot(stg, plotType = "qq"), "multiple trials detected")
expect_error(plot(stg, plotType = "qq", trial = "ph1"),
             "multiple traits detected")
p <- plot(stg, plotType = "qq", trial = "ph1", trait = "X1", output = FALSE)
expect_true(inherits(p, "ggplot"))

## Test qtl plot

# Check on random qtl data.

qtlDat <- data.frame(trait = c("X1", "X2"), effect = c(0.3, -0.2),
                     chr = c(1, 2), pos = c(2, 3), sign = c(TRUE, FALSE))

expect_error(statgenGWAS:::qtlPlot(dat = 1),
             "dat should be a data.frame")
expect_error(statgenGWAS:::qtlPlot(dat = qtlDat, map = 1),
             "map should be a data.frame")

p <- statgenGWAS:::qtlPlot(dat = qtlDat, map = map, output = FALSE)
expect_true(inherits(p, "ggplot"))
p1 <- statgenGWAS:::qtlPlot(dat = qtlDat, map = map, output = FALSE, 
                            yLab = "Testlab")
expect_equal(p1$labels$y, "Testlab")

# Check for result of GWAS.

expect_error(plot(stg, plotType = "qtl", trial = "ph1", output = FALSE),
             "No significant SNPs found. No plot can be made")

stg1 <- runSingleTraitGwas(gDataTest, thrType = "fixed", LODThr = 0.2)
expect_error(plot(stg1, plotType = "qtl"), "multiple trials detected")
p <- plot(stg1, plotType = "qtl", trial = "ph1", output = FALSE)
expect_true(inherits(p, "ggplot"))

# Check yThr
expect_error(plot(stg, plotType = "qtl", trial = "ph1", yThr = -1, 
                  output = FALSE),
             "yThr should be a single numerical value greater than 0")
expect_silent(plot(stg, plotType = "qtl", trial = "ph1", yThr = 0.2, 
                   output = FALSE))

# Check option chr.
expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", chr = 3,
                  output = FALSE),
             "Select at least one valid chromosome for plotting")
p1 <- plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", chr = 1,
           output = FALSE)
expect_equal(nrow(p1$data), 7)

# Check option normalize.

expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  normalize = 1, output = FALSE),
             "normalize should be a single logical")
p1 <- plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
           normalize = TRUE, output = FALSE)

# Check option sortData.
expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  sortData = TRUE, output = FALSE),
             "sortData should be either FALSE or a single character")
expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  sortData = "sortCol", output = FALSE),
             "dat lacks the following columns: sortCol")
expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  sortData = "trait", output = FALSE),
             "sortData should be a numerical column")

# Data contains no numerical column for sorting. 
# Add one manually.
stg1a <- stg1
stg1a$GWAResult$ph1$sortCol <- rep(c(3, 4, 5, 1, 2), each = 3)
expect_silent(p1 <- plot(stg1a, plotType = "qtl", trial = "ph1", trait = "X1", 
                         sortData = "sortCol", output = FALSE))

# Check option binPositions
binPos <- data.frame(chr = 1)
binPos1 <- data.frame(chr = 1, pos = 2)

expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  binPositions = "binPos", output = FALSE),
             "binPositions should be either NULL or an data.frame")
expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  binPositions = binPos, output = FALSE),
             "binPositions lacks the following columns: pos")
expect_silent(p1 <- plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                         binPositions = binPos1, output = FALSE))

# Check pdf output.
# Create tmpfile.
tmpPptx <- tempfile(fileext = ".pptx")
tmpPpt <- tempfile(fileext = ".ppt")

expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  exportPptx = 1, pptxName = tmpPptx, output = FALSE),
             "exportPptx should be a single logical")
expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  exportPptx = TRUE, pptxName = NULL, output = FALSE),
             "pptxName should be a single character string")
expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  exportPptx = TRUE, pptxName = tmpPpt, output = FALSE),
             "should have '.pptx' extension.")

expect_silent(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                   exportPptx = TRUE, pptxName = tmpPptx, output = FALSE))

## Test manhattan plot

# Check on random data.

map <- data.frame(chr = rep(1:2, each = 3), cumPos = 1:6)
p <- statgenGWAS:::manhattanPlot(xValues = 1:6, yValues = 3:8, map = map, 
                                 output = FALSE)
expect_true(inherits(p, "ggplot"))
p1 <- statgenGWAS:::manhattanPlot(xValues = 1:6, yValues = 3:8, map = map, 
                                  xLab = "labx", yLab = "laby", output = FALSE)
expect_equal(p1$labels$x, "labx")
expect_equal(p1$labels$y, "laby")

# Check for result of GWAS.

# Manhattan plots are always made for single traits.
expect_error(plot(stg, plotType = "manhattan", trial = "ph1"), 
             "multiple traits detected")

p <- plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1", 
          output = FALSE)
expect_true(inherits(p, "ggplot"))

# Check option chr.
expect_error(plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
                  chr = 3, output = FALSE),
             "Select at least one valid chromosome for plotting")
p1 <- plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
           chr = 1, output = FALSE)
expect_equal(nrow(p1$data), 2)

# Check option effects.
expect_error(plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
                  effects = "M5", output = FALSE),
             "All known effects should be in the map")
p1 <- plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
           effects = "M1", output = FALSE)
# For the colored effect point a new layer should be added in p1.
expect_equal(length(p$layers), length(p1$layers) - 1)

# Check option lod.
expect_error(plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
                  lod = -1, output = FALSE),
             "lod should be a single numerical value greater than 0")
p1 <- plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
           lod = .5, output = FALSE)
# Just one SNP left after sampling.
expect_equal(nrow(p1$data), 1)

# Check yThr
expect_error(plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
                  yThr = -1, output = FALSE),
             "yThr should be a single numerical value greater than 0")

p1 <- plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
           yThr = .1, output = FALSE)
# For the colored points above the new thr a new layer should be added in p1.
expect_equal(length(p$layers), length(p1$layers) - 1)

# Check combination of yThr and effect to create true neg/false pos.
p1 <- plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
           effects = c("M2", "M3"), yThr = .1, output = FALSE)
# This should add one of each: true positive, false positive and false negative.
# So three extra layers with colors should be added.
expect_equal(length(p$layers), length(p1$layers) - 3)

# Check that significant SNPs are picked up directly from GWAS output.
p1 <- plot(stg1, plotType = "manhattan", trial = "ph1", trait = "X1", 
           output = FALSE)
# One new layer should be added for the significant SNP.
expect_equal(length(p$layers), length(p1$layers) - 1)
