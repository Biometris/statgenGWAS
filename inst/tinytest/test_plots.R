load(file = "testdata.rda")

### Test plotting functions.

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

## Test QQ-plot

# Check on random p-Values.

expect_error(statgenGWAS:::qqPlot(pValues = 0:2),
             "pValues should be an numeric vector with values between 0 and 1")

pVals <- runif(n = 50)
p <- statgenGWAS:::qqPlot(pValues = pVals)
expect_true(inherits(p, "ggplot"))
p1 <- statgenGWAS:::qqPlot(pValues = pVals, title = "Test title", 
                           output = FALSE)
expect_equal(p1$labels$title, "Test title")

# Check for result of GWAS.

expect_error(plot(stg, plotType = "qq"), "multiple trials detected")
expect_error(plot(stg, plotType = "qq", trial = "ph1"),
             "multiple traits detected")
p <- plot(stg, plotType = "qq", trial = "ph1", trait = "X1")
expect_true(inherits(p, "ggplot"))

## Test qtl plot

# Check on random qtl data.

qtlDat <- data.frame(trait = c("X1", "X2"), effect = c(0.3, -0.2),
                     chr = c(1, 2), pos = c(2, 3), sign = c(TRUE, FALSE))

expect_error(statgenGWAS:::qtlPlot(dat = 1),
             "dat should be a data.frame")
expect_error(statgenGWAS:::qtlPlot(dat = qtlDat, map = 1),
             "map should be a data.frame")

p <- statgenGWAS:::qtlPlot(dat = qtlDat, map = map)
expect_true(inherits(p, "ggplot"))
p1 <- statgenGWAS:::qtlPlot(dat = qtlDat, map = map, 
                            yLab = "Testlab")
expect_equal(p1$labels$y, "Testlab")

p1 <- statgenGWAS:::qtlPlot(dat = qtlDat, map = map, title = "Test title")
expect_equal(p1$labels$title, "Test title")

# Check for result of GWAS.

expect_error(plot(stg, plotType = "qtl", trial = "ph1"),
             "No significant SNPs found. No plot can be made")

stg1 <- runSingleTraitGwas(gDataTest, thrType = "fixed", LODThr = 0.2)
expect_error(plot(stg1, plotType = "qtl"), "multiple trials detected")
p <- plot(stg1, plotType = "qtl", trial = "ph1")
expect_true(inherits(p, "ggplot"))

# Check yThr
expect_error(plot(stg, plotType = "qtl", trial = "ph1", yThr = -1),
             "yThr should be a single numerical value greater than 0")
expect_silent(plot(stg, plotType = "qtl", trial = "ph1", yThr = 0.2))

# Check option chr.
expect_error(plot(stg1, plotType = "qtl", trial = "ph1", chr = 3),
             "Select at least one valid chromosome for plotting")
p1 <- plot(stg1, plotType = "qtl", trial = "ph1", chr = 1)
expect_equal(nrow(p1$data), 9)

# Check option normalize.

expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  normalize = 1),
             "normalize should be a single logical")
p1 <- plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
           normalize = TRUE)

# Check option sortData.

expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  sortData = TRUE),
             "sortData should be either FALSE or a single character")
expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  sortData = "sortCol"),
             "dat lacks the following columns: sortCol")
expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  sortData = "trait"),
             "sortData should be a numerical column")

# Data contains no numerical column for sorting. 
# Add one manually.
stg1a <- stg1
stg1a$GWAResult$ph1$sortCol <- rep(c(3, 4, 5, 1, 2), each = 4)
expect_silent(p1 <- plot(stg1a, plotType = "qtl", trial = "ph1", trait = "X1", 
                         sortData = "sortCol"))

# Check option binPositions
binPos <- data.frame(chr = 1)
binPos1 <- data.frame(chr = 1, pos = 2)

expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  binPositions = "binPos"),
             "binPositions should be either NULL or an data.frame")
expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  binPositions = binPos),
             "binPositions lacks the following columns: pos")
expect_silent(p1 <- plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                         binPositions = binPos1))

# Check pdf output.
# Create tmpfile.
tmpPptx <- tempfile(fileext = ".pptx")
tmpPpt <- tempfile(fileext = ".ppt")

expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  exportPptx = 1, pptxName = tmpPptx),
             "exportPptx should be a single logical")
expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  exportPptx = TRUE, pptxName = NULL),
             "pptxName should be a single character string")
expect_error(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                  exportPptx = TRUE, pptxName = tmpPpt),
             "should have '.pptx' extension")

expect_silent(plot(stg1, plotType = "qtl", trial = "ph1", trait = "X1", 
                   exportPptx = TRUE, pptxName = tmpPptx))

## Test manhattan plot

expect_error(statgenGWAS:::manhattanPlot(xValues = "1"),
             "xValues should be a numerical vector")
expect_error(statgenGWAS:::manhattanPlot(xValues = 1, yValues = "1"),
             "yValues should be a numerical vector")
expect_error(statgenGWAS:::manhattanPlot(xValues = 1, yValues = 1:2),
             "xValues and yValues should be of the same length")
expect_error(statgenGWAS:::manhattanPlot(xValues = 1:2, yValues = 1:2, 
                                         xSig = 1.3),
             "xSig should be an integer vector")
expect_error(statgenGWAS:::manhattanPlot(xValues = 1:2, yValues = 1:2, 
                                         xEffects = 1.3),
             "xEffects should be an integer vector")

# Check on random data.

map <- data.frame(chr = rep(1:2, each = 3), cumPos = 1:6)
p <- statgenGWAS:::manhattanPlot(xValues = 1:6, yValues = 3:8, map = map)
expect_true(inherits(p, "ggplot"))
p1 <- statgenGWAS:::manhattanPlot(xValues = 1:6, yValues = 3:8, map = map, 
                                  xLab = "labx", yLab = "laby", 
                                  title = "Test title")
expect_equal(p1$labels$x, "labx")
expect_equal(p1$labels$y, "laby")
expect_equal(p1$labels$title, "Test title")

# Check for result of GWAS.

# Manhattan plots are always made for single traits.
expect_error(plot(stg, plotType = "manhattan", trial = "ph1"), 
             "multiple traits detected")

p <- plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1")
expect_true(inherits(p, "ggplot"))

# Check option chr.
expect_error(plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
                  chr = 3),
             "Select at least one valid chromosome for plotting")
p1 <- plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
           chr = 1)
expect_equal(nrow(p1$data), 2)

# Check options startPos/endPos.
expect_error(plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
                  chr = 1, startPos = -1),
             "startPos should be a single numerical value between 0 and 2")
expect_error(plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
                  chr = 1, endPos = -1),
             "endPos should be a single numerical value greater than 0")
expect_error(plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
                  chr = 1, startPos = 1, endPos = 0),
             "Start position should be smaller than end position")
expect_error(plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
                  chr = 1, startPos = 1.5, endPos = 1.6),
             "No SNPs in selected range")
p1 <- plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
           chr = 1, startPos = 2, endPos = 2)
expect_equal(nrow(p1$data), 1)

# Check option effects.
expect_error(plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
                  effects = "M5"),
             "All known effects should be in the map")
p1 <- plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
           effects = "M1")
# For the colored effect point a new layer should be added in p1.
expect_equal(length(p$layers), length(p1$layers) - 1)

# Check option lod.
expect_error(plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
                  lod = -1),
             "lod should be a single numerical value greater than 0")
p1 <- plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
           lod = .5)
# Just one SNP left after sampling.
expect_equal(nrow(p1$data), 2)

# Check yThr
expect_error(plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
                  yThr = -1),
             "yThr should be a single numerical value greater than 0")

p1 <- plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
           yThr = .1)
# For the colored points above the new thr a new layer should be added in p1.
expect_equal(length(p$layers), length(p1$layers) - 1)

# Check combination of yThr and effect to create true neg/false pos.
p1 <- plot(stg, plotType = "manhattan", trial = "ph1", trait = "X1",
           effects = c("M2", "M3"), yThr = .1)
# This should add one of each: true positive, false positive and false negative.
# So three extra layers with colors should be added.
expect_equal(length(p$layers), length(p1$layers) - 3)

# Check that significant SNPs are picked up directly from GWAS output.
p1 <- plot(stg1, plotType = "manhattan", trial = "ph1", trait = "X1")
# One new layer should be added for the significant SNP.
expect_equal(length(p$layers), length(p1$layers) - 1)

## Cleanup temporary files.
unlink(tmpPptx)
unlink(tmpPpt)
