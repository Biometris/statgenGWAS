load(file = "testdata.rda")

### Test plotting funtions.

## Test qq plot

# Check on random p-Values

pVals <- runif(n = 50)
p <- statgenGWAS:::qqPlot(pValues = pVals, output = FALSE)
expect_true(inherits(p, "ggplot"))
p1 <- statgenGWAS:::qqPlot(pValues = pVals, title = "Test title", output = FALSE)
expect_equal(p1$labels$title, "Test title")

# Check for result of GWAS.

stg <- runSingleTraitGwas(gDataTest)
expect_error(plot(stg, plotType = "qq"), "multiple trials detected")
expect_error(plot(stg, plotType = "qq", trial = "ph1"),
             "multiple traits detected")
p <- plot(stg, type = "qq", trial = "ph1", trait = "X1", output = FALSE)
expect_true(inherits(p, "ggplot"))

## Test qtl plot

# Check on random qtl data.

qtlDat <- data.frame(trait = c("X1", "X2"), effect = c(0.3, -0.2),
                     chr = c(1, 2), pos = c(2, 3))

p <- statgenGWAS:::qtlPlot(dat = qtlDat, map = map, output = FALSE)
expect_true(inherits(p, "ggplot"))
p1 <- statgenGWAS:::qtlPlot(dat = qtlDat, map = map, output = FALSE, 
                            yLab = "Testlab")
expect_equal(p1$labels$y, "Testlab")

# Check for result of GWAS.

stg <- runSingleTraitGwas(gDataTest, thrType = "fixed", LODThr = 0.5)
expect_error(plot(stg, plotType = "qtl"), "multiple trials detected")
p <- plot(stg, plotType = "qtl", trial = "ph1", output = FALSE)
expect_true(inherits(p, "ggplot"))

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

stg <- runSingleTraitGwas(gDataTest)
expect_error(plot(stg, plotType = "manhattan"), "multiple trials detected")
expect_error(plot(stg, plotType = "manhattan", trial = "ph1"),
             "multiple traits detected")
p <- plot(stg, type = "manhattan", trial = "ph1", trait = "X1", output = FALSE)
expect_true(inherits(p, "ggplot"))
