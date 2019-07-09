context("plots")

qtlDat <- data.frame(trait = c("X1", "X2"), effect = c(0.3, -0.2),
                     chr = c(1, 2), pos = c(2, 3))

test_that("qq plot functions properly", {
  pVals <- runif(n = 50)
  p <- qqPlot(pValues = pVals, output = FALSE)
  expect_is(p, "ggplot")
  p1 <- qqPlot(pValues = pVals, title = "Test title", output = FALSE)
  expect_equal(p1$labels$title, "Test title")
})

test_that("GWAS qq plot functions properly", {
  stg <- runSingleTraitGwas(gDataTest)
  expect_error(plot(stg, plotType = "qq"), "multiple environments detected")
  expect_error(plot(stg, plotType = "qq", environment = "ph1"),
               "multiple traits detected")
  p <- plot(stg, type = "qq", environment = "ph1", trait = "X1", output = FALSE)
  expect_is(p, "ggplot")
})

test_that("qtl plot functions properly", {
  p <- qtlPlot(data = qtlDat, map = map, output = FALSE)
  expect_is(p, "ggplot")
  p1 <- qtlPlot(data = qtlDat, map = map, output = FALSE, yLab = "Testlab")
  expect_equal(p1$labels$y, "Testlab")
})

test_that("GWAS qtl plot functions properly", {
  stg <- runSingleTraitGwas(gDataTest, thrType = "fixed", LODThr = 0.5)
  expect_error(plot(stg, plotType = "qtl"), "multiple environments detected")
  p <- plot(stg, plotType = "qtl", environment = "ph1", output = FALSE)
  expect_is(p, "ggplot")
})

test_that("manhattan plot functions properly", {
  map <- data.frame(chr = rep(1:2, each = 3), cumPos = 1:6)
  p <- manhattanPlot(xValues = 1:6, yValues = 3:8, map = map, output = FALSE)
  expect_is(p, "ggplot")
  p1 <- manhattanPlot(xValues = 1:6, yValues = 3:8, map = map, xLab = "labx",
                      yLab = "laby", output = FALSE)
  expect_equal(p1$labels$x, "labx")
  expect_equal(p1$labels$y, "laby")
})

test_that("GWAS manhattan plot functions properly", {
  stg <- runSingleTraitGwas(gDataTest)
  expect_error(plot(stg, plotType = "manhattan"), "multiple environments detected")
  expect_error(plot(stg, plotType = "manhattan", environment = "ph1"),
               "multiple traits detected")
  p <- plot(stg, type = "manhattan", environment = "ph1", trait = "X1",
            output = FALSE)
  expect_is(p, "ggplot")
})
