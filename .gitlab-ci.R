installIfNeeded <- function(pkg,
                            quiet) {
  pkgPath <- find.package(pkg, quiet = quiet)
  if (length(pkgPath) == 0) {
    message("NOTE: pkg ", pkg, " missing, installing...")
    install.packages(pkg, quiet = FALSE)
  } else {
    update.packages(pkg, ask = FALSE)
  }
}

pkgsUpdate <- function(quiet = TRUE,
                       instPkgdown = FALSE) {
  installIfNeeded(pkg = "remotes", quiet = quiet)
  if (instPkgdown) {
    installIfNeeded(pkg = "pkgdown", quiet = quiet)
  }
  remotes::install_deps(dependencies = TRUE, quiet = quiet)
  cat("INSTALLED:\n")
  instld <- as.data.frame(installed.packages())
  rownames(instld) <- NULL
  print(instld[, c("Package", "Version")])
  return(invisible(TRUE))
}



