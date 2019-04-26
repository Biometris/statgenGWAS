install_if_needed <- function(pkg, repos, quiet) {
  package_path <- find.package(pkg, quiet = quiet)
  if (length(package_path) == 0) {
    message("NOTE: pkg ", pkg, " missing, installing...")
    install.packages(pkg, repos = repos, quiet = quiet)
  }
}

gl_update_pkg_all <- function(repos = "https://cran.rstudio.com",
                              quiet = TRUE,
                              install_pkgdown = FALSE) {
  # update existing
  update.packages(ask = FALSE, repos = repos, quiet = quiet)

  install_if_needed(pkg = "devtools", repos = repos, quiet = quiet)
  if (install_pkgdown == TRUE) {
    install_if_needed(pkg = "pkgdown", repos = repos, quiet = quiet)
  }

  devtools::install_dev_deps(repos = repos, quiet = quiet, upgrade = TRUE)

  cat("INSTALLED:\n")
  instld <- as.data.frame(installed.packages())
  rownames(instld) <- NULL
  print(instld[, c("Package", "Version")])

  return(invisible(TRUE))
}
