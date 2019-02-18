.onLoad <- function(libname = find.package("gwas"),
                    pkgname = "gwas"){
  ## CRAN Note avoidance.
  if (getRversion() >= "2.15.1") {
    utils::globalVariables("j")
  }
  invisible()
}

.onUnload <- function(libpath) {
  library.dynam.unload("gwas", libpath)
}
