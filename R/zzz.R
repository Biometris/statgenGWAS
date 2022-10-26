.onUnload <- function(libpath) {
  library.dynam.unload("statgenGWAS", libpath)
}
