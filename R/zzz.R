#' @import shiny
#' @import shinydashboard
#' @export
.onLoad <- function(libname, pkgname) {
  packageStartupMessage("Loading GUIdedRNA package...")
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("GUIdedRNA: Guided RNA-seq Analysis Pipeline")
  packageStartupMessage("Use launch_GUIdedRNA() to start the application")
}