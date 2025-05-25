#' Launch GUIdedRNA Shiny Application
#'
#' @description Launches the GUIdedRNA Shiny application for guided RNA-seq analysis
#' @param port Port number for the application (default: 3838)
#' @param host Host address (default: "0.0.0.0")
#' @param launch.browser Whether to launch browser automatically (default: TRUE)
#' @return Starts the Shiny application
#' @export
#' @examples
#' \dontrun{
#' launch_GUIdedRNA()
#' }
launch_GUIdedRNA <- function(port = 3838, host = "0.0.0.0", launch.browser = TRUE) {
  app_dir <- system.file("shiny-app", package = "GUIdedRNA")
  
  if (app_dir == "") {
    stop("Could not find Shiny app directory. Try re-installing GUIdedRNA.", call. = FALSE)
  }
  
  message("Starting GUIdedRNA application...")
  shiny::runApp(
    appDir = app_dir,
    port = port,
    host = host,
    launch.browser = launch.browser
  )
}