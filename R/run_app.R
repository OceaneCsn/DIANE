#' Run the Shiny Application
#'
#' @param ... A series of options to be used inside the app.
#' @param server_version TRUE if the app is deployed on web server, 
#' else (default), FALSE
#' @return shiny application
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(server_version = FALSE,
  ...
) {
  golem::with_golem_options(
    app = shinyApp(
      ui = app_ui, 
      server = app_server,
      options = list(host = "0.0.0.0")
    ),
    golem_opts = list("server_version" = server_version , ...)
  )
}
