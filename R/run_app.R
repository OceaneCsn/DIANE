#' Run the Shiny Application
#'
#' @param ... A series of options to be used inside the app.
#' @return shiny application
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(
  ...
) {
  with_golem_options(
    app = shinyApp(
      ui = app_ui, 
      server = app_server,
      options = list(host = "0.0.0.0")
    ), 
    
    golem_opts = list( ...)
  )
}
