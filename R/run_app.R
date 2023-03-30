#' Run the Shiny Application
#'
#' @param ... A series of options to be used inside the app.
#'
#' @param server_version TRUE if the app is deployed on web server,
#' else (default), FALSE
#' @param seed seed for random state to ensure reproducible runs
#' along DIANE's pipeline
#' @param host IP address of the app 
#' (e.g "127.0.0.1" for localhost, "0.0.0.0", ...)
#'
#' @return shiny application
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <-
  function(server_version = FALSE,
           seed = round(runif(n = 1, min = 0, max = 2 ^ 7)),
           host = "127.0.0.1",
           ...) {
    golem::with_golem_options(
      app = shinyApp(ui = app_ui,
                     server = app_server, 
                     options = list(host = host)),
      golem_opts = list("server_version" = server_version,
                        "seed" = seed, ...)
    )
  }