#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function (input, output, session) {
  # List the first level callModules here
  r <- shiny::reactiveValues(
    PATH_TO_DEMO = "D:/These/QuantifFiles",
    raw_counts = NULL,
    normalized_counts = NULL,
    normalized_counts_pre_filter = NULL,
    norm_factor = NULL,
    conditions = NULL,
    design = NULL
  )
  
  
  shiny::callModule(mod_context_server, "context_ui_1")
  shiny::callModule(mod_import_data_server, "import_data_ui_1", r)
  shiny::callModule(mod_normalisation_server, "normalisation_ui_1", r)
}
