#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function( input, output, session ) {
  # List the first level callModules here
  r <- reactiveValues()
  r$PATH_TO_DEMO = "D:/These/QuantifFiles"
  
  
  callModule(mod_import_data_server, "import_data_ui_1", r)
  callModule(mod_context_server, "context_ui_1")
}

