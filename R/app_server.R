#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function( input, output, session ) {
  # List the first level callModules here

  
  callModule(mod_context_server, "context_ui_1")
  callModule(mod_import_data_server, "import_data_ui_1")

}
