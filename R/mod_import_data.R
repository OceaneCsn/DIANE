library(shinipsum)
#' import_data UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_import_data_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h2("A Random Text 2"),
          tableOutput(ns("text")))
}

#' import_data Server Function
#'
#' @noRd
mod_import_data_server <- function(input, output, session) {
  ns <- session$ns
  
  output$text <- renderText({
    random_text(nwords = 100)
    
  })
  
}

## To be copied in the UI
# mod_import_data_ui("import_data_ui_1")

## To be copied in the server
# callModule(mod_import_data_server, "import_data_ui_1")
