library(shinipsum)
#' context UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_context_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    
      h2("A Random Text 1"),
      tableOutput(ns("text2"))
      
 
  )
}
    
#' context Server Function
#'
#' @noRd 
mod_context_server <- function(input, output, session){
  ns <- session$ns
  output$text2 <- renderText({
    random_text(nwords = 50)

  })
}
    
## To be copied in the UI
# mod_context_ui("context_ui_1")
    
## To be copied in the server
# callModule(mod_context_server, "context_ui_1")
 
