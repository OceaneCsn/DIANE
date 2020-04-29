#' context UI Function
#'
#' @description A shiny Module for the application context.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_context_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    
      h2("Dashboard for the Inference and Analysis of Networks from Expression data"),

      tags$img(src="/www/hex-DIANE.png", width = 300)
 
  )
}
    
#' context Server Function
#'
#' @noRd 
mod_context_server <- function(input, output, session){
  ns <- session$ns
  
}
    
## To be copied in the UI
# mod_context_ui("context_ui_1")
    
## To be copied in the server
# callModule(mod_context_server, "context_ui_1")
 
