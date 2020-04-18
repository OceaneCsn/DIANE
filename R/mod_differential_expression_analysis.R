#' differential_expression_analysis UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_differential_expression_analysis_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' differential_expression_analysis Server Function
#'
#' @noRd 
mod_differential_expression_analysis_server <- function(input, output, session){
  ns <- session$ns
 
}
    
## To be copied in the UI
# mod_differential_expression_analysis_ui("differential_expression_analysis_ui_1")
    
## To be copied in the server
# callModule(mod_differential_expression_analysis_server, "differential_expression_analysis_ui_1")
 
