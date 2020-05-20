#' module_analysis UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_module_analysis_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::h1("Communities analysis")
 
  )
}
    
#' module_analysis Server Function
#'
#' @noRd 
mod_module_analysis_server <- function(input, output, session){
  ns <- session$ns
 
}
    
## To be copied in the UI
# mod_module_analysis_ui("module_analysis_ui_1")
    
## To be copied in the server
# callModule(mod_module_analysis_server, "module_analysis_ui_1")
 
