#' network_analysis UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_network_analysis_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::h1("Network analysis")
 
  )
}
    
#' network_analysis Server Function
#'
#' @noRd 
mod_network_analysis_server <- function(input, output, session){
  ns <- session$ns
 
}
    
## To be copied in the UI
# mod_network_analysis_ui("network_analysis_ui_1")
    
## To be copied in the server
# callModule(mod_network_analysis_server, "network_analysis_ui_1")
 
