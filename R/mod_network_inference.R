#' network_inference UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_network_inference_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::h1("Network inference")
 
  )
}
    
#' network_inference Server Function
#'
#' @noRd 
mod_network_inference_server <- function(input, output, session){
  ns <- session$ns
 
}
    
## To be copied in the UI
# mod_network_inference_ui("network_inference_ui_1")
    
## To be copied in the server
# callModule(mod_network_inference_server, "network_inference_ui_1")
 
