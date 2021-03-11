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
      shiny::includeMarkdown(system.file("extdata", "welcome.md", package = "DIANE")),
      shinydashboardPlus::socialButton(
        href = "https://github.com/OceaneCsn/DIANE",
        icon = shiny::icon("github")
      )
 
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
 
