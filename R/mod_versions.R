#' versions UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_versions_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::includeMarkdown(system.file("extdata", "NEWS.md", package = "DIANE"))
  )
}
    
#' versions Server Functions
#'
#' @noRd 
mod_versions_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_versions_ui("versions_ui_1")
    
## To be copied in the server
# mod_versions_server("versions_ui_1")
 