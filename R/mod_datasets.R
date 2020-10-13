#' datasets UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_datasets_ui <- function(id){
  ns <- NS(id)
  tagList(
    
   shiny::uiOutput(ns("examples"))
 
  )
}
    
#' datasets Server Function
#'
#' @noRd 
mod_datasets_server <- function(input, output, session){
  ns <- session$ns
  
  
  output$examples <- shiny::renderUI({
    
    if(golem::get_golem_options("server_version"))
      shiny::includeMarkdown(
        system.file("extdata", "ex_datasets.md", package = "DIANE"))
    else
      h4("Dataset examples on human, mouse, drosophilia and more are 
         available for download in this tab at diane.bpmp.inrae.fr")
    
  })
 
}
    
## To be copied in the UI
# mod_datasets_ui("datasets_ui_1")
    
## To be copied in the server
# callModule(mod_datasets_server, "datasets_ui_1")
 
