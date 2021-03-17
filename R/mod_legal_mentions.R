#' legal_mentions UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_legal_mentions_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::includeMarkdown(system.file("extdata", "legal_mentions.md", package = "DIANE"))
  )
}
    
#' legal_mentions Server Function
#'
#' @noRd 
mod_legal_mentions_server <- function(input, output, session){
  ns <- session$ns
 
}
    
## To be copied in the UI
# mod_legal_mentions_ui("legal_mentions_ui_1")
    
## To be copied in the server
# callModule(mod_legal_mentions_server, "legal_mentions_ui_1")
 
