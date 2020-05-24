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
    shiny::h1("Network analysis"),
    
    shiny::hr(),
    
    shinyalert::useShinyalert(),
    
    shinybusy::add_busy_spinner(
      spin = "self-building-square",
      position = 'top-left',
      margins = c(70, 1200)
    ),
    
    shinyWidgets::switchInput(
      inputId = ns("louvain_color_switch"),
      value = FALSE,
      onLabel = "Communities",
      offLabel = "Gene type",
      label = "Nodes color"
    ),
    
    
    
#   ____________________________________________________________________________
#   network view                                                            ####
  col_6(
    visNetwork::visNetworkOutput(ns("network_view"), height = "900px")
  )
 
  )
}
    
#' network_analysis Server Function
#'
#' @noRd 
mod_network_analysis_server <- function(input, output, session, r){
  ns <- session$ns
  
  
#   ____________________________________________________________________________
#   network view                                                            ####

  output$network_view <- visNetwork::renderVisNetwork({
    
    shiny::req(r$current_network, r$networks)
    shiny::req(r$networks[[r$current_network]]$nodes)
    shiny::req(r$networks[[r$current_network]]$edges)

    nodes <- r$networks[[r$current_network]]$nodes
    
    if(input$louvain_color_switch ) nodes$group <- nodes$community
    else nodes$group <- nodes$gene_type
    
    
    draw_network(nodes = nodes,
                 edges = r$networks[[r$current_network]]$edges)
  })

}
    
## To be copied in the UI
# mod_network_analysis_ui("network_analysis_ui_1")
    
## To be copied in the server
# callModule(mod_network_analysis_server, "network_analysis_ui_1")
 
