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
    
    
    #   ____________________________________________________________________________
    #   network view                                                            ####
    
    shiny::fluidRow(col_2(shinyWidgets::switchInput(
      inputId = ns("louvain_color_switch"),
      value = FALSE,
      onLabel = "Communities",
      offLabel = "Gene type",
      label = "Nodes color"
    )),
    
    col_10(shiny::uiOutput(ns("network_summary")))),
    
    shiny::hr(),

  column(width = 5,
    visNetwork::visNetworkOutput(ns("network_view"), height = "900px")
  ),
  
  
  
#   ____________________________________________________________________________
#   network infos                                                           ####

  column(width = 7,
    shinydashboard::tabBox(
      title = "Network informations",
      width = 12,
      
      
      shiny::tabPanel(
        title = "Gene ranking",
        DT::dataTableOutput(ns("gene_ranking"))
        
      ),
      shiny::tabPanel(
        title = "Degree distributions",
        shiny::plotOutput(ns("distributions"), height = "700px")
        
      ),
      shiny::tabPanel(
        title = "Per module expression profiles",
        shiny::h5("Topolocical clusters expression profiles :"),
        shiny::plotOutput(ns("profiles"), height = "750px")
        
      )
      
    )
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
  
  
  
#   ____________________________________________________________________________
#   network summaries                                                       ####

  output$network_summary <- shiny::renderUI({
    shiny::req(r$current_network, r$networks)
    shiny::req(r$networks[[r$current_network]]$nodes)
    shiny::req(r$networks[[r$current_network]]$edges)
    
    
    nodes <- r$networks[[r$current_network]]$nodes
    n_genes <- dim(nodes)[1]
    n_tfs <- dim(nodes[nodes$gene_type == "Regulator",])[1]
    n_edges <- dim(r$networks[[r$current_network]]$edges)[1]
    tagList(
      shiny::fluidRow(
        col_4(shinydashboardPlus::descriptionBlock(
          number = n_genes,
          number_color = "primary",
          text = "Genes",
          right_border = TRUE
        )),
        col_4(shinydashboardPlus::descriptionBlock(
          number = n_tfs,
          number_color = "green",
          text = "Regulators",
          right_border = TRUE
        )),
        col_4(shinydashboardPlus::descriptionBlock(
          number = n_edges,
          number_color = "navy",
          text = "Edges",
          right_border = FALSE
        )
      ))
    )
    
  })
  
  
#   ____________________________________________________________________________
#   table                                                                   ####

  
  output$gene_ranking <- DT::renderDataTable({
    
    shiny::req(r$current_network, r$networks)
    shiny::req(r$networks[[r$current_network]]$nodes)
    
    data <- r$networks[[r$current_network]]$nodes
    
    columns <- c("label", "gene_type", "degree", "community")
    if (!is.null(r$gene_info)) {
      columns <- unique(c(colnames(r$gene_info), columns))
    }
    
    data[order(-data$degree), columns]
  })
  
  

#   ____________________________________________________________________________
#   densities                                                               ####

  output$distributions <- shiny::renderPlot({
    
    shiny::req(r$current_network, r$networks)
    shiny::req(r$networks[[r$current_network]]$nodes)
    
    draw_network_degrees(nodes = r$networks[[r$current_network]]$nodes,
                         graph = r$networks[[r$current_network]]$graph)
  })
  
  
  output$profiles <- shiny::renderPlot({
    shiny::req(r$normalized_counts, r$networks)
    shiny::req(r$networks[[r$current_network]]$membership)
    shiny::req(r$networks[[r$current_network]]$conditions)
    
    draw_profiles(data = r$normalized_counts,
                  membership = r$networks[[r$current_network]]$membership,
                  conds = r$networks[[r$current_network]]$conditions)
  })

}
    
## To be copied in the UI
# mod_network_analysis_ui("network_analysis_ui_1")
    
## To be copied in the server
# callModule(mod_network_analysis_server, "network_analysis_ui_1")
 
