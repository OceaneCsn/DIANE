#' cluster_exploration UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_cluster_exploration_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    shiny::h1("Analyse the genes of a specific cluster"),
    shiny::hr(),
    shiny::uiOutput(ns("cluster_to_explore_choice")),
    shiny::hr(),
    
#   ____________________________________________________________________________
#   Profiles column                                                         ####

  col_4(
    boxPlus( width = 12, title = "Expression profiles",
      shiny::plotOutput(ns("profiles_to_explore"), height = "700px"))
  ),


#   ____________________________________________________________________________
#   Clusters characteristics                                                ####


  col_8(
    shinydashboard::tabBox(title = "Genes in that clusters", width = 12,
           shiny::tabPanel(title = "Genes table",
                           DT::dataTableOutput(ns("genes_to_explore"))),
           shiny::tabPanel(title = "Ontologies enrichment"),
           shiny::tabPanel(title = "GLM for factors effect")
    )
  )
 
  )
}
    
#' cluster_exploration Server Function
#'
#' @noRd 
mod_cluster_exploration_server <- function(input, output, session, r){
  ns <- session$ns
  
  membership <- shiny::reactive({
    req(r$clusterings[[r$current_comparison]])
    r$clusterings[[r$current_comparison]]$membership
  })
    
  output$cluster_to_explore_choice <- shiny::renderUI({
    shiny::req(membership())
    shinyWidgets::radioGroupButtons(
      inputId = ns("cluster_to_explore"),
      label = "Cluster to explore",
      choices = unique(membership()),
      justified = TRUE,
      checkIcon = list(
        yes = shiny::icon("ok", 
                          lib = "glyphicon"))
    )
  })
  
  output$profiles_to_explore <- shiny::renderPlot({
    draw_profiles(data = r$normalized_counts,
                  membership = membership(),
                  k = input$cluster_to_explore)
    })
  
  output$genes_to_explore <- DT::renderDataTable({
    req(r$top_tags, r$current_comparison, membership())
    
    table <- r$top_tags[[r$current_comparison]]
    table[table$genes %in% get_genes_in_cluster(membership = membership(), cluster = input$cluster_to_explore),
          c("logFC", "logCPM", "FDR")]
  })
  
 
}
    
## To be copied in the UI
# mod_cluster_exploration_ui("cluster_exploration_ui_1")
    
## To be copied in the server
# callModule(mod_cluster_exploration_server, "cluster_exploration_ui_1")
 
