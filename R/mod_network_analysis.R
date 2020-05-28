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
    # TODO group selection in network
    
    # TODO gene zoom in network
    
    #   ____________________________________________________________________________
    #   network view                                                            ####
    
    shiny::fluidRow(col_2(shinyWidgets::switchInput(
      inputId = ns("louvain_color_switch"),
      value = FALSE,
      onLabel = "Communities",
      offLabel = "Gene type",
      label = "Nodes color"
    )),
    
    col_6(shiny::uiOutput(ns("cluster_to_explore_choice"))),
    col_4(shiny::uiOutput(ns("network_summary")))),
    
    shiny::hr(),
    
    

  column(width = 5,
         shiny::textInput(ns("gene_to_zoom"), label = "Gene ID to focus on :"),
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
        shiny::h4("Topolocical clusters correspond to the network structural communitities"),
        shiny::plotOutput(ns("profiles"), height = "750px")
        
      ),
      shiny::tabPanel(
        title = "Per module GO enrichment",
        col_4(shinyWidgets::actionBttn(
          ns("go_enrich_btn"),
          label = "Start GO enrichment analysis for this gene community",
          color = "success",
          style = 'bordered'
        )),
        
        col_4(shinyWidgets::switchInput(
          inputId = ns("draw_go"),
          value = TRUE,
          onLabel = "Plot",
          offLabel = "Data table",
          label = "Result type"
        )),
        col_4(shiny::uiOutput(ns("max_go_choice"))),
        
        shiny::hr(),
        
        shiny::fluidRow(col_12(shiny::uiOutput(ns("go_results"))))
        
      )
      
    )
  )

 
  )
}
    
#' network_analysis Server Function
#' @import visNetwork
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
                 edges = r$networks[[r$current_network]]$edges) %>% 
      visNetwork::visOptions(
                   highlightNearest = list(enabled = TRUE,
                                           degree = 0),
                   collapse = FALSE,
                  ) %>%
                 visEvents(click = "function(nodes){
                  Shiny.onInputChange('network_analysis_ui_1-click', nodes.nodes);
                  ;}") 
  })
  
  output$module <- shiny::renderPrint({
    print(paste(input$click, input$select))
  })
  
  
  
#   ____________________________________________________________________________
#   zoom                                                                    ####

  shiny::observe({
    
    shiny::req(r$current_network, r$networks)
    shiny::req(r$networks[[r$current_network]]$nodes)
    
    nodes <- r$networks[[r$current_network]]$nodes
    
    visNetwork::visNetworkProxy(ns("network_view")) %>% 
      visNetwork::visFocus(id = input$gene_to_zoom, scale = 1) %>% 
      visNetwork::visSelectNodes(id = input$gene_to_zoom)
    
  })
  
  #   ____________________________________________________________________________
  #   module  selection                                                       ####
  
  output$cluster_to_explore_choice <- shiny::renderUI({
    shiny::req(r$current_network, r$networks)
    shiny::req(r$networks[[r$current_network]]$membership)
    
    shinyWidgets::radioGroupButtons(
      inputId = ns("cluster_to_explore"),
      label = "Cluster to explore",
      choices = c("All", unique(r$networks[[r$current_network]]$membership)),
      justified = TRUE, selected = "All",
      checkIcon = list(yes = shiny::icon("ok",
                                         lib = "glyphicon"))
    )
  })
  
  shiny::observe({
    
    shiny::req(r$current_network, r$networks)
    shiny::req(r$networks[[r$current_network]]$nodes)
    shiny::req(input$cluster_to_explore)
    shiny::req(input$cluster_to_explore != "All")
    
    nodes <- r$networks[[r$current_network]]$nodes
    
    visNetwork::visNetworkProxy(ns("network_view")) %>% 
      visNetwork::visSelectNodes(id = get_genes_in_cluster(r$networks[[r$current_network]]$membership, 
                                                           cluster = input$cluster_to_explore))
    
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
    data <- data[order(-data$degree),]
    if(input$cluster_to_explore == "All")
      data[, columns]
    else
      data[data$community == input$cluster_to_explore, columns]
  })
  
  

#   ____________________________________________________________________________
#   densities                                                               ####

  output$distributions <- shiny::renderPlot({
    
    shiny::req(r$current_network, r$networks)
    shiny::req(r$networks[[r$current_network]]$nodes)
    
    draw_network_degrees(nodes = r$networks[[r$current_network]]$nodes,
                         graph = r$networks[[r$current_network]]$graph)
  })
  
  
  
#   ____________________________________________________________________________
#   profiles                                                                ####

  output$profiles <- shiny::renderPlot({
    shiny::req(r$normalized_counts, r$networks)
    shiny::req(r$networks[[r$current_network]]$membership)
    shiny::req(r$networks[[r$current_network]]$conditions)
    shiny::req(input$cluster_to_explore)
    
    if(input$cluster_to_explore == "All"){
      draw_profiles(data = r$normalized_counts,
                    membership = r$networks[[r$current_network]]$membership,
                    conds = r$networks[[r$current_network]]$conditions)
    }
    else{
      draw_profiles(data = r$normalized_counts,
                    membership = r$networks[[r$current_network]]$membership,
                    conds = r$networks[[r$current_network]]$conditions, 
                    k = input$cluster_to_explore)
    }
    
  })
  
  r_mod <- shiny::reactiveValues(go = NULL)
  #   ____________________________________________________________________________
  #   GO enrich                                                               ####
  
  
  shiny::observeEvent((input$go_enrich_btn), {
    shiny::req(r$normalized_counts)
    shiny::req(r$networks[[r$current_network]]$membership)
    
    if (r$organism == "Other") {
      shinyalert::shinyalert("For now, only Arabidopsis thaliana and 
        Homo sapiens are supported for GO analysis", 
                             "Did you correctly set your organism in the 
                               Data import tab?",
                             type = "error")
    }
    shiny::req(r$organism != "Other")
    
    
    if (input$cluster_to_explore == "All") {
      shinyalert::shinyalert("Please specify a module to perform the analysis on", 
                            type = "error")
    }
    shiny::req(r$organism != "Other")
    
    shiny::req(input$cluster_to_explore != "All")
    
    
    
    
    genes <- get_genes_in_cluster(membership = 
                                    r$networks[[r$current_network]]$membership,
                                  cluster = input$cluster_to_explore)
    
    background <- rownames(r$normalized_counts)
    
    if (r$splicing_aware){
      genes <- get_locus(genes)
      background <- get_locus(background)
    }
    
    if (r$organism == "Arabidopsis thaliana") {
      genes <- convert_from_agi(genes)
      background <- convert_from_agi(background)
      org = org.At.tair.db::org.At.tair.db
    }
    
    if(r$organism == "Homo sapiens"){
      genes <- convert_from_ensembl(genes)
      background <- convert_from_ensembl(background)
      org = org.Hs.eg.db::org.Hs.eg.db
    }
    
    # TODO add check if it is entrez with regular expression here
    shiny::req(length(genes) > 0, length(background) > 0)
    
    r_mod$go <- enrich_go(genes, background, org = org)
  })
  
  #   ____________________________________________________________________________
  #   go results                                                              ####
  
  output$go_table <- DT::renderDataTable({
    shiny::req(r_mod$go)
    r_mod$go[,c("Description", "GeneRatio", "BgRatio", "p.adjust")]
  })
  
  output$max_go_choice <- shiny::renderUI({
    shiny::req(r_mod$go)
    shiny::numericInput(ns("n_go_terms"), 
                        label = "Top number of GO terms to plot :", 
                        min = 1, value = dim(r_mod$go)[1])
  })
  
  output$go_plot <- plotly::renderPlotly({
    shiny::req(r_mod$go)
    max = ifelse(is.na(input$n_go_terms), dim(r_mod$go)[1],input$n_go_terms )
    draw_enrich_go(r_mod$go, max_go = max)
  })
  
  output$go_results <- shiny::renderUI({
    
    if(r$organism == "Other")
      shiny::h4("GO analysis is only supported for Arabidopsis and human (for now!)")
    
    shiny::req(r$organism != "Other")
    
    shiny::req(r_mod$go)
    if (!input$draw_go){
      DT::dataTableOutput(ns("go_table"))
    }
    else{
      plotly::plotlyOutput(ns("go_plot"), height = "700px")
    }
  })

}
    
## To be copied in the UI
# mod_network_analysis_ui("network_analysis_ui_1")
    
## To be copied in the server
# callModule(mod_network_analysis_server, "network_analysis_ui_1")
 
