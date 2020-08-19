#' network_analysis UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_network_analysis_ui <- function(id) {
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
    
    shiny::fluidRow(
      col_2(
        shinyWidgets::switchInput(
          inputId = ns("louvain_color_switch"),
          value = FALSE,
          onLabel = "Communities",
          offLabel = "Gene type",
          label = "Nodes color"
        )
      ),
      
      col_6(shiny::uiOutput(ns(
        "cluster_to_explore_choice"
      ))),
      col_4(shiny::uiOutput(ns(
        "network_summary"
      )))
    ),
    
    shiny::hr(),
    
    
    
    column(
      width = 5,
      shiny::textInput(ns("gene_to_zoom"), label = "Gene ID to focus on :"),
      visNetwork::visNetworkOutput(ns("network_view"), height = "900px")
    ),
    
    
    
    #   ____________________________________________________________________________
    #   network infos                                                           ####
    
    column(
      width = 7,
      shinydashboard::tabBox(
        width = 12,
        
        
        shiny::tabPanel(
          title = "Degree-ranked gene list",
          DT::dataTableOutput(ns("gene_ranking")),
          shinyWidgets::downloadBttn(
            outputId = ns("download_node_table"),
            label = "Download nodes as csv table",
            style = "bordered",
            color = "success"
          ),
          shinyWidgets::downloadBttn(
            outputId = ns("download_edges_table"),
            label = "Download edges as csv table",
            style = "bordered",
            color = "success"
          )
          
        ),
        shiny::tabPanel(
          title = "Correlated regulators network",
          shiny::h5(
            "The super nodes (dark green squares) of the network are detailed here.
        Each edge represents a correlation above the specified threshold between
        the regulators. The community detection in this graph was used to group
        highly correlated variables before network inference with GENIE3."
          ),
          visNetwork::visNetworkOutput(ns("cor_tfs_network"), height = "700px")
          
        ),
        shiny::tabPanel(title = "In-Out degree distributions",
                        shiny::plotOutput(ns("distributions"), height = "700px")),
        shiny::tabPanel(
          title = "Modules expression profiles",
          shiny::h4(
            "Topolocical clusters correspond to the network structural communitities."
          ),
          shiny::h5(
            "The expression profiles of genes within a community can be positively correlated,
                  but also, unlike with expression based clustering, negatively correlated."
          ),
          shiny::plotOutput(ns("profiles"), height = "750px")
          
        ),
        shiny::tabPanel(
          title = "Modules GO enrichment",
          col_4(
            shinyWidgets::actionBttn(
              ns("go_enrich_btn"),
              label = "Start GO enrichment analysis for this gene community",
              color = "success",
              style = 'bordered'
            )
          ),
          
          col_4(
            shinyWidgets::radioGroupButtons(
              ns("draw_go"),
              choices = c("Dot plot", "Enrichment map", "Data table"),
              selected = "Dot plot",
              justified = TRUE,
              direction = "vertical",
              checkIcon = list(yes = icon("ok",
                                          lib = "glyphicon"))
            )
          ),
          col_4(
            shinyWidgets::radioGroupButtons(
              ns("go_type"),
              choiceNames = c(
                "Biological process",
                "Cellular component",
                "Molecular function"
              ),
              choiceValues = c("BP", "CC", "MF"),
              selected = "BP",
              justified = TRUE,
              direction = "vertical",
              checkIcon = list(yes = icon("ok",
                                          lib = "glyphicon"))
            ),
            shiny::uiOutput(ns("max_go_choice"))
          ),
          
          shiny::hr(),
          
          shiny::fluidRow(col_12(shiny::uiOutput(ns(
            "go_results"
          ))))
          
        )
        
      )
    )
  )
}

#' network_analysis Server Function
#' @import visNetwork
#' @noRd
mod_network_analysis_server <- function(input, output, session, r) {
  ns <- session$ns
  
  
  #   ____________________________________________________________________________
  #   network view                                                            ####
  
  output$network_view <- visNetwork::renderVisNetwork({
    shiny::req(r$current_network, r$networks)
    shiny::req(r$networks[[r$current_network]]$nodes)
    shiny::req(r$networks[[r$current_network]]$edges)
    
    nodes <- r$networks[[r$current_network]]$nodes
    
    if (input$louvain_color_switch)
      nodes$group <- nodes$community
    else
      nodes$group <- nodes$gene_type
    
    
    draw_network(nodes = nodes,
                 edges = r$networks[[r$current_network]]$edges) %>%
      visNetwork::visOptions(highlightNearest = list(enabled = TRUE,
                                                     degree = 0),
                             collapse = FALSE,) %>%
      visEvents(click = "function(nodes){
                  Shiny.onInputChange('network_analysis_ui_1-click', nodes.nodes);
                  ;}")
  })
  
  output$module <- shiny::renderPrint({
    print(paste(input$click, input$select))
  })
  
  
  
  #   ____________________________________________________________________________
  #   Node description                                                        ####
  
  shiny::observeEvent(input$click, {
    shiny::req(r$current_network, r$networks)
    shiny::req(r$networks[[r$current_network]]$nodes)
    data <- r$networks[[r$current_network]]$graph
    
    shiny::showModal(
      shiny::modalDialog(
        title = "Gene description",
        size = 'l',
        shiny::htmlOutput(ns("node_details")),
        
        shiny::hr(),
        
        shiny::plotOutput(ns("node_profile")),
        
        shiny::hr(),
        
        shiny::uiOutput(ns("node_reg")),
        
        shiny::hr(),
        
        shiny::uiOutput(ns("node_targ")),
        
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  
  output$node_profile <- shiny::renderPlot({
    if (r$splicing_aware) {
      data <- r$aggregated_normalized_counts
    }
    else{
      data <- r$normalized_counts
    }
    
    if (sum(grepl("mean_",
                  r$networks[[r$current_network]]$nodes$id)) > 0) {
      data <- r$grouped_normalized_counts
    }
    
    
    draw_expression_levels(data,
                           genes = c(input$click),
                           conds = r$networks[[r$current_network]]$conditions)
  })
  
  node_descr <- shiny::reactive({
    shiny::req(r$current_network, r$networks)
    shiny::req(r$networks[[r$current_network]]$nodes)
    return(describe_node(r$networks[[r$current_network]]$graph, input$click))
  })
  
  output$node_reg <- shiny::renderUI({
    if (length(node_descr()$regulators) > 0) {
      tagList(shiny::h3("Regulators :"),
              DT::dataTableOutput(ns("node_regulators")))
    }
    else
      shiny::h3("No regulators found")
  })
  
  output$node_targ <- shiny::renderUI({
    if (length(node_descr()$targets) > 0) {
      tagList(shiny::h3("Targets :"),
              DT::dataTableOutput(ns("node_targets")))
    }
    else
      shiny::h3("No targets found")
  })
  
  
  output$node_regulators <- DT::renderDataTable({
    columns <- c("label", "gene_type", "degree", "community")
    if (!is.null(r$gene_info)) {
      columns <- unique(c(colnames(r$gene_info), columns))
    }
    regulators <- node_descr()$regulators
    
    data <- r$networks[[r$current_network]]$nodes
    data[regulators, columns]
    
  })
  
  output$node_targets <- DT::renderDataTable({
    columns <- c("label", "gene_type", "degree", "community")
    if (!is.null(r$gene_info)) {
      columns <- unique(c(colnames(r$gene_info), columns))
    }
    targets <- node_descr()$targets
    
    data <- r$networks[[r$current_network]]$nodes
    data[targets, columns]
  })
  
  output$node_details <- shiny::renderText({
    shiny::req(r$current_network, r$networks)
    shiny::req(r$networks[[r$current_network]]$nodes)
    data <- r$networks[[r$current_network]]$nodes
    
    if ("label" %in% colnames(data))
      label <- data[input$click, "label"]
    else
      label <- "-"
    
    if ("description" %in% colnames(data)) {
      if (stringr::str_detect(input$click, "mean_")) {
        tfs <-
          unlist(strsplit(stringr::str_remove(input$click, 'mean_'), '-'))
        if (!is.null(r$gene_info) &
            "description" %in% colnames(r$gene_info)) {
          description <- ""
          for (tf in tfs) {
            description <-
              paste(description, '<br>', tf, ':', r$gene_info[tf, "description"])
          }
        }
        else
          description <- "-"
      }
      else{
        description <- data[input$click, "description"]
      }
    }
    else
      description <- "-"
    
    descr <- paste(
      "<b> AGI : </b>",
      input$click,
      '<br>',
      "<b> Common name : </b>",
      label,
      '<br>',
      "<b> Description : </b>",
      description
    )
    descr
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
      justified = TRUE,
      selected = "All",
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
    n_tfs <- dim(nodes[nodes$gene_type == "Regulator", ])[1]
    n_edges <- dim(r$networks[[r$current_network]]$edges)[1]
    tagList(shiny::fluidRow(
      col_4(
        shinydashboardPlus::descriptionBlock(
          number = n_genes,
          numberColor = "primary",
          text = "Genes",
          rightBorder = TRUE
        )
      ),
      col_4(
        shinydashboardPlus::descriptionBlock(
          number = n_tfs,
          numberColor = "green",
          text = "Regulators",
          rightBorder = TRUE
        )
      ),
      col_4(
        shinydashboardPlus::descriptionBlock(
          number = n_edges,
          numberColor = "navy",
          text = "Edges",
          rightBorder = FALSE
        )
      )
    ))
    
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
    data <- data[order(-data$degree), ]
    if (input$cluster_to_explore == "All")
      data[, columns]
    else
      data[data$community == input$cluster_to_explore, columns]
  })
  
  
  #   ____________________________________________________________________________
  #   TFs network                                                             ####
  
  
  output$cor_tfs_network <- visNetwork::renderVisNetwork({
    shiny::req(r$current_network, r$networks)
    shiny::req(r$networks[[r$current_network]]$nodes)
    shiny::req(r$networks[[r$current_network]]$edges)
    shiny::req(r$cor_network)
    
    nodes <- r$cor_network$nodes
    nodes$label <-
      r$gene_info[match(nodes$id, rownames(r$gene_info)), "label"]
    
    visNetwork::visNetwork(nodes, r$cor_network$edges) visNetwork::%>% 
      visNetwork::visNodes(font = list("size" = 35), 
               color = list("background" = "#1C5435", "border" = "#FFFFCC"),
               shape = "square", size = 30)
  })
  
  output$module <- shiny::renderPrint({
    print(paste(input$click, input$select))
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
    
    if (r$splicing_aware) {
      data <- r$aggregated_normalized_counts
    }
    else{
      data <- r$normalized_counts
    }
    
    if (sum(grepl("mean_",
                  r$networks[[r$current_network]]$nodes$id)) > 0) {
      data <- r$grouped_normalized_counts
    }
    
    if (input$cluster_to_explore == "All") {
      draw_profiles(
        data = data,
        membership = r$networks[[r$current_network]]$membership,
        conds = r$networks[[r$current_network]]$conditions
      )
    }
    else{
      draw_profiles(
        data = data,
        membership = r$networks[[r$current_network]]$membership,
        conds = r$networks[[r$current_network]]$conditions,
        k = input$cluster_to_explore
      )
    }
    
  })
  
  r_mod <- shiny::reactiveValues(go = NULL)
  #   ____________________________________________________________________________
  #   GO enrich                                                               ####
  
  
  shiny::observeEvent((input$go_enrich_btn), {
    shiny::req(r$normalized_counts)
    shiny::req(r$networks[[r$current_network]]$membership)
    
    if (r$organism == "Other") {
      shinyalert::shinyalert(
        "For now, only Arabidopsis thaliana and
        Homo sapiens are supported for GO analysis",
        "Did you correctly set your organism in the
                               Data import tab?",
        type = "error"
      )
    }
    shiny::req(r$organism != "Other")
    
    
    if (input$cluster_to_explore == "All") {
      shinyalert::shinyalert("Please specify a module to perform the analysis on",
                             type = "error")
    }
    shiny::req(r$organism != "Other")
    
    shiny::req(input$cluster_to_explore != "All")
    
    
    
    
    genes <- get_genes_in_cluster(
      membership =
        r$networks[[r$current_network]]$membership,
      cluster = input$cluster_to_explore
    )
    
    background <- rownames(r$normalized_counts)
    
    # spreads the grouped regulators
    if (sum(grepl("mean_", genes)) > 0) {
      individuals <- genes[!grepl("means_", genes)]
      groups <- setdiff(genes, individuals)
      for (group in groups) {
        individuals <- c(individuals,
                         strsplit(stringr::str_split_fixed(group, "_", 2)[, 2]),
                         '-')
      }
      genes <- individuals
    }
    
    
    
    if (r$splicing_aware) {
      genes <- get_locus(genes)
      background <- get_locus(background)
    }
    
    if (r$organism == "Arabidopsis thaliana") {
      genes <- convert_from_agi(genes)
      background <- convert_from_agi(background)
      org = org.At.tair.db::org.At.tair.db
    }
    
    if (r$organism == "Homo sapiens") {
      genes <- convert_from_ensembl(genes)
      background <- convert_from_ensembl(background)
      org = org.Hs.eg.db::org.Hs.eg.db
    }
    
    if (r$organism == "Mus musculus") {
      genes <- convert_from_ensembl_mus(genes)
      background <- convert_from_ensembl_mus(background)
      org = org.Mm.eg.db::org.Mm.eg.db
    }
    
    # TODO add check if it is entrez with regular expression here
    shiny::req(length(genes) > 0, length(background) > 0)
    
    r_mod$go <-
      enrich_go(genes,
                background,
                org = org,
                GO_type = input$go_type)
  })
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### download GO                                                             ####
  
  
  output$download_go_table <- shiny::downloadHandler(
    filename = function() {
      paste(paste0(
        "enriched_GOterms_module_",
        input$cluster_to_explore,
        ".csv"
      ))
    },
    content = function(file) {
      write.csv(r_mod$go, file = file, quote = FALSE)
    }
  )
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### download nodes                                                          ####
  
  
  output$download_node_table <- shiny::downloadHandler(
    filename = function() {
      paste("network_nodes.csv")
    },
    content = function(file) {
      write.csv(r$networks[[r$current_network]]$nodes, file = file, quote = FALSE)
    }
  )
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### download edges                                                          ####
  
  
  output$download_edges_table <- shiny::downloadHandler(
    filename = function() {
      paste("network_edges.csv")
    },
    content = function(file) {
      write.csv(r$networks[[r$current_network]]$edges, file = file, quote = FALSE)
    }
  )
  
  #   ____________________________________________________________________________
  #   go results                                                              ####
  
  output$go_table <- DT::renderDataTable({
    shiny::req(r_mod$go)
    r_mod$go[, c("Description", "GeneRatio", "BgRatio", "p.adjust")]
  })
  
  output$max_go_choice <- shiny::renderUI({
    shiny::req(r_mod$go)
    shiny::numericInput(
      ns("n_go_terms"),
      label = "Top number of GO terms to plot :",
      min = 1,
      value = dim(r_mod$go)[1]
    )
  })
  
  output$go_plot <- plotly::renderPlotly({
    shiny::req(r_mod$go)
    max = ifelse(is.na(input$n_go_terms),
                 dim(r_mod$go)[1],
                 input$n_go_terms)
    draw_enrich_go(r_mod$go, max_go = max)
  })
  
  output$go_map_plot <- shiny::renderPlot({
    shiny::req(r_mod$go)
    draw_enrich_go_map(r_mod$go)
  })
  
  output$go_results <- shiny::renderUI({
    if (r$organism == "Other")
      shiny::h4("GO analysis is only supported for Arabidopsis and human (for now!)")
    
    shiny::req(r$organism != "Other")
    
    shiny::req(r_mod$go)
    if (nrow(r_mod$go) == 0) {
      shinyalert::shinyalert(
        "No enriched GO terms were found",
        "It can happen if input gene list is not big enough",
        type = "error"
      )
    }
    
    shiny::req(nrow(r_mod$go) > 0)
    
    
    if (input$draw_go == "Data table") {
      tagList(
        DT::dataTableOutput(ns("go_table")),
        shinyWidgets::downloadBttn(
          outputId = ns("download_go_table"),
          label = "Download enriched GO term as a csv table",
          style = "bordered",
          color = "success"
        )
      )
    }
    else{
      if (input$draw_go == "Enrichment map") {
        shiny::plotOutput(ns("go_map_plot"), height = "800px")
      }
      else
        plotly::plotlyOutput(ns("go_plot"), height = "800px")
    }
  })
  
}

## To be copied in the UI
# mod_network_analysis_ui("network_analysis_ui_1")

## To be copied in the server
# callModule(mod_network_analysis_server, "network_analysis_ui_1")
