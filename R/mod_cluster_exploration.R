#' cluster_exploration UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_cluster_exploration_ui <- function(id) {
  ns <- NS(id)
  tagList(
    shinybusy::add_busy_spinner(
      spin = "self-building-square",
      position = 'top-left',
      margins = c(70, 1100)
    ),
    
    shiny::h1("Analyse the genes of a specific cluster"),
    shiny::hr(),
    shiny::fluidRow(col_8(shiny::uiOutput(ns("cluster_to_explore_choice"))),
                    col_4(valueBoxOutput(ns(
                      "gene_number_cluster"
                    )))
    ),
    shiny::hr(),
    
    #   ____________________________________________________________________________
    #   Profiles column                                                         ####
    
    col_6(
      boxPlus(
        width = 12,
        closable = FALSE,
        title = "Expression profiles",
        shiny::plotOutput(ns("profiles_to_explore"), height = "700px")
      )
    ),
    
    # TODO download data from one cluster
    
    # TODO true reset when toggle demo data
    #   ____________________________________________________________________________
    #   Clusters characteristics                                                ####
    
    
    col_6(
      shinydashboard::tabBox(
        title = "Genes in that clusters",
        width = 12,
        shiny::tabPanel(
          title = "Genes table",
          DT::dataTableOutput(ns("genes_to_explore"))
          
        ),
        
        
#   ____________________________________________________________________________
#   GO                                                                      ####

        shiny::tabPanel(title = "Gene Ontologies enrichment",
                        
                        
                        col_6(shinyWidgets::actionBttn(
                          ns("go_enrich_btn"),
                          label = "Start GO enrichment analysis",
                          color = "success",
                          style = 'bordered'
                        )),
                        
                        col_6(shinyWidgets::switchInput(
                          inputId = ns("draw_go"),
                          value = TRUE,
                          onLabel = "Plot",
                          offLabel = "Data table",
                          label = "Result type"
                        )),
                        
                        shiny::hr(),
                        
                        shiny::fluidRow(col_12(shiny::uiOutput(ns("go_results"))))
                        ),
        
#   ____________________________________________________________________________
#   glm                                                                     ####

        shiny::tabPanel(title = "GLM for factors effect",
                        col_2(shinyWidgets::dropdownButton(
                          size = 'xs',
                          shiny::verbatimTextOutput(ns("glm_summary")),
                          circle = FALSE,
                          status = "success",
                          width = "600px",
                          label = "Glm summary")
                        ),
                        shiny::plotOutput(ns("glm_plot")),
                        shiny::hr(),
                        shiny::h5("The absolute value of a coefficient gives information about the intensity of
                           its effect on gene expression. The highest coefficient(s) thus are the one(s) 
                           driving the profiles in a specific cluster. The genes in this cluster are potentially
                           involved in the response to that factor.
                           
                           The sign of a coefficient gives information about the way it imapcts expression.
                           If it is positive it increases the expression when the facot is in its perturbation
                           level. If negative, it decreases it.")
                      )
      )
    )
    
  )
}

# TODO : faire choisir le clustering préfait (comparaison) plutot que le dernier deja fait

#' cluster_exploration Server Function
#'
#' @noRd
mod_cluster_exploration_server <-
  function(input, output, session, r) {
    ns <- session$ns
    
    
    #   ____________________________________________________________________________
    #   Reactive memberships and conds                                          ####
    
    
    membership <- shiny::reactive({
      shiny::req(r$clusterings, r$current_comparison)
      req(r$clusterings[[r$current_comparison]])
      r$clusterings[[r$current_comparison]]$membership
    })
    
    
    conditions <- shiny::reactive({
      shiny::req(r$clusterings, r$current_comparison)
      req(r$clusterings[[r$current_comparison]])
      r$clusterings[[r$current_comparison]]$conditions
    })
    
    r_clust <- shiny::reactiveValues(go = NULL)
    
    #   ____________________________________________________________________________
    #   cluster to explore choice                                               ####
    
    
    output$cluster_to_explore_choice <- shiny::renderUI({
      shiny::req(membership())
      shinyWidgets::radioGroupButtons(
        inputId = ns("cluster_to_explore"),
        label = "Cluster to explore",
        choices = unique(membership()),
        justified = TRUE,
        checkIcon = list(yes = shiny::icon("ok",
                                           lib = "glyphicon"))
      )
    })
    
    
    #   ____________________________________________________________________________
    #   profiles                                                                ####
    
    
    output$profiles_to_explore <- shiny::renderPlot({
      # to reset the go analysis if new cluster
      r_clust$go <- NULL
      shiny::req(input$cluster_to_explore)
      draw_profiles(
        data = r$normalized_counts,
        membership = membership(),
        k = input$cluster_to_explore,
        conds = conditions()
      )
     
    })
    
    
    #   ____________________________________________________________________________
    #   table                                                                   ####
    
    
    output$genes_to_explore <- DT::renderDataTable({
      req(r$top_tags, r$current_comparison, membership())
      
      
      table <- r$top_tags[[r$current_comparison]]
      
      columns <- c("logFC", "logCPM", "FDR")
      if (!is.null(r$gene_info)) {
        
        if (r$splicing_aware) ids <- get_locus(rownames(table), unique = FALSE)
        else ids <- rownames(table)
        
        columns <- c(colnames(r$gene_info), columns)
        table[,colnames(r$gene_info)] <- r$gene_info[match(ids, rownames(r$gene_info)),]
      }
      
      table[table$genes %in% get_genes_in_cluster(membership = membership(),
                                                  cluster = input$cluster_to_explore),
            columns]
    })
    
    
    output$gene_number_cluster <- renderValueBox({
      valueBox(
        value = length(
          get_genes_in_cluster(
            membership = membership(),
            cluster = input$cluster_to_explore
          )
        ),
        subtitle = "genes in this cluster",
        color = "olive"
      )
    })
    
    
    
#   ____________________________________________________________________________
#   glm                                                                     ####

    
    glm <- shiny::reactive({
      shiny::req(r$design, r$normalized_counts, membership())
      fit_glm(normalized_counts = r$normalized_counts,
              genes = get_genes_in_cluster(membership = membership(),
                                           cluster = input$cluster_to_explore),
              design = r$design)
    })
    
    output$glm_summary <- shiny::renderPrint({
      print(summary(glm()))
    })
    
    output$glm_plot <- shiny::renderPlot({
      draw_glm(glm())
    })
    
    
    
    #   ____________________________________________________________________________
    #   GO enrich                                                               ####
    
    
    shiny::observeEvent((input$go_enrich_btn), {
      shiny::req(r$normalized_counts)
      shiny::req(membership())

      
      
      # for now, other orgs will come hopefully
      shiny::req(r$organism != "Other")
      
      
      genes <- get_genes_in_cluster(membership = membership(),
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
      
      r_clust$go <- enrich_go(genes, background, org = org)
    })
    
    #   ____________________________________________________________________________
    #   go results                                                              ####
    
    output$go_table <- DT::renderDataTable({
      shiny::req(r_clust$go)
      r_clust$go[,c("Description", "GeneRatio", "BgRatio", "p.adjust")]
    })
    
    output$go_plot <- plotly::renderPlotly({
      shiny::req(r_clust$go)
      draw_enrich_go(r_clust$go)
    })
    
    output$go_results <- shiny::renderUI({
      
      if(r$organism == "Other")
        shiny::h4("GO analysis is only supported for Arabidopsis and human (for now!)")
      
      shiny::req(r$organism != "Other")
      
      shiny::req(r_clust$go)
      if (!input$draw_go){
        DT::dataTableOutput(ns("go_table"))
      }
      else{
        plotly::plotlyOutput(ns("go_plot"), height = "700px")
      }
    })
  }

## To be copied in the UI
# mod_cluster_exploration_ui("cluster_exploration_ui_1")

## To be copied in the server
# callModule(mod_cluster_exploration_server, "cluster_exploration_ui_1")
