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
    shiny::fluidRow(
      col_8(shiny::uiOutput(ns(
        "cluster_to_explore_choice"
      ))),
      col_4(shinydashboard::valueBoxOutput(ns(
        "gene_number_cluster"
      )))
    ),
    shiny::hr(),
    
    #   ____________________________________________________________________________
    #   Profiles column                                                         ####
    
    col_6(
      shinydashboardPlus::boxPlus(
        width = 12,
        closable = FALSE,
        title = "Expression profiles",
        shiny::plotOutput(ns("profiles_to_explore"), height = "700px"),
        shiny::br(),
        shiny::fluidRow(col_12(
          shinyWidgets::downloadBttn(
            outputId = ns("download_genes_in_cluster"),
            label = "Download genes in this cluster as a csv table",
            style = "material-flat",
            color = "success"
          )
        ))
      )
    ),
    
    # TODO true reset when toggle demo data
    #   ____________________________________________________________________________
    #   Clusters characteristics                                                ####
    
    
    col_6(
      shinydashboard::tabBox(
        title = "Genes in that clusters",
        width = 12,
        shiny::tabPanel(title = "Genes table",
                        DT::dataTableOutput(ns(
                          "genes_to_explore"
                        ))),
        
        
        #   ____________________________________________________________________________
        #   GO                                                                      ####
        
        shiny::tabPanel(
          title = "Gene Ontologies enrichment",
          
          
          col_4(
            shinyWidgets::actionBttn(
              ns("go_enrich_btn"),
              label = "Start GO enrichment analysis",
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
          shiny::uiOutput(ns("custom_data_go")),
          
          shiny::hr(),
          
          shiny::fluidRow(col_12(shiny::uiOutput(ns(
            "go_results"
          ))))
        ),
        
        #   ____________________________________________________________________________
        #   glm                                                                     ####
        
        shiny::tabPanel(
          title = "GLM for factors effect",
          col_2(
            shinyWidgets::dropdownButton(
              size = 'xs',
              shiny::verbatimTextOutput(ns("glm_summary")),
              circle = FALSE,
              status = "success",
              width = "600px",
              label = "Glm summary"
            )
          ),
          
          col_2(
            shinyWidgets::dropdownButton(
              size = 'xs',
              shiny::includeMarkdown(system.file("extdata", "pglm.md", package = "DIANE")),
              circle = TRUE,
              status = "success",
              icon = shiny::icon("question"),
              width = "600px",
              tooltip = shinyWidgets::tooltipOptions(title = "More details")
            )
          ),
          
          shiny::plotOutput(ns("glm_plot"), height = "700px"),
          shiny::hr(),
          shiny::h5(
            "The absolute value of a coefficient gives information about the intensity of
                           its effect on gene expression. The highest coefficient(s) thus are the one(s)
                           driving the profiles in a specific cluster. The genes in this cluster are potentially
                           involved in the response to that factor.

                           The sign of a coefficient gives information about the way it impacts expression.
                           If it is positive it increases the expression when the facot is in its perturbation
                           level. If negative, it decreases it."
          )
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
    
    table <- shiny::reactive({
      req(r$top_tags, r$current_comparison, membership())
      
      genes <- get_genes_in_cluster(membership = membership(),
                                    cluster = input$cluster_to_explore)
      table <- data.frame(Genes = genes)
      
      if (!is.null(r$gene_info)) {
        if (r$splicing_aware)
          ids <- get_locus(genes, unique = FALSE)
        else
          ids <- genes
        table[, colnames(r$gene_info)] <-
          r$gene_info[match(ids, rownames(r$gene_info)), ]
      }
      else{
        table
      }
    })
    
    
    output$genes_to_explore <- DT::renderDataTable({
      table()
    })
    
    
    output$gene_number_cluster <- shinydashboard::renderValueBox({
      shinydashboard::valueBox(
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
    
    
    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
    ### download table                                                          ####
    
    
    output$download_genes_in_cluster <- shiny::downloadHandler(
      filename = function() {
        paste(paste0("genes_cluster_", input$cluster_to_explore, ".csv"))
      },
      content = function(file) {
        write.csv(table(), file = file, quote = FALSE)
      }
    )
    
    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
    ### download GO                                                             ####
    
    
    output$download_go_table <- shiny::downloadHandler(
      filename = function() {
        paste(paste0(
          "enriched_GOterms_cluster_",
          input$cluster_to_explore,
          ".csv"
        ))
      },
      content = function(file) {
        write.csv(r_clust$go, file = file, quote = FALSE)
      }
    )
    
    #   ____________________________________________________________________________
    #   glm                                                                     ####
    
    
    glm <- shiny::reactive({
      shiny::req(r$design, r$normalized_counts, membership())
      fit_glm(
        normalized_counts = r$normalized_counts,
        genes = get_genes_in_cluster(
          membership = membership(),
          cluster = input$cluster_to_explore
        ),
        design = r$design,
        factors = get_factors_from_conditions(conditions(), r$design)
      )
    })
    
    output$glm_summary <- shiny::renderPrint({
      print(summary(glm()))
    })
    
    output$glm_plot <- shiny::renderPlot({
      draw_glm(glm())
    })
    #   ____________________________________________________________________________
    #   custom go                                                               ####
    
    output$custom_data_go <- shiny::renderUI({
      shiny::req(r$organism == "Other")
      shiny::req(is.null(r$custom_go))
      
      tagList(
        col_2(
          shinyWidgets::dropdownButton(
            size = 'xs',
            shiny::includeMarkdown(system.file("extdata", "custom_go.md", package = "DIANE")),
            circle = TRUE,
            status = "success",
            icon = shiny::icon("question"),
            width = "600px",
            tooltip = shinyWidgets::tooltipOptions(title = "More details")
          )
        ),
        col_10(
          shiny::h4(
            "Your organism is not known to DIANE, but you can provide a matching between
         gene IDs and GO IDs."
          )
        ),
        
        
        col_6(shiny::radioButtons(
          ns('sep'),
          
          'Separator : ',
          c(
            Comma = ',',
            Semicolon = ';',
            Tab = '\t'
          ),
          inline = TRUE
        )),
        
        col_6(
          shiny::fileInput(
            ns('go_data'),
            'Choose CSV/TXT GO terms file',
            accept = c(
              'text/csv',
              'text/comma-separated-values,text/plain',
              '.csv',
              '.txt'
            )
          )
        )
      )
      
    })
    
    
    
    #   ____________________________________________________________________________
    #   GO enrich                                                               ####
    
    
    shiny::observeEvent((input$go_enrich_btn), {
      shiny::req(r$normalized_counts)
      shiny::req(membership())
      
      if (r$organism == "Other") {
        if (is.null(r$custom_go)) {
          if (!is.null(input$go_data)) {
            pathName = input$go_data$datapath
            d <- read.csv(
              sep = input$sep,
              file = pathName,
              header = TRUE,
              stringsAsFactors = FALSE
            )
            print(ncol(d))
            r$custom_go <- d
          }
          else{
            shinyalert::shinyalert(
              "Please input Gene to GO term file. ",
              "For now, only Arabidopsis thaliana and
        Homo sapiens are supported, but you can input your own gene - GO terms matching.",
              type = "error"
            )
          }
        }
        shiny::req(r$custom_go)
        if (ncol(r$custom_go) != 2) {
          r$custom_go <- NULL
          shinyalert::shinyalert(
            "Invalid file",
            "It must contain two columns as described.
            Did you correctly set the separator?",
            type = "error"
          )
        }
        
        shiny::req(ncol(r$custom_go) == 2)
        
        GOs <- r$custom_go
        genes <- get_genes_in_cluster(membership = membership(),
                                      cluster = input$cluster_to_explore)
        universe <-
          intersect(rownames(r$normalized_counts), GOs[, 1])
        
        r_clust$go <- enrich_go_custom(genes, universe, GOs)
        
      } else{
        genes <- get_genes_in_cluster(membership = membership(),
                                      cluster = input$cluster_to_explore)
        
        background <- rownames(r$normalized_counts)
        
        if (r$splicing_aware) {
          genes <- get_locus(genes)
          background <- get_locus(background)
        }
        
        if (r$organism == "Lupinus albus") {
          GOs <- DIANE:::lupine$go_list
          universe <- intersect(background, GOs[, 1])
          r_clust$go <- enrich_go_custom(genes, universe, GOs)
        }
        
        else{
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
          
          if (r$organism == "Drosophilia melanogaster") {
            genes <- convert_from_ensembl_dm(genes)
            background <- convert_from_ensembl_dm(background)
            org = org.Dm.eg.db::org.Dm.eg.db
          }
          
          if (r$organism == "Caenorhabditis elegans") {
            genes <- convert_from_ensembl_ce(genes)
            background <- convert_from_ensembl_ce(background)
            org = org.Ce.eg.db::org.Ce.eg.db
          }
          
          
          if (r$organism == "Escherichia coli") {
            genes <- convert_from_ensembl_eck12(genes)
            background <- convert_from_ensembl_eck12(background)
            org = org.EcK12.eg.db::org.EcK12.eg.db
          }
          
          # TODO add check if it is entrez with regular expression here
          shiny::req(length(genes) > 0, length(background) > 0)
          
          r_clust$go <-
            enrich_go(genes,
                      background,
                      org = org,
                      GO_type = input$go_type)
          
        }
        
      }
      
      if (golem::get_golem_options("server_version"))
        loggit::loggit(
          custom_log_lvl = TRUE,
          log_lvl = r$session_id,
          log_msg = "GO enrichment cluster"
        )
    })
    
    #   ____________________________________________________________________________
    #   go results                                                              ####
    
    output$go_table <- DT::renderDataTable({
      shiny::req(r_clust$go)
      r_clust$go[, c("Description", "GeneRatio", "BgRatio", "p.adjust")]
    })
    
    output$max_go_choice <- shiny::renderUI({
      shiny::req(r_clust$go)
      shiny::numericInput(
        ns("n_go_terms"),
        label = "Top number of GO terms to plot :",
        min = 1,
        value = dim(r_clust$go)[1]
      )
    })
    
    output$go_plot <- plotly::renderPlotly({
      shiny::req(r_clust$go)
      max = ifelse(is.na(input$n_go_terms),
                   dim(r_clust$go)[1],
                   input$n_go_terms)
      draw_enrich_go(r_clust$go, max_go = max)
    })
    
    output$go_map_plot <- shiny::renderPlot({
      shiny::req(r_clust$go)
      draw_enrich_go_map(r_clust$go)
    })
    
    output$go_results <- shiny::renderUI({
      shiny::req(r_clust$go)
      
      if (nrow(r_clust$go) == 0) {
        shinyalert::shinyalert(
          "No enriched GO terms were found",
          "It can happen if input gene list is not big enough",
          type = "error"
        )
      }
      shiny::req(nrow(r_clust$go) > 0)
      
      if (input$draw_go == "Data table") {
        tagList(
          DT::dataTableOutput(ns("go_table")),
          shinyWidgets::downloadBttn(
            outputId = ns("download_go_table"),
            label = "Download enriched GO term as a csv table",
            style = "material-flat",
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
# mod_cluster_exploration_ui("cluster_exploration_ui_1")

## To be copied in the server
# callModule(mod_cluster_exploration_server, "cluster_exploration_ui_1")
