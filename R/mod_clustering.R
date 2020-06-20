#' clustering UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_clustering_ui <- function(id) {
  ns <- NS(id)
  tagList(
    shiny::h1("Expression based clustering"),
    shiny::hr(),
    shinyalert::useShinyalert(),
    shinybusy::add_busy_spinner(
      spin = "self-building-square",
      position = 'top-left',
      margins = c(70, 1200)
    ),
    #   ____________________________________________________________________________
    #   Clustering settings                                                     ####
    
    col_4(
      boxPlus(
        title = "Settings",
        solidHeader = FALSE,
        status = "success",
        collapsible = TRUE,
        closable = FALSE,
        width = 12,
        
        col_10(shiny::h4("Poisson mixture model clustering")),
        
        col_2(
          shinyWidgets::dropdownButton(
            size = 'xs',
            shiny::includeMarkdown(
              system.file("extdata", "coseq.md", package = "DIANE")
            ),
            circle = TRUE,
            status = "success",
            icon = shiny::icon("question"),
            width = "600px",
            tooltip = shinyWidgets::tooltipOptions(title = "More details")
          )
        ),
        
        
        shiny::br(),
        shiny::hr(),
        col_12(shiny::uiOutput(
          ns("input_genes")
        )),
        
        col_12(shiny::uiOutput(
          ns("input_conditions_choice")
        )),
        
        shiny::hr(),
        shiny::fluidRow(col_12(
          shinyWidgets::sliderTextInput(
            inputId = ns("cluster_range"),
            label = "Min mand max number of clusters to test :",
            choices = seq(3, 15),
            selected = c(6, 9),
            from_min = 3,
            from_max = 7,
            to_min = 8,
            to_max = 15
          )
        )),
        shiny::br(),
        
        shiny::fluidRow(col_12(
          shinyWidgets::actionBttn(
            ns("launch_coseq_btn"),
            label = "Launch clustering",
            color = "success",
            style = 'bordered'
          )
        )),
        
        
        shiny::hr(),
        shiny::uiOutput(ns("coseq_summary"))
        
      )
    ),
    
    #   ____________________________________________________________________________
    #   Visualisation of the results                                            ####
    
    col_8(
      shinydashboard::tabBox(
        title = "Results",
        width = 12,
        
        
        shiny::tabPanel(
          title = "Cluster profiles",
          shiny::uiOutput(ns("profiles_clusters_choice")),
          shiny::plotOutput(ns("clusters_profiles")
                            , height = "800px")
        ),
        shiny::tabPanel(
          title = "Clustering quality",
          shiny::fluidRow(col_6(shiny::plotOutput(
            ns("plot_coseq_icl")
          )),
          col_6(shiny::plotOutput(
            ns("plot_coseq_barplots")
          ))),
          
          shiny::includeMarkdown(system.file("extdata", "ICL.md", package = "DIANE"))
        ),
        
        shiny::tabPanel(
          title = "Coseq summary",
          shiny::verbatimTextOutput(ns("coseq_run_summary")),
          shiny::h4(
            "The last line informs you of the number of clusters chosen to maximise 
            the ICL, the clustering quality criteria."
          )
        )
      )
    )
  )
}

#' clustering Server Function
#'
#' @noRd
mod_clustering_server <- function(input, output, session, r) {
  ns <- session$ns
  
  #   ____________________________________________________________________________
  #   degs input select                                                       ####
  
  
  # output$input_genes <- shiny::renderUI({
  #   shiny::req(r$DEGs)
  #   shinyWidgets::pickerInput(
  #     inputId = ns('input_deg_genes'),
  #     label = "Differentially expressed genes :",
  #     choices = names(r$DEGs),
  #     choicesOpt = list(subtext = paste(lengths(r$DEGs), "genes"))
  #   )
  # })
  
  output$input_genes <- shiny::renderUI({
    shiny::req(r$DEGs)
    
    shinyWidgets::checkboxGroupButtons(
      inputId = ns('input_deg_genes'),
      label = "Input genes for clustering :",
      choiceValues = names(r$DEGs),
      justified = TRUE,
      checkIcon = list(yes = shiny::icon("ok",
                                         lib = "glyphicon")),
      direction = "vertical",
      choiceNames = paste(names(r$DEGs), paste(lengths(r$DEGs), "genes"))
    )
  })
  
  input_genes_conditions <- shiny::reactive({
    req(input$input_deg_genes)
    print(paste(input$input_deg_genes, collapse = ' + '))
    return(paste(input$input_deg_genes, collapse = ' + '))
  })
  
  
  #   ____________________________________________________________________________
  #   conditions selection                                                    ####
  
  
  output$input_conditions_choice <- shiny::renderUI({
    shiny::req(r$DEGs, r$conditions)
    
    shinyWidgets::checkboxGroupButtons(
      inputId = ns('input_conditions'),
      label = "Conditions to perform clustering on :",
      choices = unique(r$conditions),
      justified = TRUE,
      checkIcon = list(yes = shiny::icon("ok",
                                         lib = "glyphicon")),
      selected = unique(r$conditions),
      direction = "vertical"
    )
  })
  
  
  
  #   ____________________________________________________________________________
  #   Clusters choices ui                                                     ####
  
  
  
  output$profiles_clusters_choice <- shiny::renderUI({
    shiny::req(r$clusterings, input$input_deg_genes)
    shiny::req(r$clusterings[[input_genes_conditions()]]$membership)
    tagList(
      shinyWidgets::checkboxGroupButtons(
        inputId = ns("clusters"),
        label = NULL,
        choices = unique(r$clusterings[[input_genes_conditions()]]$membership),
        justified = TRUE,
        checkIcon = list(yes = shiny::icon("ok",
                                           lib = "glyphicon")),
        selected = unique(r$clusterings[[input_genes_conditions()]]$membership)
      ),
    )
  })
  
  
  #   ____________________________________________________________________________
  #   coseq summary ui                                                        ####
  
  
  output$coseq_summary <- shiny::renderUI({
    shiny::req(r$clusterings)
    shiny::req(r$clusterings[[input_genes_conditions()]])
    if (is.null(r$clusterings[[input_genes_conditions()]]$model)) {
      number_color = "orange"
      number = "Needed"
      header = ""
      number_icon = "fa fa-times"
    }
    else{
      number_color = "olive"
      number = "Done"
      number_icon = "fa fa-check"
      header = "See Coseq summary tab for more details"
    }
    shinydashboardPlus::descriptionBlock(
      number = number,
      number_color = number_color,
      number_icon = number_icon,
      header = header,
      right_border = FALSE
    )
  })
  
  
  #   ____________________________________________________________________________
  #   Bttn reactive                                                           ####
  
  shiny::observeEvent((input$launch_coseq_btn), {
    shiny::req(r$normalized_counts, r$conditions, input_genes_conditions())
    print(input_genes_conditions())
    
    genes_conditions <- unique(as.vector(str_split_fixed(input$input_deg_genes, ' ',2)))
    print("genes conditions :")
    print(genes_conditions)
    
    if (!genes_conditions %in% input$input_conditions) {
      shinyalert::shinyalert(
        paste0(
          "The conditions used for clustering should contain the conditions
          used to compute the input differentially expressed genes. In that case : ",
          paste(input$input_deg_genes, collapse = ', ')
        ),
        type = "error"
      )
    }

    shiny::req(genes_conditions %in% input$input_conditions)
    # union of all the input comparisons
    genes <- unique(unlist(r$DEGs[input$input_deg_genes]))
    
    print(head(genes))
    run <-
      run_coseq(
        data = r$normalized_counts,
        genes = genes,
        conds = input$input_conditions,
        K = seq(input$cluster_range[1], input$cluster_range[2])
      )
    r$clusterings[[input_genes_conditions()]]$model <- run$model
    r$clusterings[[input_genes_conditions()]]$membership <-
      run$membership
    r$clusterings[[input_genes_conditions()]]$conditions <-
      input$input_conditions
    r$current_comparison <- input_genes_conditions()
    
  })

  
  
  #   ____________________________________________________________________________
  #   results                                                                 ####
  output$coseq_run_summary <- shiny::renderPrint({
    shiny::req(r$DEGs)
    shiny::req(r$clusterings)
    shiny::req(r$clusterings[[input_genes_conditions()]]$model)
    print(r$clusterings[[input_genes_conditions()]]$model)
  })
  
  output$plot_coseq_icl <- shiny::renderPlot({
    shiny::req(r$DEGs)
    shiny::req(r$clusterings)
    shiny::req(r$clusterings[[input_genes_conditions()]]$model)
    draw_coseq_run(r$clusterings[[input_genes_conditions()]]$model, plot = "ICL")
  })
  
  output$plot_coseq_barplots <- shiny::renderPlot({
    shiny::req(r$DEGs)
    shiny::req(r$clusterings, input$input_deg_genes)
    shiny::req(r$clusterings[[input_genes_conditions()]]$model)
    draw_coseq_run(r$clusterings[[input_genes_conditions()]]$model, plot = "barplots")
  })
  
  output$clusters_profiles <- shiny::renderPlot({
    shiny::req(r$DEGs)
    shiny::req(r$clusterings, input$input_deg_genes)
    shiny::req(r$clusterings[[input_genes_conditions()]]$membership, r$DEGs)
    shiny::req(input$clusters)
    draw_profiles(
      data = r$normalized_counts,
      conds = r$clusterings[[input_genes_conditions()]]$conditions,
      membership = r$clusterings[[input_genes_conditions()]]$membership,
      k = input$clusters
    )
  })
  
}



## To be copied in the UI
# mod_clustering_ui("clustering_ui_1")

## To be copied in the server
# callModule(mod_clustering_server, "clustering_ui_1")
