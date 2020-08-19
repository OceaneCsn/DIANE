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
      shinydashboardPlus::boxPlus(
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
          col_6(
            shinyWidgets::numericInputIcon(ns("min_k"),
                             label = "Min number of clusters :",
                             value = 6, min = 0, max = 20,
                             help_text = "Minimum cluster number to test"
                          )
          ),
          col_6(
            shinyWidgets::numericInputIcon(ns("max_k"),
                                           label = "Max number of clusters :",
                                           value = 9, min = 0, max = 20,
                                           help_text = "Maximum cluster number to test"
            )
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
        shiny::uiOutput(ns("coseq_summary")),
        
        
        shiny::uiOutput(ns("dl_bttns"))
        
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
    ),
    shiny::actionButton(ns("browser"), "backdoor")
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
      numberColor = "orange"
      number = "Needed"
      header = ""
      numberIcon = "fa fa-times"
    }
    else{
      numberColor = "olive"
      number = "Done"
      numberIcon = "fa fa-check"
      header = "See Coseq summary tab for more details"
    }
    shinydashboardPlus::descriptionBlock(
      number = number,
      numberColor = numberColor,
      numberIcon = numberIcon,
      header = header,
      rightBorder = FALSE
    )
  })
  
  
  #   ____________________________________________________________________________
  #   Bttn reactive                                                           ####
  
  shiny::observeEvent((input$launch_coseq_btn), {
    shiny::req(r$normalized_counts, r$conditions)
    
    print(input$input_deg_genes)
    
    if(is.null(input$input_deg_genes)){
      shinyalert::shinyalert(
        "Please select an input gene list in the menu above",
        type = "error"
      )
    }
    
    shiny::req(input$input_deg_genes)

    genes_conditions <- unique(as.vector(
      stringr::str_split_fixed(input$input_deg_genes, ' ',2)))
    
    if (sum(genes_conditions %in% input$input_conditions) < length(genes_conditions)) {
      shinyalert::shinyalert(
        paste0( "Please select at least the conditions ",
          
          paste0(genes_conditions, collapse = ', ')
        ),
        "The conditions for clustering should contain the conditions
          used to compute the input differentially expressed genes.",
        type = "error"
      )
    }

    shiny::req(sum(genes_conditions %in% input$input_conditions) == length(genes_conditions))
    # union of all the input comparisons
    genes <- unique(unlist(r$DEGs[input$input_deg_genes]))
    
    run <-
      run_coseq(
        data = r$normalized_counts,
        genes = genes,
        conds = input$input_conditions,
        K = seq(input$min_k, input$max_k)
      )
    r$clusterings[[input_genes_conditions()]]$model <- run$model
    r$clusterings[[input_genes_conditions()]]$membership <-
      run$membership
    r$clusterings[[input_genes_conditions()]]$conditions <-
      input$input_conditions
    r$clusterings[[input_genes_conditions()]]$genes <- 
      input_genes_conditions()
    r$current_comparison <- input_genes_conditions()
    
  })
  
  
#   ____________________________________________________________________________
#   dl button                                                               ####

  
  output$dl_bttns <- shiny::renderUI({
    shiny::req(r$DEGs)
    shiny::req(r$clusterings)
    shiny::req(r$clusterings[[input_genes_conditions()]]$model)
    tagList(
      shiny::hr(),
      shiny::downloadButton(
          ns("report"), "Generate html report")
      )
  })

  
  
  #   ____________________________________________________________________________
  #   report                                                                  ####
  
  output$report <- shiny::downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "clustering_report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "clustering_report.Rmd")
      tempImage <- file.path(tempdir(), "favicon.ico")
      file.copy(system.file("extdata", "clustering_report.Rmd", package = "DIANE"), 
                tempReport, overwrite = TRUE)
      file.copy(system.file("extdata", "favicon.ico", package = "DIANE"),
                tempImage, overwrite = TRUE)      
      # Set up parameters to pass to Rmd document
      params <- list(r = r$clusterings[[input_genes_conditions()]],
                     input = input, normalized_counts = r$normalized_counts)
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  
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
  
  shiny::observeEvent(input$browser, {
    browser()
  })
  
}



## To be copied in the UI
# mod_clustering_ui("clustering_ui_1")

## To be copied in the server
# callModule(mod_clustering_server, "clustering_ui_1")
