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
    
    shiny::fluidRow(
      shinydashboardPlus::box(
        title = "Settings",
        solidHeader = FALSE,
        status = "success",
        collapsible = TRUE,
        closable = FALSE,
        width = 4,
        
        col_10(shiny::h4("Mixture models clustering")),
        
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
          ),
          
        shiny::br(),
        
        
        col_6(shinyWidgets::switchInput(
          inputId = ns("coseq_model"),
          label = "Mixture Model to use",
          value = TRUE,
          onLabel = "Poisson",
          offLabel = "Normal",
          onStatus = 'success',
          offStatus = 'success'
        )),
        
        col_6(shiny::uiOutput(ns("transfo_ui"))),
        
        shiny::fluidRow(shiny::column(12, align="center",
          shinyWidgets::actionBttn(
            ns("launch_coseq_btn"),
            label = "Launch clustering",
            style = "material-flat",
            color = "success"
          )
        )),
        
        
        shiny::hr(),
        shiny::uiOutput(ns("coseq_summary")),
        
        shiny::fluidRow(shiny::column(12, align="center",
                                      shiny::uiOutput(ns("dl_bttns"))
        ))
        
      ),
    
    #   ____________________________________________________________________________
    #   Visualisation of the results                                            ####
    
      shinydashboard::tabBox(
        title = "Results",
        width = 8,
        
        
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
    if(length(r$DEGs) > 0){
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
    }
    else{
      shinydashboardPlus::descriptionBlock(
        number = "Please perform one or more differential expression analysis before clustering",
        numberColor = "orange",
        rightBorder = FALSE
      )
    }
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
  #   Transfos choices ui                                                     ####
  
  
  
  output$transfo_ui <- shiny::renderUI({
    shiny::req(r$normalized_counts)
    shiny::req(!input$coseq_model)
    tagList(
      shiny::selectInput(
        inputId = ns("transfo"),
        label = "Transformation prior to Gaussian Mixtures:",
        choices = c("voom", "arcsin", "logit", 
                    "logMedianRef", "profile", 
                    "logclr", "clr", "alr", 
                    "ilr", "none"),
        selected = "arcsin"
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
      numberIcon = shiny::icon('times')
    }
    else{
      numberColor = "olive"
      number = "Done"
      numberIcon = shiny::icon('check')
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
  
  
  # rep <- shiny::reactiveValues(
  #   coseq = shiny::repeatable(DIANE::run_coseq, seed = 123)
  # )
  
  
  #   ____________________________________________________________________________
  
  #   Bttn reactive                                                           ####
  
  shiny::observeEvent((input$launch_coseq_btn), {
    shiny::req(r$normalized_counts, r$conditions)
    

    if(is.null(input$input_deg_genes)){
      shinyalert::shinyalert(
        "Please select an input gene list in the menu above",
        type = "error"
      )
    }
    
    shiny::req(input$input_deg_genes)

    genes_conditions <- unique(unlist(
      stringr::str_extract_all(
        unlist(stringr::str_extract_all(input$input_deg_genes, pattern = "[^()+ ]+")),
        pattern = "[^\\s\\+]+"
      )
    ))
    
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
    
    if(input$coseq_model){
      mod <- "Poisson"
      transfo <- "none"
    }
    else{
      mod <- "Normal"
      transfo <- input$transfo
    }
    
    run <-
      run_coseq(
        data = r$normalized_counts,
        genes = genes,
        conds = input$input_conditions,
        transfo = transfo,
        model = mod,
        K = seq(input$min_k, input$max_k),
        seed = r$seed
      )
    
    r$clusterings[[input_genes_conditions()]]$model <- run$model
    r$clusterings[[input_genes_conditions()]]$membership <-
      run$membership
    r$clusterings[[input_genes_conditions()]]$conditions <-
      input$input_conditions
    r$clusterings[[input_genes_conditions()]]$genes <- 
      input_genes_conditions()
    r$current_comparison <- input_genes_conditions()
    
    if(golem::get_golem_options("server_version"))
      loggit::loggit(custom_log_lvl = TRUE,
                   log_lvl = r$session_id,
                   log_msg = "clustering")
    
  })
  
  
#   ____________________________________________________________________________
#   dl button                                                               ####

  
  output$dl_bttns <- shiny::renderUI({
    shiny::req(r$DEGs)
    shiny::req(r$clusterings)
    shiny::req(r$clusterings[[input_genes_conditions()]]$model)
    tagList(
      shiny::hr(),
      
      shinyWidgets::downloadBttn(
        ns("report"), "Generate html report",
        style = "material-flat", color = "default")
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
      params <- list(r_clust = r$clusterings[[input_genes_conditions()]],
                     r = r,
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
  

}



## To be copied in the UI
# mod_clustering_ui("clustering_ui_1")

## To be copied in the server
# callModule(mod_clustering_server, "clustering_ui_1")
