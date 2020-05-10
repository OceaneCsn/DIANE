#' clustering UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_clustering_ui <- function(id){
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
    
    col_3(
      boxPlus(
        title = "Settings",
        solidHeader = F,
        status = "success",
        collapsible = T,
        closable = F,
        width = 12,
        
        col_10(shiny::h4("Poisson mixture model clustering")),
        
        col_2(shinyWidgets::dropdownButton(
          size = 'xs',
          shiny::includeMarkdown(system.file("extdata", "normalisation.md", package = "DIANE")),
          circle = TRUE,
          status = "success",
          icon = shiny::icon("question"),
          width = "600px",
          tooltip = shinyWidgets::tooltipOptions(title = "More details")
        )),
        
        
        shiny::br(),
        shiny::hr(),
        shiny::fluidRow(col_10(shiny::uiOutput(ns("input_genes")))),
        shiny::hr(),
        shiny::fluidRow(col_12(shiny::uiOutput(ns("input_conditions_choice")))),

        shiny::hr(),  
        shiny::fluidRow(col_12(shinyWidgets::sliderTextInput(
          inputId = ns("cluster_range"),
          label = "Min mand max number of clusters to test :", 
          choices = seq(3,15),
          selected = c(6, 9),
          from_min = 3, 
          from_max = 7,
          to_min = 8,
          to_max = 15
        ))),
        shiny::br(),
        
        shiny::fluidRow(
          col_12(shinyWidgets::actionBttn(
            ns("launch_coseq_btn"),
            label = "Launch clustering",
            color = "success",
            style = 'bordered'
          ))),
        
        
        shiny::hr(),
        shiny::uiOutput(ns("coseq_summary"))
        
      )
    ),
    
    #   ____________________________________________________________________________
    #   Visualisation of the results                                            ####
    
    col_8(
      shinydashboard::tabBox(title = "Results", width = 12,


                               shiny::tabPanel(title = "Cluster profiles", 
                                               shiny::uiOutput(ns("profiles_clusters_choice")),
                                               shiny::plotOutput(ns("clusters_profiles")
                                                                             ,height = "800px"
                                                                             )),
                               shiny::tabPanel(title = "Clustering quality", 
                                               shiny::fluidRow(
                                                 col_6(shiny::plotOutput(ns("plot_coseq_icl"))),
                                                  col_6(shiny::plotOutput(ns("plot_coseq_barplots")))),
                                               
                                               shiny::includeMarkdown(system.file("extdata", "ICL.md", package = "DIANE"))),
                             
                               shiny::tabPanel(title = "Coseq summary",
                                shiny::verbatimTextOutput(ns("coseq_run_summary")),
                                shiny::h4("The last line informs you of the number of clusters chosen to maximise the ICL, the clustering 
                                   quality criteria."))
      )
    )
  )
}
    
#' clustering Server Function
#'
#' @noRd 
mod_clustering_server <- function(input, output, session, r){
  ns <- session$ns
  
#   ____________________________________________________________________________
#   degs input select                                                       ####

  
  output$input_genes <- shiny::renderUI({
    shiny::req(r$DEGs)
    shinyWidgets::pickerInput(
      inputId = ns('input_deg_genes'),
      label = "Differentially expressed genes :", 
      choices = names(r$DEGs),
      choicesOpt = list(
        subtext = paste(lengths(r$DEGs), "genes"))
    )
  })
  
  
#   ____________________________________________________________________________
#   conditions selection                                                    ####

  
  output$input_conditions_choice <- shiny::renderUI({
    shiny::req(r$DEGs)
    shinyWidgets::checkboxGroupButtons(
      inputId = ns('input_conditions'),
      label = "Conditions to perform clustering on :", 
      choices = unique(r$conditions),
      justified = TRUE,
      checkIcon = list(
        yes = shiny::icon("ok", 
                          lib = "glyphicon")),
      selected = unique(r$conditions)
    )
  })
  
  

  #   ____________________________________________________________________________
  #   Clusters choices ui                                                     ####
  
  
  
  output$profiles_clusters_choice <- shiny::renderUI({
    shiny::req(r$clusterings[[input$input_deg_genes]]$membership)
    tagList(
      
      shinyWidgets::checkboxGroupButtons(
        inputId = ns("clusters"),
        label = NULL,
        choices = unique(r$clusterings[[input$input_deg_genes]]$membership),
        justified = TRUE,
        checkIcon = list(
          yes = shiny::icon("ok", 
                            lib = "glyphicon")),
        selected = unique(r$clusterings[[input$input_deg_genes]]$membership)
      ),
    )
  })
  

  #   ____________________________________________________________________________
  #   coseq summary ui                                                        ####
  
  
  output$coseq_summary <- shiny::renderUI({
    if (is.null(r$clusterings[[input$input_deg_genes]]$model)) {
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
    shiny::req(r$normalized_counts, input$input_deg_genes, r$conditions)
    
    genes_conditions = as.vector(str_split_fixed(input$input_deg_genes, ' ', 2))
    if (length(intersect(genes_conditions, input$input_conditions)) < 2) {
      shinyalert::shinyalert(paste0("The conditions used for clustering should contain the conditions used to compute the input 
                             differentially expressed genes. In that case : ", input$input_deg_genes), type = "error")
    }
    shiny::req(length(intersect(genes_conditions, input$input_conditions)) == 2)
    run <- run_coseq(data = r$normalized_counts, genes = r$DEGs[[input$input_deg_genes]],
                         conds = input$input_conditions, K = seq(input$cluster_range[1], input$cluster_range[2]))
    r$clusterings[[input$input_deg_genes]]$model <- run$model
    r$clusterings[[input$input_deg_genes]]$membership <- run$membership
    r$clusterings[[input$input_deg_genes]]$conditions <- input$input_conditions
    r$current_comparison <- input$input_deg_genes
    
  })
  
  
  
    
  #   ____________________________________________________________________________
  #   results                                                                 ####
  output$coseq_run_summary <- shiny::renderPrint({
    shiny::req(r$DEGs)
    shiny::req(r$clusterings[[input$input_deg_genes]]$model, r$DEGs)
    print(r$clusterings[[input$input_deg_genes]]$model)
  })
  
  output$plot_coseq_icl <- shiny::renderPlot({
    shiny::req(r$DEGs)
    shiny::req(r$clusterings[[input$input_deg_genes]]$model, r$DEGs)
    draw_coseq_run(r$clusterings[[input$input_deg_genes]]$model)
  })
  
  output$plot_coseq_barplots <- shiny::renderPlot({
    shiny::req(r$DEGs)
    shiny::req(r$clusterings[[input$input_deg_genes]]$model, r$DEGs)
    draw_coseq_run(r$clusterings[[input$input_deg_genes]]$model, plot = "barplot")
  })
  
  output$clusters_profiles <- shiny::renderPlot({
    shiny::req(r$DEGs)
    shiny::req(r$clusterings[[input$input_deg_genes]]$membership, r$DEGs)
    draw_profiles(data = r$normalized_counts, conds = r$clusterings[[input$input_deg_genes]]$conditions, 
                  membership = r$clusterings[[input$input_deg_genes]]$membership,
                  k = input$clusters)
  })

}


    
## To be copied in the UI
# mod_clustering_ui("clustering_ui_1")
    
## To be copied in the server
# callModule(mod_clustering_server, "clustering_ui_1")
 
