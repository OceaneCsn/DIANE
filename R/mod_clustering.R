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
        
        shiny::h4("Input list of genes : "),
        
        col_2(shinyWidgets::dropdownButton(
          size = 'xs',
          shiny::includeMarkdown(system.file("extdata", "normalisation.md", package = "DIANE")),
          circle = TRUE,
          status = "success",
          icon = shiny::icon("question"),
          width = "600px",
          tooltip = shinyWidgets::tooltipOptions(title = "More details")
        )),
        
        shiny::fluidRow(shiny::uiOutput(ns("input_genes"))),

          
        col_12(shinyWidgets::sliderTextInput(
          inputId = ns("cluster_range"),
          label = "Min mand max number of clusters to test :", 
          choices = seq(3,15),
          selected = c(6, 12),
          from_min = 3, 
          from_max = 7,
          to_min = 8,
          to_max = 15
        )),
        
        shiny::fluidRow(
          col_12(shinyWidgets::actionBttn(
            ns("launch_coseq_btn"),
            label = "Launch clustering",
            color = "success",
            style = 'bordered'
          )),
        
        shiny::fluidRow(
        shiny::br(),
        shiny::br(),
        shiny::hr(),
        shiny::uiOutput(ns("coseq_summary")),
        shiny::hr()
      )))
    ),
    
    #   ____________________________________________________________________________
    #   Visualisation of the results                                            ####
    
    col_8(
      shinydashboard::tabBox(title = "Results", width = 12,
                             # shiny::tabPanel(title = "Results table",
                             #                 DT::dataTableOutput(ns("deg_table"))),
                             # shiny::tabPanel(title = "MA - Vulcano plots",
                             #                 
                             #                 shinyWidgets::switchInput(
                             #                   inputId = ns("MA_vulcano_switch"),
                             #                   value = TRUE,
                             #                   onLabel = "MA",
                             #                   offLabel = "Vulcano"
                             #                 ),
                             #                 
                             #                 
                             #                 shiny::plotOutput(ns("ma_vulcano"), height = "700px")
                             #                 
                             # ),
                             shiny::tabPanel(title = "Clustering quality", col_6(shiny::plotOutput(ns("plot_coseq_icl"))),
                                             col_6(shiny::plotOutput(ns("plot_coseq_barplots")))),
                             shiny::tabPanel(title = "Coseq summary",
                              shiny::verbatimTextOutput(ns("coseq_run_summary")))
                        
      )
      
    )
  )
}
    
#' clustering Server Function
#'
#' @noRd 
mod_clustering_server <- function(input, output, session, r){
  ns <- session$ns
  
  
  r_coseq <- shiny::reactiveValues(model = NULL,
                                   membership=NULL)
  
#   ____________________________________________________________________________
#   degs input select                                                       ####

  
  output$input_genes <- shiny::renderUI({
    shiny::req(r$DEGs)
    shinyWidgets::pickerInput(
      inputId = ns('input_deg_genes'),
      label = "Differentially expressed genes", 
      choices = names(r$DEGs),
      choicesOpt = list(
        subtext = paste(lengths(r$DEGs), "genes"))
    )
  })
  
  
  
  #   ____________________________________________________________________________
  #   coseq summary ui                                                        ####
  
  
  output$coseq_summary <- shiny::renderUI({
    if (is.null(r_coseq$model)) {
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
    print(unique(r$conditions))
    print(input$cluster_range)
    print(r$DEGs[[input$input_deg_genes]])
    run <- run_coseq(data = r$normalized_counts, genes = r$DEGs[[input$input_deg_genes]],
                         conds = unique(r$conditions), K = seq(input$cluster_range[1], input$cluster_range[2]))
    r_coseq$model <- run$model
    r_coseq$membership <- run$membership
    print("Done")
  })
  
  
  
    
  #   ____________________________________________________________________________
  #   results                                                                 ####
  output$coseq_run_summary <- shiny::renderPrint({
    shiny::req(r_coseq$model)
    print(r_coseq$model)
  })
  
  output$plot_coseq_icl <- shiny::renderPlot({
    shiny::req(r_coseq$model)
    draw_coseq_run(r_coseq$model)
  })
  
  output$plot_coseq_barplots <- shiny::renderPlot({
    shiny::req(r_coseq$model)
    draw_coseq_run(r_coseq$model, plot = "barplot")
  })
    
}


    
## To be copied in the UI
# mod_clustering_ui("clustering_ui_1")
    
## To be copied in the server
# callModule(mod_clustering_server, "clustering_ui_1")
 
