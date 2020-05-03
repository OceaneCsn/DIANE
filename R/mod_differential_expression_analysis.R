#' differential_expression_analysis UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_differential_expression_analysis_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::h1("Differential expression analysis - Now we're talking!"),
    shiny::hr(),
    
    shinybusy::add_busy_spinner(
      spin = "self-building-square",
      position = 'top-left',
      margins = c(800, 200)
    ),
#   ____________________________________________________________________________
#   Dispersion estimation                                                   ####

    col_3(
      boxPlus(
        title = "Parameters",
        solidHeader = F,
        status = "primary",
        collapsible = T,
        closable = F,
        width = 12,
        
        shiny::h4("Estimation of disperion : "),

        shinyWidgets::dropdownButton(
          size = 'xs',
          shiny::includeMarkdown(system.file("extdata", "normalisation.md", package = "DIANE")),
          circle = TRUE,
          status = "primary",
          icon = shiny::icon("question"),
          width = "600px",
          tooltip = shinyWidgets::tooltipOptions(title = "More details")
        ),
        
        shinyWidgets::actionBttn(
          ns("estimate_disp_btn"),
          label = "Launch estimation",
          color = "primary",
          style = 'bordered'
        ),
        
        
        shiny::hr(),
        shiny::uiOutput(ns("disp_estimate_summary")),
        shiny::hr(),
        
#   ____________________________________________________________________________
#   DEG parameters                                                          ####

        
        shiny::h4("Choose the conditions to compare for differential analysis : "),
        
        shinyWidgets::dropdownButton(
          size = 'xs',
          shiny::includeMarkdown(system.file("extdata", "normalisation.md", package = "DIANE")),
          circle = TRUE,
          status = "primary",
          icon = icon("question"),
          width = "600px",
          tooltip = shinyWidgets::tooltipOptions(title = "More details")
        ),
        
        
        shiny::uiOutput(ns("condition_choices")),

        
        col_8(shiny::numericInput(
          ns("dea_fdr"),
          min = 0,
          max = 1,
          value = 0.01,
          label = "Adjusted pvalue (fdr)"
        )),
        
          shinyWidgets::actionBttn(
            ns("deg_test_btn"),
            label = "Detect differentially expressed genes",
            color = "primary",
            style = 'bordered'
          ),
        
        shiny::hr(),
        shiny::uiOutput(ns("deg_test_summary")),
        shiny::hr()
        
      )
    ),

#   ____________________________________________________________________________
#   Visualisation of the results                                            ####

    col_8(
      shinydashboard::tabBox(title = "Results", width = 12,
                             shiny::tabPanel(title = "EdgeR summary",
                      shiny::verbatimTextOutput(ns("edgeR_summary")))
      )
      
    )
  )
}



#   __________________________________________________________________________________________________________________________________
#   Server                                                                                                                        ####

    
#' differential_expression_analysis Server Function
#'
#' @noRd 
mod_differential_expression_analysis_server <- function(input, output, session, r){
  ns <- session$ns
  
  r_dea <- shiny::reactiveValues(
    fit = NULL,
    top_tags = NULL
  )
  
  
#   ____________________________________________________________________________
#   Condition choices ui                                                    ####

  
  
  output$condition_choices <- shiny::renderUI({
    tagList(
      shinyWidgets::radioGroupButtons(
      inputId = ns("reference"),
      label = "Reference",
      choices = unique(r$conditions),
      justified = TRUE,
      checkIcon = list(
        yes = shiny::icon("ok", 
                          lib = "glyphicon"))
    ),
    
    shinyWidgets::radioGroupButtons(
      inputId = ns("perturbation"),
      label = "Perturbation",
      choices = unique(r$conditions),
      justified = TRUE,
      checkIcon = list(
        yes = shiny::icon("ok", 
                          lib = "glyphicon"))
    )
    )
  })
  
  shiny::observe({
    shiny::req(r$conditions)
    shinyWidgets::updateRadioGroupButtons(session, ns("reference"),
                            choices = ) 
    shinyWidgets::updateRadioGroupButtons(session, ns("perturbation"),
                            choices = r$conditions)
  })
  
#   ____________________________________________________________________________
#   Buttons reactives                                                       ####

  
  shiny::observeEvent((input$estimate_disp_btn), {
    shiny::req(r$tcc)
    r_dea$fit <- estimateDispersion(r$tcc)
    # the filtering needs to be done again if previously made, so :
    #r$normalized_counts <- NULL
    
  })
  
  shiny::observeEvent((input$deg_test_btn), {
    shiny::req(r_dea$fit, input$dea_fdr, input$reference, input$perturbation)
    shiny::req(input$reference != input$perturbation)
    print("before tests")
    r_dea$top_tags <- estimateDEGs(r_dea$fit, reference = input$reference, 
                                 perturbation = input$perturbation, 
                                 fdr = input$dea_fdr)
    
  })
  
  output$edgeR_summary <- shiny::renderPrint({
    shiny::req(r_dea$fit)
    print(r_dea$fit)
  })
  
  
#   ____________________________________________________________________________
#   Summaries                                                               ####

  output$disp_estimate_summary <- shiny::renderUI({
    if (is.null(r_dea$fit)) {
      number_color = "orange"
      number = "Needed"
      header = ""
      number_icon = "fa fa-times"
    }
    else{
      number_color = "olive"
      number = "Done"
      number_icon = "fa fa-check"
      header = "See EdgeR summary tab for more details"
    }
    shinydashboardPlus::descriptionBlock(
      number = number, 
      number_color = number_color, 
      number_icon = number_icon,
      header = header,
      right_border = FALSE
    )
  })
  
  output$deg_test_summary <- shiny::renderUI({
    if (is.null(r_dea$fit)) {
      number_color = "red"
      number = "Dispersion estimation needed"
      header = ""
      number_icon = "fa fa-times"
    }
    else{
      if (is.null(r_dea$top_tags)) {
        number_color = "orange"
        number = "Test can be performed"
        header = ""
        number_icon = "fa fa-times"
      }
      else{
        number_color = "olive"
        number = "Done"
        number_icon = "fa fa-check"
        header = "See plots and tables for more details"
      }
    }
    shinydashboardPlus::descriptionBlock(
      number = number, 
      number_color = number_color, 
      number_icon = number_icon,
      header = header,
      right_border = FALSE
    )
  })
}
    
## To be copied in the UI
# mod_differential_expression_analysis_ui("differential_expression_analysis_ui_1")
    
## To be copied in the server
# callModule(mod_differential_expression_analysis_server, "differential_expression_analysis_ui_1")
 
