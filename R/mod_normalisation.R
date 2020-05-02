source("R/fct_normalisation.R")
source("R/fct_heatmap.R")
library(shinybusy)

#' normalisation UI Function
#'
#' @description A shiny Module for data filtering and normalisation.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @importFrom shinyWidgets actionBttn dropdownButton switchInput
#' @import shinydashboard
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_normalisation_ui <- function(id) {
  ns <- NS(id)
  tagList(
    
    shinybusy::add_busy_spinner(
      spin = "self-building-square",
      position = 'top-left',
      margins = c(600, 800)
    ),
    
    
    shiny::h1("Data filtering and normalisation"),
    shiny::hr(),
    shiny::h2(
      "Because low count genes and differences in sequencing depths are a true source of bias, we want to perform some data cleaning and transformation."
    ),
    shiny::hr(),
    col_3(
      boxPlus(
        title = "Parameters",
        solidHeader = F,
        status = "primary",
        collapsible = T,
        closable = F,
        width = 12,
        
        
        h2("Normalisation"),
        
        
        shinyWidgets::dropdownButton(
          size = 'xs',
          shiny::includeMarkdown(system.file("extdata", "normalisation.md", package = "DIANE")),
          circle = TRUE,
          status = "primary",
          icon = icon("question"),
          width = "600px",
          tooltip = shinyWidgets::tooltipOptions(title = "More details")
        ),
        
        h4("Prior removal of differentially expressed genes:"),

        shiny::fluidRow(col_12(shinyWidgets::switchInput(
          inputId = ns("prior_removal"),
          value = TRUE,
          onLabel = "ON",
          offLabel = "OFF",
          inline = T
        ))),
        
        
        
        shiny::fluidRow(
          
        col_8(shinyWidgets::awesomeRadio(
          inputId = "norm_method", label = "Normalisation :",
          choices = c("tmm", "deseq"),inline = T,
          selected = "tmm"
        )),
        
        col_4(shinyWidgets::actionBttn(
          ns("normalize_btn"),
          label = "Normalize",
          color = "primary",
          style = 'bordered'
        ))
        ),
        
        shiny::hr(),
        shiny::uiOutput(ns("norm_summary")),
        shiny::hr(),
     
        h2("Low counts filtering"),
        
        shinyWidgets::dropdownButton(
          size = 'xs',
          shiny::includeMarkdown(system.file("extdata", "normalisation.md", package = "DIANE")),
          circle = TRUE,
          status = "primary",
          icon = icon("question"),
          width = "600px",
          tooltip = shinyWidgets::tooltipOptions(title = "More details")
        ),

        
        shiny::h5("Minimal gene count sum accross conditions : "),
        col_8(shiny::numericInput(
          ns("low_counts_filter"),
          min = 0,
          value = 40,
          label = NULL
        )),
        col_4(
          shinyWidgets::actionBttn(
            ns("use_SumFilter"),
            label = "Filter",
            color = "primary",
            style = 'bordered'
          )
        ),
        
        
      shiny::hr(),
        
      # shinyWidgets::actionBttn(
      #   ns("use_HTSFilter"),
      #   label = "Use HST Filter, an entropy based filtering method",
      #   style = 'bordered',
      #   color = 'success'
      # ),
      shiny::hr(),
      shiny::hr(),
      shiny::uiOutput(ns("filtering_summary")),
      shiny::hr()
      
      )
    ),
    
    
    col_8(
      
      shinydashboard::tabBox(title = "Data exploration",
                             width = 12, height = "750px",
                             
        shiny::tabPanel(
          
          title = "Samples distributions",
          col_4(
            shinyWidgets::switchInput(
              inputId = ns("preview_norm_heatmap"),
              value = TRUE,
              onLabel = "normalized",
              offLabel = "raw",
              onStatus = "success",
              offStatus = "danger"
            )
          ),
          col_4(
            shinyWidgets::switchInput(
              inputId = ns("violin_preview"),
              value = TRUE,
              onLabel = "Boxplots",
              offLabel = "Violin"
            )
          ),
          shiny::plotOutput(ns('heatmap_preview_norm'), height = "600px")
        ),
        
        shiny::tabPanel(title = "MDS plot",
                        shiny::plotOutput(ns('mds_plot'), height = "600px"))
      )

    ),
    shiny::br()
    #DT::dataTableOutput(ns("norm_factor"))
  )
}

#' normalisation Server Function
#'
#' @noRd
mod_normalisation_server <- function(input, output, session, r) {
  ns <- session$ns
  
  shiny::observeEvent((input$normalize_btn), {
    req(r$raw_counts)
    norm <- normalize(r$raw_counts, r$conditions, norm_method = input$norm_method,
                      iteration = input$prior_removal)
    r$normalized_counts_pre_filter <- norm$normalized.counts
    r$norm_factors <- norm$norm_factors
    
  })
  
  
  output$norm_factor <- DT::renderDataTable({
    req(r$norm_factors)
    data.frame(t(round(r$norm_factors, 3)))
  })
  
  # shiny::observeEvent((input$use_HTSFilter), {
  #   shiny::req(r$normalized_counts_pre_filter)
  #   r$normalized_counts <-
  #     filter_hts(r$normalized_counts_pre_filter, conditions = r$conditions)
  # })
  
  
  shiny::observeEvent((input$use_SumFilter), {
    shiny::req(r$normalized_counts_pre_filter)
    r$normalized_counts <-
      filter_sum(r$normalized_counts_pre_filter, thr = input$low_counts_filter)
  })
  
  
  output$norm_summary <- shiny::renderUI({
    if (is.null(r$normalized_counts_pre_filter)) {
      number_color = "orange"
      number = "Normalisation needed"
      header = ""
      number_icon = "fa fa-times"
    }
    else{
      number_color = "olive"
      number = "Done"
      number_icon = "fa fa-check"
      header =  paste(dim(r$normalized_counts_pre_filter)[1],
                      " genes before filtering")
    }
    shinydashboardPlus::descriptionBlock(
      number = number, 
      number_color = number_color, 
      number_icon = number_icon,
      header = header,
      right_border = FALSE
    )
  })
  
  output$filtering_summary <- shiny::renderUI({
    if (is.null(r$normalized_counts_pre_filter)) {
      number_color = "red"
      number = "Normalisation needed"
      header = ""
      number_icon = "fa fa-times"
    }
    else{
      if (is.null(r$normalized_counts)) {
        number_color = "orange"
        number = "Filtering needed"
        header = ""
        number_icon = "fa fa-times"
      }
      else{
        number_color = "olive"
        number = "Done"
        number_icon = "fa fa-check"
        header = paste(dim(r$normalized_counts)[1],
                       " genes after filtering")
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
  
  
  output$heatmap_preview_norm <- shiny::renderPlot({
    shiny::req(r$normalized_counts_pre_filter)
    
    if (!input$preview_norm_heatmap) {
      # raw data
      d <- r$raw_counts 
    }
    else{
      if (is.null(r$normalized_counts)) {
        d <- r$normalized_counts_pre_filter
      }
      else{
        d <- r$normalized_counts
      }
    }
    draw_distributions(d, boxplot = input$violin_preview)
  })
  
  output$mds_plot <- shiny::renderPlot({
    shiny::req(r$normalized_counts)
    draw_MDS(r$normalized_counts)
  })
  
}

## To be copied in the UI
# mod_normalisation_ui("normalisation_ui_1")

## To be copied in the server
# callModule(mod_normalisation_server, "normalisation_ui_1")
