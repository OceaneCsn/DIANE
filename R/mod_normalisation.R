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
    shiny::h1("Data filtering and normalisation"),
    shiny::h2(
      "Because low count genes and differences in sequencing depths are a true source of bias, we want to perform some data cleaning and transformation."
    ),
    shiny::hr(),
    col_6(
      boxPlus(
        title = "Sample-wise Normalisation",
        solidHeader = F,
        height = "600px",
        status = "primary",
        collapsible = T,
        closable = F,
        col_6(h3("TCC package : ")),
        
        
        col_4(
          shinyWidgets::actionBttn(
            ns("normalize_btn"),
            label = "Normalize",
            color = "primary",
            style = 'bordered'
          )
        ),
        col_2(
          shinyWidgets::dropdownButton(
            size = 'xs',
            tags$h3("Normalisation method :"),
            h5(
              "TCC performs normalisation to correct for different secencing depth between the samples, so they can be comparable in terms of transcript counts.
            Most normalisation methods rely on the hypothesis that very few genes are differentially expressed. TCC proposes
             a prior step to first detect differentially expressed genes, and remove them to perform the final normalisation. If you suspect you have a lot of differentially
             expressed genes in your data, keep \"remove differentially expressed genes first\" enabled. Else, standard normalisation will be performed, using either TMM or DESeq2
             "
            ),
            circle = TRUE,
            status = "info",
            icon = icon("question"),
            width = "300px",
            tooltip = tooltipOptions(title = "More details")
          )
        ),
        
        shiny::hr(),
        
        
        shinybusy::add_busy_spinner(
          spin = "self-building-square",
          position = 'top-left',
          margins = c(500, 500)
        ),
        
        valueBoxOutput(ns("norm_summary"), width = 12)
      ),
      
      boxPlus(
        title = "Low counts filtering",
        solidHeader = F,
        height = "600px",
        status = "primary",
        collapsible = T,
        closable = F,
        
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
        
        
        shiny::br(),
        
        shinyWidgets::actionBttn(
          ns("use_HTSFilter"),
          label = "Use HST Filter, an entropy based filtering method",
          style = 'bordered',
          color = 'success'
        ),
        shiny::hr(),
        valueBoxOutput(ns("filtering_summary"), width = 12)
        #footer = "This might help you visualize the different sequencing depths of your conditions."
      )
    ),
    col_6(
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
          inputId = ns("log_preview"),
          value = TRUE,
          onLabel = "Log(Value+1)",
          offLabel = "Value"
        )
      ),
      shiny::plotOutput(ns('heatmap_preview_norm'))
      
    ),
    shiny::br(),
    DT::dataTableOutput(ns("norm_factor"))
    
    
  )
}

#' normalisation Server Function
#'
#' @noRd
mod_normalisation_server <- function(input, output, session, r) {
  ns <- session$ns
  
  shiny::observeEvent((input$normalize_btn), {
    norm <- normalize(r$raw_counts, r$conditions)
    r$normalized_counts_pre_filter <- norm$normalized.counts
    r$norm_factors <- norm$norm_factors
    shinyWidgets::sendSweetAlert(session = session,
                                 title = "Normalisation completed, you can proceed to filtering now",
                                 type = "success")
    
  })
  
  
  output$norm_factor <- DT::renderDataTable({
    req(r$norm_factors)
    data.frame(t(round(r$norm_factors, 3)))
  })
  
  # 
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
  
  
  output$norm_summary <- renderValueBox({
    if (is.null(r$normalized_counts_pre_filter)) {
      color = "orange"
      value = "Normalisation needed"
      text = ""
    }
    else{
      color = "olive"
      text = paste(dim(r$normalized_counts_pre_filter)[1],
                   " genes before filtering")
      value = "Done"
    }
    valueBox(value,
             text,
             color = color)
  })
  
  output$filtering_summary <- renderValueBox({
    if (is.null(r$normalized_counts_pre_filter)) {
      color = "red"
      value = "Normalisation needed"
      text = ""
    }
    else{
      if (is.null(r$normalized_counts)) {
        color = "orange"
        text = ""
        value = "Filetring needed"
      }
      else{
        color = "olive"
        text = paste(dim(r$normalized_counts)[1], " genes after filtering")
        value = "Done"
      }
    }
    valueBox(value,
             text,
             color = color)
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
    print(typeof(d))
    print(d)
    draw_heatmap(d, title = "", log = input$log_preview)
  })
  
}

## To be copied in the UI
# mod_normalisation_ui("normalisation_ui_1")

## To be copied in the server
# callModule(mod_normalisation_server, "normalisation_ui_1")
