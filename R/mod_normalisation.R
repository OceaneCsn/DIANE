# FIXME silent melt in reshape2

# TODO markdown dea, coseq, glm


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
      margins = c(70, 1000)
    ),
    
    shiny::h1("Data filtering and normalisation"),
    shiny::hr(),
    # shiny::h2(
    #   "Because low count genes and differences in sequencing depths are a 
    #   source of bias, we want to perform some data cleaning and transformation."
    # ),
    # shiny::hr(),
    
    #   ____________________________________________________________________________
    #   Normalisation settings                                                  ####
    
    col_3(
      boxPlus(
        title = "Settings",
        solidHeader = FALSE,
        status = "success",
        collapsible = TRUE,
        closable = FALSE,
        width = 12,
        
        
        shiny::h2("Normalisation"),
        
        
        shinyWidgets::dropdownButton(
          size = 'xs',
          shiny::includeMarkdown(
            system.file("extdata", "normalisation.md", package = "DIANE")
          ),
          circle = TRUE,
          status = "success",
          icon = shiny::icon("question"),
          width = "600px",
          tooltip = shinyWidgets::tooltipOptions(title = "More details")
        ),
        
        shiny::h4("Prior removal of differentially expressed genes:"),
        
        shiny::fluidRow(col_12(
          shinyWidgets::switchInput(
            inputId = ns("prior_removal"),
            value = FALSE,
            onLabel = "ON",
            offLabel = "OFF",
            inline = TRUE,
            onStatus = "success"
          )
        )),
        
        
        
        shiny::fluidRow(
          col_8(
            shinyWidgets::awesomeRadio(
              inputId = ns("norm_method"),
              label = "Normalisation method:",
              choices = c("tmm", "deseq"),
              inline = TRUE,
              selected = "tmm",
              status = "success"
            )
          ),
          
          col_4(
            shinyWidgets::actionBttn(
              ns("normalize_btn"),
              label = "Normalize",
              color = "success",
              style = 'bordered'
            )
          )
        ),
        
        shiny::hr(),
        shiny::uiOutput(ns("norm_summary")),
        shiny::hr(),
        
        
        #   ____________________________________________________________________________
        #   filtering settings                                                      ####
        
        
        shiny::h2("Low counts filtering"),
        
        shinyWidgets::dropdownButton(
          size = 'xs',
          shiny::includeMarkdown(system.file("extdata", "filtering.md", package = "DIANE")),
          circle = TRUE,
          status = "success",
          icon = shiny::icon("question"),
          width = "600px",
          tooltip = shinyWidgets::tooltipOptions(title = "More details")
        ),
        
        
        shiny::h5("Minimal gene count sum accross conditions : "),
        col_8(shiny::uiOutput(ns(
          "filter_proposition"
        ))),
        col_4(
          shinyWidgets::actionBttn(
            ns("use_SumFilter"),
            label = "Filter",
            color = "success",
            style = 'bordered'
          )
        ),
        
        
        shiny::hr(),
        
        shiny::hr(),
        shiny::hr(),
        shiny::uiOutput(ns("filtering_summary")),
        shiny::hr(),
        
        shiny::br(),
        
        shiny::uiOutput(ns("dl_bttns"))
      )
    ),
    
    
    #   ____________________________________________________________________________
    #   plot results ui                                                         ####
    
    col_8(
      shinydashboard::tabBox(
        title = "Data exploration",
        width = 12,
        height = "750px",
        
        shiny::tabPanel(
          title = "Samples distributions",
          col_4(
            shinyWidgets::switchInput(
              inputId = ns("preview_norm_distr"),
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
              offLabel = "Distributions",
              onStatus = "success"
            )
          ),
          shiny::plotOutput(ns('heatmap_preview_norm'), height = "900px")
        ),
        shiny::tabPanel(title = "PCA",
                        shiny::plotOutput(ns('pca_plot'), height = "700px")),
        
        shiny::tabPanel(title = "MDS plot",
                        shiny::plotOutput(ns('mds_plot'), height = "600px")),
        
        shiny::tabPanel(title = "Summary",
                        shiny::verbatimTextOutput(ns("tcc_summary")))
      )
      
    ),
    shiny::br()
  )
}

#' normalisation Server Function
#' @importFrom TCC getNormalizedData
#' @importFrom utils write.csv
#' @noRd
mod_normalisation_server <- function(input, output, session, r) {
  ns <- session$ns
  
  
  output$filter_proposition <- shiny::renderUI({
    shiny::numericInput(
      ns("low_counts_filter"),
      min = 0,
      value = 10 * length(r$conditions),
      label = NULL
    )
  })
  
  
  #   ____________________________________________________________________________
  #   buttn reactives                                                         ####
  
  shiny::observeEvent(input$normalize_btn, {
    shiny::req(r$raw_counts)
    r$tcc <-
      normalize(
        r$raw_counts,
        r$conditions,
        norm_method = input$norm_method,
        iteration = input$prior_removal
      )
    r$normalized_counts_pre_filter <- TCC::getNormalizedData(r$tcc)
    # the filtering needs to be done again if previously made, so :
    r$normalized_counts <- NULL
    
  })
  
  shiny::observeEvent((input$use_SumFilter), {
    shiny::req(r$normalized_counts_pre_filter)
    r$tcc <- filter_low_counts(r$tcc, thr = input$low_counts_filter)
    r$normalized_counts <- TCC::getNormalizedData(r$tcc)
    
  })
  
  
  #   ____________________________________________________________________________
  #   summaries                                                               ####
  
  output$norm_summary <- shiny::renderUI({
    if (is.null(r$normalized_counts_pre_filter)) {
      numberColor = "orange"
      number = "Normalisation needed"
      header = ""
      numberIcon = "fa fa-times"
    }
    else{
      numberColor = "olive"
      number = "Done"
      numberIcon = "fa fa-check"
      header =  paste(dim(r$normalized_counts_pre_filter)[1],
                      " genes before filtering")
    }
    shinydashboardPlus::descriptionBlock(
      number = number,
      numberColor = numberColor,
      numberIcon = numberIcon,
      header = header,
      rightBorder = FALSE
    )
  })
  
  output$tcc_summary <- shiny::renderPrint({
    print(r$tcc)
  })
  
 # toDownload <- shiny::reactiveVal()
  
  
  output$filtering_summary <- shiny::renderUI({
    if (is.null(r$normalized_counts_pre_filter)) {
      numberColor = "red"
      number = "Normalisation needed"
      header = ""
      numberIcon = "fa fa-times"
    }
    else{
      if (is.null(r$normalized_counts)) {
        numberColor = "orange"
        number = "Filtering needed"
        header = ""
        numberIcon = "fa fa-times"
      }
      else{
        numberColor = "olive"
        number = "Done"
        numberIcon = "fa fa-check"
        header = paste(dim(r$normalized_counts)[1],
                       " genes after filtering")
        #toDownload <<- round(r$normalized_counts, 2)
      }
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
  #   distribution plots                                                      ####
  
  # IDEA also implement PCA, maybe 3D, in data exploration
  output$heatmap_preview_norm <- shiny::renderPlot({
    shiny::req(r$raw_counts)
    
    if (!input$preview_norm_distr |
        is.null(r$normalized_counts_pre_filter)) {
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
  
  
#   ____________________________________________________________________________
#   mds                                                                     ####

  
  output$mds_plot <- shiny::renderPlot({
    shiny::req(r$normalized_counts)
    draw_MDS(r$normalized_counts)
  })
  
  
#   ____________________________________________________________________________
#   pca                                                                     ####

  output$pca_plot <- shiny::renderPlot({
    shiny::req(r$normalized_counts)
    draw_PCA(r$normalized_counts)
  })
  
  #   ____________________________________________________________________________
  #   download buttons                                                        ####
  
  output$dl_bttns <- shiny::renderUI({
    shiny::req(r$normalized_counts)
    tagList(
    shiny::fluidRow(col_12(
      shinyWidgets::downloadBttn(
        outputId = ns("download_normalized_counts_csv"),
        label = "Download normalized counts as .csv",
        style = "bordered",
        color = "success"
      ),
      
      shinyWidgets::downloadBttn(
        outputId = ns("download_normalized_counts_RData"),
        label = "Download normalized counts as .RData",
        style = "bordered",
        color = "success"
      )
    ))
    )
    
  })
  
  output$download_normalized_counts_RData <- shiny::downloadHandler(
    filename = function() {
      paste("normalized_counts.RData")
    },
    content = function(file) {
      save(round(r$normalized_counts, 2), file = file)
    }
  )
  
  output$download_normalized_counts_csv <- shiny::downloadHandler(
    filename = function() {
      paste("normalized_counts.csv")
    },
    content = function(file) {
      write.csv(round(r$normalized_counts, 2), file = file, quote = FALSE)
    }
  )
}

## To be copied in the UI
# mod_normalisation_ui("normalisation_ui_1")

## To be copied in the server
# callModule(mod_normalisation_server, "normalisation_ui_1")
