library(stringr)
library(shinyWidgets)

source("R/fct_heatmap.R")

#' import_data UI Function
#'
#' @description A shiny Module to import expression data.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @importFrom shinydashboardPlus boxPlus
#' @importFrom shinydashboard valueBoxOutput
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_import_data_ui <- function(id) {
  ns <- NS(id)
  tagList(
    
    ######################### Title and text
    
    shiny::h1("Upload expression data and experimental design"),
    shiny::hr(),
    

#   ____________________________________________________________________________
#   File upload                                                             ####

    boxPlus(
      title = "Expression file upload",
      width = 4,
      solidHeader = F,
      status = "success",
      collapsible = T,
      closable = F,
      shinyWidgets::prettyCheckbox(ns("use_demo"),
                     "Use demo data",
                     value = TRUE,
                     fill = T,
                     thick = TRUE,
                     status = "success",
                     animation = "smooth",
                     icon = NULL,
                     bigger=TRUE,
                     width = "200%"), 
      radioButtons(
        ns('sep'),
          
        'Separator : ',
        c(
          Comma = ',',
          Semicolon = ';',
          Tab = '\t'
        ),
        inline = T
      ),
      
      shinyWidgets::dropdownButton(
        size = 'xs',
        label = "Input file requirements",
        shiny::includeMarkdown(system.file("extdata", "expressionFile.md", package = "DIANE")),
        circle = TRUE,
        status = "primary",
        icon = shiny::icon("question"),
        width = "600px",
        tooltip = shinyWidgets::tooltipOptions(title = "More details")
      ),
    
    
      shiny::fileInput(
        ns('raw_data'),
        'Choose CSV/TXT expression file',
        accept = c(
          'text/csv',
          'text/comma-separated-values,text/plain',
          '.csv',
          '.txt'
        )
      ),
    
#   ____________________________________________________________________________
#   design upload                                                           ####

    
      shinyWidgets::dropdownButton(
        size = 'xs',
        label = "Design file requirements",
        shiny::includeMarkdown(system.file("extdata", "designFile.md", package = "DIANE")),
        circle = TRUE,
        status = "primary",
        icon = shiny::icon("question"),
        width = "1200px",
        tooltip = shinyWidgets::tooltipOptions(title = "More details")
      ),
      shiny::fileInput(
        ns('design'),
        'Choose CSV/TXT design file',
        accept = c(
          'text/csv',
          'text/comma-separated-values,text/plain',
          '.csv',
          '.txt'
        )
      ),
     
      valueBoxOutput(ns("data_dim")),
      valueBoxOutput(ns("conditions")),
      valueBoxOutput(ns("samples"))
    ),
    

#   ____________________________________________________________________________
#   Previews                                                                ####

    
    boxPlus(
      title = "Preview of the expression matrix",
      width = 4,
      solidHeader = F,
      status = "success",
      collapsible = T,
      closable = F,
      shiny::plotOutput(ns("heatmap_preview"), height = 550),
      footer = "This might help you visualize the different sequencing depths of your conditions."
    ),
    
    boxPlus(
      title = "Preview of the design",
      width = 3,
      solidHeader = F,
      status = "success",
      collapsible = T,
      closable = F,
      DT::dataTableOutput(ns("design_preview"), height = 550),
      footer = "Describe the levels of each factors for your conditions"
    ),
    
  
    shiny::hr(),
    DT::dataTableOutput(ns("raw_data_preview"))
  )
}

#' import_data Server Function
#' @importFrom utils read.csv
#' @importFrom utils head
#' @importFrom stats heatmap
#' @importFrom shinydashboard renderValueBox
#' @importFrom shinydashboard valueBox
#' @noRd
mod_import_data_server <- function(input, output, session, r) {
  ns <- session$ns
  
  
#   ____________________________________________________________________________
#   expression file                                                         ####

  raw_data <- shiny::reactive({
    if (input$use_demo) {
      load(system.file("extdata", "raw_counts_At.RData", package = "DIANE"))
      d <- raw_data_At
    }
    else{ 
      
      # reset the global reactive variables that were maybe already created :
      
      r$raw_counts = NULL
      r$normalized_counts = NULL
      r$normalized_counts_pre_filter = NULL
      r$conditions = NULL
      r$design = NULL
      r$DEGs = list()
      r$tcc = NULL
      
     
      req(input$raw_data)
      path = input$raw_data$datapath

      d <-
        read.csv(
          path,
          sep = input$sep,
          header = T,
          stringsAsFactors = F
        )
      
      if ("Gene" %in% colnames(d)) {
        d <-
          read.csv(
            path,
            sep = input$sep,
            header = T,
            stringsAsFactors = F,
            row.names = "Gene"
          )
      }
      else{
        #bug here
        shinyalert::shinyalert("Invalid input data...", "Did you correctly set the separator? 
                               Does your data contains a column named \"Gene\"?", type = "error")
        stop()
      }
    }
    r$raw_counts <- d
    d
    })
   


#   ____________________________________________________________________________
#   design loading                                                          ####

  design <- shiny::reactive({
    if (input$use_demo) {
      load(system.file("extdata", "design_At.RData", package = "DIANE"))
      d <- design_At
    }
    else{
      req(input$design)
      path = input$design$datapath
      d <- read.csv(sep = input$sep,
        path,
        header = T,
        stringsAsFactors = F,
        row.names = "Condition"
      )
    }
    r$design <- d
  })
  
  ########### table view
  
  output$raw_data_preview <- DT::renderDataTable({
    head(raw_data(), n = 8)
  })
  
  ########## matrix preview
  output$heatmap_preview <- shiny::renderPlot({
    shiny::req(r$raw_counts)
    d <- raw_data()[rowSums(raw_data()) > 0, ]
    draw_heatmap(d)
  })
  

  
#   ____________________________________________________________________________
#   ValueBoxes summaries                                                    ####

  output$data_dim <- renderValueBox({
    valueBox(
      
      value = dim(raw_data())[1],
      subtitle = "genes",
      color = "aqua",
      width = 4
    )
  })
  output$conditions <- renderValueBox({
    r$conditions <- stringr::str_split_fixed(colnames(raw_data()), "_", 2)[,1]
    valueBox(
      value = length((unique(r$conditions))),
      subtitle = "conditions",
      color = "teal"
    )
  })  
  
  output$samples <- renderValueBox({
    valueBox(
      value = length(colnames(raw_data())),
      subtitle = "samples",
      color = "olive"
    )
  }) 
  
  ######### render design
  output$design_preview <- DT::renderDataTable({
    DT::datatable(design())
  })
  
}

## To be copied in the UI
# mod_import_data_ui("import_data_ui_1")

## To be copied in the server
# callModule(mod_import_data_server, "import_data_ui_1")
