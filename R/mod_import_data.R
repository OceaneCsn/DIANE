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
    
    h1("Upload expression data and experimental design"),
    shiny::hr(),
    
    ######################### File upload ###################
    boxPlus(
      title = "Expression file upload",
      width = 4,
      solidHeader = F,
      status = "primary",
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
          Comma = ','
          #Semicolon = ';',
          #Tab = '\t'
        ),
        inline = T
      ),
      
      shinyWidgets::dropdownButton(
        size = 'xs',
        label = "Input file requirements",
        shiny::includeMarkdown('markdown/expressionFile.md'),
        circle = TRUE,
        status = "primary",
        icon = icon("question"),
        width = "600px",
        tooltip = shinyWidgets::tooltipOptions(title = "More details")
      ),
    
    
      fileInput(
        ns('raw_data'),
        'Choose CSV/TXT expression file',
        accept = c(
          'text/csv',
          'text/comma-separated-values,text/plain',
          '.csv',
          '.txt'
        )
      ),
    
    
      shinyWidgets::dropdownButton(
        size = 'xs',
        label = "Design file requirements",
        shiny::includeMarkdown('markdown/designFile.md'),
        circle = TRUE,
        status = "primary",
        icon = icon("question"),
        width = "1200px",
        tooltip = shinyWidgets::tooltipOptions(title = "More details")
      ),
      fileInput(
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
      valueBoxOutput(ns("conditions"))
    ),
    
    ###################### heatmap
    
    boxPlus(
      title = "Preview of the expression matrix",
      width = 4,
      solidHeader = F,
      status = "primary",
      collapsible = T,
      closable = F,
      shiny::plotOutput(ns("heatmap_preview"), height = 550),
      footer = "This might help you visualize the different sequencing depths of your conditions."
    ),
    
    boxPlus(
      title = "Preview of the design",
      width = 3,
      solidHeader = F,
      status = "primary",
      collapsible = T,
      closable = F,
      DT::dataTableOutput(ns("design_preview"), height = 550),
      footer = "footer?"
    ),
    
  
    hr(),
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
  
  raw_data <- shiny::reactive({
    
    if (input$use_demo) {
      path = paste0(r$PATH_TO_DEMO,"/rawData_At.csv")
    }
    else{ 
      req(input$raw_data)
      path = input$raw_data$datapath
    }

      attempt::attempt(
        expr = {d <-
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
        else{stop()}
        }
      )
      r$raw_counts <- d
      d
    })
   

  ######### load the design  
  design <- shiny::reactive({
    if (input$use_demo) {
      path = paste0(r$PATH_TO_DEMO, "/design_At.csv")
    }
    else{
      req(input$design)
      path = input$design$datapath
    }
      
    d <- read.csv(
      path,
      header = T,
      stringsAsFactors = F,
      row.names = "Condition"
    )
    r$design <- d
  })
  
  ########### table view
  
  output$raw_data_preview <- DT::renderDataTable({
    head(raw_data(), n = 8)
  })
  
  ########## matrix preview
  output$heatmap_preview <- shiny::renderPlot({
    shiny::validate(need(expr = raw_data(), message = "Data in the wrong format. Wrong separator? No column named Gene?"))
    d <- raw_data()[rowSums(raw_data()) > 0, ]
    draw_heatmap(d)
  })
  
  ########### data summary
  output$data_dim <- renderValueBox({
    valueBox(
      "Your data",
      paste0(dim(raw_data())[1], " genes, ", dim(raw_data())[2], 
             " samples"),
      color = "teal",
      width=4
    )
  })
  output$conditions <- renderValueBox({
    r$conditions <- stringr::str_split_fixed(colnames(raw_data()), "_", 2)[,1]
    valueBox(
      "Conditions",
      paste0("Your ", length((unique(r$conditions))), " conditions are : ", paste(unique(r$conditions), collapse = ' ')),
      color = "olive",
      width=NULL
    )
  }) 
  
  ######### render design
  output$design_preview <- DT::renderDataTable({
    design()
  })
  
}

## To be copied in the UI
# mod_import_data_ui("import_data_ui_1")

## To be copied in the server
# callModule(mod_import_data_server, "import_data_ui_1")
