library(shinipsum)
#' import_data UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_import_data_ui <- function(id) {
  ns <- NS(id)
  tagList(
    
    h1("Upload expression data and experimental design"),
    br(),
    p("It should be a matrix like file, with genes
       as rownames and conditions as columns. 
       
       Please speficy the
       replicates unsing the notation _i for the relicate i, placed after the condition name"),
    
    boxPlus(title = "Expression file upload",
                width = 4,
            solidHeader = F,
                status = "primary", 
                #boxToolSize = "xs", 
                collapsible=T,
                closable= F,
                fileInput(
                  ns('raw_data'),
                  'Choose CSV/TXT file', 
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values,text/plain',
                    '.csv',
                    '.txt'
                  )
                ),
                radioButtons(
                  ns('sep'),
                  'Separator : ',
                  c(
                    Comma = ',',
                    Semicolon = ';',
                    Tab = '\t'
                  ),
                  inline = T
                )
    ),
    
    hr(),
    DT::dataTableOutput(ns("raw_data_preview"))
  )
}

#' import_data Server Function
#'
#' @noRd
mod_import_data_server <- function(input, output, session, r) {
  ns <- session$ns
  
  
  raw_data <- reactive({
    validate(need(expr = input$raw_data, message = "No data provided"))
    d <- read.csv(input$raw_data$datapath, sep = input$sep, header = T)
    r$raw_data <- d
  })
  
  output$raw_data_preview <- DT::renderDataTable({
    
    head(raw_data())
    
  })
  
}

## To be copied in the UI
# mod_import_data_ui("import_data_ui_1")

## To be copied in the server
# callModule(mod_import_data_server, "import_data_ui_1")
