#' network_inference UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_network_inference_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::h1("Network inference"),
    
    col_3(
      boxPlus(
        title = "Settings",
        solidHeader = FALSE,
        status = "success",
        collapsible = TRUE,
        closable = FALSE,
        width = 12,
        
        col_10(shiny::h4("GENIE3 Gene regulatory network inferenceinsta")),
        
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
    )
 
  )
}
    
#' network_inference Server Function
#'
#' @noRd 
mod_network_inference_server <- function(input, output, session){
  ns <- session$ns
 
}
    
## To be copied in the UI
# mod_network_inference_ui("network_inference_ui_1")
    
## To be copied in the server
# callModule(mod_network_inference_server, "network_inference_ui_1")
 
