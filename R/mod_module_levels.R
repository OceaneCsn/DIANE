#' module_levels UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_module_levels_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    shiny::h1("Visualise normalized gene expression profiles"),
    
    shiny::hr(),
    
    shinyalert::useShinyalert(),
    
    shinybusy::add_busy_spinner(
      spin = "self-building-square",
      position = 'top-left',
      margins = c(70, 1200)
    ),
    
    boxPlus(
      title = "Plot settings",
      solidHeader = FALSE,
      status = "success",
      collapsible = TRUE,
      closable = FALSE,
      width = 12,
      

    shiny::uiOutput(ns("gene_choice")),
      
     
    shiny::uiOutput(ns("condition_choice"))
      
    ),
    boxPlus(solidHeader = FALSE,
            status = "success",
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            shiny::plotOutput(ns("expression_plot"), height = "700px"))
 
  )
}
    
#' module_levels Server Function
#'
#' @noRd 
mod_module_levels_server <- function(input, output, session, r){
  ns <- session$ns
  
  
  output$condition_choice <- shiny::renderUI({
    shiny::req(r$normalized_counts, r$conditions)
    
    shinyWidgets::checkboxGroupButtons(
      inputId = ns('input_conditions'),
      label = "Conditions to include to the expression levels plot:",
      choices = unique(r$conditions),
      justified = TRUE,
      checkIcon = list(yes = shiny::icon("ok",
                                         lib = "glyphicon")),
      selected = unique(r$conditions)
      #direction = "vertical"
    )
  })
  
  output$gene_choice <- shiny::renderUI({
    shiny::req(r$normalized_counts)
    
    shiny::textInput(ns("genes"), 
                     label = "Genes to plot, as identified in the Gene 
                     column of expression data, comma separated for several genes :", 
                     width = '100%',
                     value = paste0(sample(rownames(r$normalized_counts), 4), collapse = ','))
  })
  
  
  output$expression_plot <- shiny::renderPlot({
    
    shiny::req(r$normalized_counts, r$conditions, input$genes)
    
    genes <- unlist(strsplit(input$genes, ','))
    shiny::req(length(genes) > 0)
    shiny::req(length(genes) < 10)
    
    shiny::req(sum(genes %in% rownames(r$normalized_counts)) > 0)
    
  
  draw_expression_levels(data.frame(r$normalized_counts),
                           conds = input$input_conditions,
                           genes = genes, gene.name.size = 22)
  })
 
}
    
## To be copied in the UI
# mod_module_levels_ui("module_levels_ui_1")
    
## To be copied in the server
# callModule(mod_module_levels_server, "module_levels_ui_1")
 
