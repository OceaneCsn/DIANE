source("R/fct_normalisation.R")
source("R/fct_heatmap.R")
library(shinybusy)

#' normalisation UI Function
#'
#' @description A shiny Module for data filtering and normalisation.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @importFrom shinyWidgets actionBttn
#' @import shinydashboard
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_normalisation_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    shiny::h1("Data filtering and normalisation"),
    shiny::h2("Because low count genes and differences in sequencing depths are a true source of bias, we want to perform some data cleaning and transformation."),
    shiny::hr(),
    col_6(
      boxPlus(
        title = "Normalisation",
        solidHeader = F,
        status = "primary",
        collapsible = T,
        closable = F,
        h3("Normalisation step using TCC"),
        
        shinyWidgets::actionBttn(ns("normalize_btn"), label = "Normalize", color = "primary",
                                 style = 'bordered'))
        
      ,
      boxPlus(
        title = "Low counts filtering",
        solidHeader = F,
        status = "primary",
        collapsible = T,
        closable = F,
        
        shiny::h5( "Minimal gene count sum accross conditions : "),
        col_8(shiny::numericInput(ns("low_counts_filter"), 
                                  min = 0, value = 40, 
                                  label = NULL)),
        col_4(shinyWidgets::actionBttn(ns("use_SumFilter"), label = "Filter", color = "primary",
                                       style = 'bordered')),
        shiny::br(),
        
        shinyWidgets::actionBttn(ns("use_HTSFilter"), label = "Use HST Filter, an entropy based filtering method",
                                 style = 'bordered', color = 'success')
        #valueBoxOutput(ns("filtering_summary")),
        #footer = "This might help you visualize the different sequencing depths of your conditions."
      )
    ),
    col_6(h3("Plots"),
          DT::dataTableOutput(ns("norm_factor") )
    )
    
    
  )
}
    
#' normalisation Server Function
#'
#' @noRd 
mod_normalisation_server <- function(input, output, session, r){
  ns <- session$ns
  
  shiny::observeEvent((input$use_HTSFilter), {
    
      print("HTS")
  })
  
  
  
  shiny::observeEvent((input$use_SumFilter), {
      print("Sum")
  })
  
  shiny::observeEvent((input$normalize_btn), {
    print("Normalize!!")
    
    shinybusy::show_modal_spinner(text = "TCC normalisation ongoing...", spin= "self-building-square") # show the modal window
    
    
    # progressSweetAlert(
    #   session = session, id = "myprogress",
    #   title = "Work in progress",
    #   display_pct = TRUE, value = 50
    # )
    norm <- normalize(r$raw_counts, r$conditions)
    r$normalized_counts <- norm$normalized.counts
    r$norm_factors <- norm$norm_factors
    
    #closeSweetAlert(session = session)
    shinybusy::remove_modal_spinner()
    sendSweetAlert(
      session = session,
      title =" Normalisation completed !",
      type = "success"
    )
  })
  
  
  output$norm_factor <- DT::renderDataTable({
    req(r$norm_factors)
    data.frame(r$norm_factor)
  })
  
  
  # output$data_dim <- renderValueBox({
  #   req(r$filtered_raw_counts)
  #   valueBox(
  #     "Your data",
  #     paste0(dim(raw_data())[1], " genes, ", dim(raw_data())[2], 
  #            " samples"),
  #     color = "teal",
  #     width=4
  #   )
  # })
 
}
    
## To be copied in the UI
# mod_normalisation_ui("normalisation_ui_1")
    
## To be copied in the server
# callModule(mod_normalisation_server, "normalisation_ui_1")
 
