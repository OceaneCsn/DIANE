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
        
        col_10(shiny::h4("GENIE3 Gene regulatory network inference")),
        
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
        
        shiny::fluidRow(col_10(shiny::uiOutput(ns("input_genes_net")))),
        shiny::hr(),
        shiny::fluidRow(col_12(shiny::uiOutput(ns("input_conditions_choice_net")))),
        
        shiny::hr(),
        
        shinyWidgets::pickerInput(
          inputId = ns('regulators_picker'),
          label = "Available regulators lists :",
          choices = c("Arabidopsis thaliana - 2192 regulators" = "At"),
        ),
        
        shiny::h5("Or choose your own :"),
        shiny::radioButtons(
          ns('sep_reg'),
          'Separator : ',
          c(
            Comma = ',',
            Semicolon = ';',
            Tab = '\t'
          ),
          inline = TRUE
        ),
        shiny::fileInput(
          ns('TFs_list_input'),
          'Upload custom CSV/TXT regulators list',
          accept = c(
            'text/csv',
            'text/comma-separated-values,text/plain',
            '.csv',
            '.txt'
          )
        ),
        
        shiny::uiOutput(ns("regulators_summary")),
        
        shiny::hr(),
        
        shiny::numericInput(ns("n_cores"), label = "Number of cores available for multithreaded inference :", min = 1, value = 1),
        shiny::numericInput(ns("n_trees"), label = "Number of trees for GENIE3 Random Forests :", min = 1, value = 1),
        
        shiny::fluidRow(
          col_12(shinyWidgets::actionBttn(
            ns("launch_genie_btn"),
            label = "Launch Network Inference",
            color = "success",
            style = 'bordered'
          ))),
        
        
        shiny::hr()
        #shiny::uiOutput(ns("GENIE3_summary"))
        
      )
    )
 
  )
}
    
#' network_inference Server Function
#'
#' @noRd 
mod_network_inference_server <- function(input, output, session, r){
  ns <- session$ns
  
  
  
#   ____________________________________________________________________________
#   deg input select                                                        ####

  output$input_genes_net <- shiny::renderUI({
    shiny::req(r$DEGs)
    shinyWidgets::pickerInput(
      inputId = ns('input_deg_genes_net'),
      label = "Genes composing the network :",
      choices = names(r$DEGs),
      choicesOpt = list(subtext = paste(lengths(r$DEGs), "genes"))
    )
  })
  
  
#   ____________________________________________________________________________
#   conditions selection                                                    ####

  
  output$input_conditions_choice_net <- shiny::renderUI({
    shiny::req(r$conditions)
    
    shinyWidgets::checkboxGroupButtons(
      inputId = ns('input_conditions_net'),
      label = "Conditions to used to infer network edges :",
      choices = unique(r$conditions),
      justified = TRUE,
      checkIcon = list(yes = shiny::icon("ok",
                                         lib = "glyphicon")),
      selected = unique(r$conditions)
    )
  })
  
  
#   ____________________________________________________________________________
#   regulators setting                                                      ####

  shiny::observe({
    shiny::req(r$raw_counts)
    if (r$use_demo) {
      data("demo_data_At", package = "DIANE")
      r$regulators <- demo_data_At$regulators
    }
    else{
      
      if(is.null(input$TFs_list_input)){
        
        organism <- input$regulators_picker
        if(organism == "At"){
          r$regulators <- demo_data_At$regulators
        }
        # else, gérer les autres organismes ici
      }
      else{
        shiny::req(input$TFs_list_input)
        path = input$TFs_list_input$datapath
        
        d <-
          read.csv(
            path,
            sep = input$sep_reg,
            header = FALSE,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        
        r$regulators <- as.vector(d)
        
      }
      
      if (sum(r$regulator %in% length(r$raw_counts)) == 0){
       
        shinyalert::shinyalert(
          "Something is wrong with the chosen regulators :",
          "No regulators were found in the expression data rownames",
          type = "error"
        )
        r$regulators = NULL
        stop()
      }
    }
  })
  #   ____________________________________________________________________________
  #   summaries                                                               ####
  
  
  
  output$regulators_summary <- shiny::renderUI({
    
    if (is.null(r$regulators)) {
      number_color = "orange"
      number = "Please provide regulators"
      header = ""
      number_icon = "fa fa-times"
    }
    else{
      number_color = "olive"
      number = length(unique(r$regulators))
      number_icon = "fa fa-check"
      header = "regulators"
    }
    shinydashboardPlus::descriptionBlock(
      number = number,
      number_color = number_color,
      number_icon = number_icon,
      header = header,
      right_border = FALSE
    )
  })
  
  
#   ____________________________________________________________________________
#   bttn reactive                                                           ####

  
  shiny::observeEvent((input$launch_genie3_btn), {
    shiny::req(r$normalized_counts, input$input_deg_genes_net, r$regulators)
    
    # sets r$networks$genes-conds$mat, 
    
  })
  
}
    
## To be copied in the UI
# mod_network_inference_ui("network_inference_ui_1")
    
## To be copied in the server
# callModule(mod_network_inference_server, "network_inference_ui_1")
 
