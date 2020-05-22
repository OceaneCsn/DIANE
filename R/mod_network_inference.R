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
    shinyalert::useShinyalert(),
    
    shinybusy::add_busy_spinner(
      spin = "self-building-square",
      position = 'top-left',
      margins = c(70, 1200)
    ),
    
    shiny::h1("Network inference"),
    col_3(
      
      
#   ____________________________________________________________________________
#   inference settings                                                      ####

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
        
#   ____________________________________________________________________________
#   regulators input                                                        ####

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
        
        
        shiny::fluidRow(
          col_12(shinyWidgets::actionBttn(
            ns("load_regulators_btn"),
            label = "Load regulators",
            color = "success",
            style = 'bordered'
          ))),
        
        
        shiny::fluidRow(col_4(shiny::uiOutput(ns("input_summary"))),
        col_8(shiny::uiOutput(ns("regulators_intersect_summary")))),

        shiny::hr(),
        
#   ____________________________________________________________________________
#   genie3 launch                                                           ####
        shiny::uiOutput(ns("n_cores_choice")),
       
        shiny::numericInput(ns("n_trees"), 
                            label = "Number of trees for 
                            GENIE3 Random Forests :", 
                            min = 500, value = 1000),
        
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
    ),
shiny::actionButton(ns("browser"), "backdoor"),
 
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

  shiny::observeEvent(input$load_regulators_btn, {
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
      #print(r$regulators)
      #print(head(r$raw_counts))
      #print(row.names(r$raw_counts))
      print(sum(r$regulators %in% row.names(r$raw_counts)))
      if (sum(r$regulators %in% row.names(r$raw_counts)) == 0){
       
        shinyalert::shinyalert(
          "Something is wrong with the chosen regulators :",
          "No regulators were found in the expression data rownames",
          type = "error"
        )
        r$regulators = NULL
        print("done")
        stop()
      }
      
    }
    
  })
  
  #   ____________________________________________________________________________
  #   summaries                                                               ####
  
  
  output$input_summary <- shiny::renderUI({
    shiny::req(input$input_deg_genes_net, r$DEGs)
    if (is.null(input$input_deg_genes_net)) {
      number_color = "orange"
      number = "Please input genes"
      header = ""
      number_icon = "fa fa-times"
    }
    else{
      number_color = "olive"
      number = length(r$DEGs[[input$input_deg_genes_net]])
      number_icon = "fa fa-check"
      header = "input genes"
    }
    shinydashboardPlus::descriptionBlock(
      number = number,
      number_color = number_color,
      text = header,
      right_border = TRUE
    )
  })
  
  output$regulators_intersect_summary <- shiny::renderUI({
    shiny::req(input$input_deg_genes_net, r$regulators, r$DEGs)
    
    tfs <- intersect(r$DEGs[[input$input_deg_genes_net]], r$regulators)

    shinydashboardPlus::descriptionBlock(
      number = length(tfs),
      number_color = "blue",
      text = "Regulators among the input genes",
      right_border = FALSE
    )
  })
  
  
  
#   ____________________________________________________________________________
#   cores choice                                                            ####

  output$n_cores_choice <- shiny::renderUI({
    
    # assigns either one core if detection fails,
    # either the total number of cores minus one as max
    cpus <- parallel::detectCores()
    if(is.na(cpus)) cpus <- 1
    else cpus <- max(cpus-1, 1)
    
    shinyWidgets::sliderTextInput(
      inputId = ns("n_cores"),
      label = "Number of cores available for 
                            multithreaded inference :",
      choices = seq(1, cpus),
      grid = TRUE,
      selected = cpus)
    
  })
  
  
  
  
  
  shiny::observeEvent(input$browser, {
    browser()
  })
  
  
#   ____________________________________________________________________________
#   bttn reactive                                                           ####

  
  shiny::observeEvent((input$launch_genie_btn), {
    shiny::req(r$normalized_counts, input$input_deg_genes_net, r$regulators, r$DEGs)
    
    targets <- r$DEGs[[input$input_deg_genes_net]]
    mat <- network_inference(r$normalized_counts, targets = targets, 
                      regressors = intersect(targets, r$regulators),
                      nTrees = input$n_trees,
                      nCores = input$n_cores)
    
    r$networks[[input$input_deg_genes_net]]$mat <- mat
    print(head(mat))
  })
  
}
    
## To be copied in the UI
# mod_network_inference_ui("network_inference_ui_1")
    
## To be copied in the server
# callModule(mod_network_inference_server, "network_inference_ui_1")
 
