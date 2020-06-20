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
    shiny::hr(),
    
    col_4(
      
      
#   ____________________________________________________________________________
#   inference settings                                                      ####

      boxPlus(
        title = "Inference Settings",
        solidHeader = FALSE,
        status = "success",
        collapsible = TRUE,
        closable = FALSE,
        width = 12,
        
        col_10(shiny::h4("GENIE3 Gene regulatory network inference")),
        
        col_2(shinyWidgets::dropdownButton(
          size = 'xs',
          shiny::includeMarkdown(system.file("extdata", "genie3.md", package = "DIANE")),
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

        shiny::uiOutput(ns("reuglators_organism_summary")),
        
        shiny::h5("Or choose your own :"),
        col_2(shinyWidgets::dropdownButton(
          size = 'xs',
          shiny::includeMarkdown(system.file("extdata", "regulatorsFile.md", package = "DIANE")),
          circle = TRUE,
          status = "success",
          icon = shiny::icon("question"),
          width = "600px",
          tooltip = shinyWidgets::tooltipOptions(title = "More details")
        )),
       
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


#   ____________________________________________________________________________
#   thresholding options                                                    ####

    col_6(boxPlus(
      title = "Thresholding settings",
      solidHeader = FALSE,
      status = "success",
      collapsible = TRUE,
      closable = FALSE,
      width = 12,
      
      shiny::fluidRow(col_6(shiny::h4("Thresholding regulatory links?")),
      
      col_6(shinyWidgets::dropdownButton(
        size = 'xs',
        shiny::includeMarkdown(system.file("extdata", "threshold.md", package = "DIANE")),
        circle = TRUE,
        status = "success",
        icon = shiny::icon("question"),
        width = "600px",
        tooltip = shinyWidgets::tooltipOptions(title = "More details")
      ))),
      
      shiny::br(),
      
      shiny::uiOutput(ns("inference_summary")),
      
      shiny::hr(),
      
      
      shiny::fluidRow(col_8(shiny::uiOutput(ns("n_edges_choice"))),
      
      
        col_4(shinyWidgets::actionBttn(
          ns("thr_btn"),
          label = "Threshold",
          color = "success",
          style = 'bordered'
        ))),
      
      shiny::hr(),
      
      shiny::uiOutput(ns("thr_summary")),
      
      visNetwork::visNetworkOutput(ns("net_preview"), height = "650px")
      
    ))


  
  ,shiny::actionButton(ns("browser"), "backdoor")


 
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
  
  
 

  regulators <- shiny::reactive({
    shiny::req(r$raw_counts, r$organism)
    d <- NULL
    if (r$organism != "Other") {
      data("regulators_per_organism", package = "DIANE")
      d <- regulators_per_organism[[r$organism]]
      print(d)
      print(r$organism)
    }
    

    if(!is.null(input$TFs_list_input)){
      path = input$TFs_list_input$datapath
      
      d <-
        read.csv(
          path,
          header = FALSE,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      
      d <- as.vector(d[,1])
    }
    


    if(r$splicing_aware){
      r$aggregated_normalized_counts <- 
        aggregate_splice_variants(data.frame(r$normalized_counts, 
                                             check.names = FALSE))
    }
    
      else {
        if(!is.null(d)){
          
        if (sum(d %in% row.names(r$raw_counts)) == 0){
        
        shinyalert::shinyalert(
          "Something is wrong with the chosen regulators",
          "No regulators were found in the rownames of the expression data",
          type = "error"
        )
        d = NULL
      }
      }
    }

    
    d
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
  
  
  output$reuglators_organism_summary <- shiny::renderUI({

    r$regulators <- regulators()
    
    if(is.null(r$regulators)){
      number_color = "orange"
      number = "Please provide a regulators list"
      number_icon = "fa fa-check"
      header = ""
    }
    else{
      if (r$organism != "Other"){
        number_color = "teal"
        number = length(r$regulators)
        number_icon = "fa fa-check"
        header = paste("regulators provided for", r$organism)
      }
      else{
        number_color = "teal"
        number = length(r$regulators)
        number_icon = "fa fa-check"
        header = "Custom regulators provided"
      }
      
    }
    tagList(shinydashboardPlus::descriptionBlock(
      number = number,
      number_color = number_color,
      text = header,
      right_border = FALSE
    ))
    
  })
  
  output$regulators_intersect_summary <- shiny::renderUI({
    shiny::req(input$input_deg_genes_net, r$regulators, r$DEGs,
               r$networks)
    shiny::req(r$regulators)
    
    if(r$splicing_aware) 
      genes <- get_locus(r$DEGs[[input$input_deg_genes_net]])
    else
      genes <- r$DEGs[[input$input_deg_genes_net]]
      
    tfs <- intersect(genes, r$regulators)

    shinydashboardPlus::descriptionBlock(
      number = length(tfs),
      number_color = "blue",
      text = "Regulators among the input genes",
      right_border = FALSE
    )
  })
  
  output$inference_summary <- shiny::renderUI({
    shiny::req(r$normalized_counts)
    shiny::req(r$DEGs)
    shiny::req(r$DEGs[[input$input_deg_genes_net]])
    shiny::req(r$networks[[input$input_deg_genes_net]])
    shiny::req(input$input_deg_genes_net, r$regulators, r$DEGs)

    if (is.null(r$networks[[input$input_deg_genes_net]]$mat)) {
      number_color = "orange"
      number = "Inference not performed yet"
      header = ""
      number_icon = "fa fa-times"
    }
    else{
      number_color = "olive"
      number = "Inference successfully completed"
      number_icon = "fa fa-check"
      header = "You can now proceed to thresholding"
    }
    shinydashboardPlus::descriptionBlock(
      number = number,
      number_color = number_color,
      text = header,
      right_border = TRUE
    )
  })
  
  output$thr_summary <- shiny::renderUI({
    shiny::req(r$normalized_counts)
    shiny::req(r$DEGs)
    shiny::req(r$DEGs[[input$input_deg_genes_net]])
    shiny::req(r$networks[[input$input_deg_genes_net]])
    shiny::req(r$networks[[input$input_deg_genes_net]]$graph)
    shiny::req(r$networks[[input$input_deg_genes_net]]$nodes)
    shiny::req(r$networks[[input$input_deg_genes_net]]$edges)
    
    number_color = "olive"
    number = "Your network is ready"
    number_icon = "fa fa-check"
    header = "You can explore it in the next tab"
    
    shinydashboardPlus::descriptionBlock(
      number = number,
      number_color = number_color,
      text = header,
      right_border = TRUE
    )
  })
  
  
  
#   ____________________________________________________________________________
#   cores choice                                                            ####

  output$n_cores_choice <- shiny::renderUI({
    
    # assigns either one core if detection fails,
    # either the total number of cores minus one as max
    cpus <- parallel::detectCores()
    if(is.na(cpus)) cpus <- 1
    
    
    shinyWidgets::sliderTextInput(
      inputId = ns("n_cores"),
      label = "Number of cores available for 
                            multithreaded inference :",
      choices = seq(1, cpus),
      grid = TRUE,
      selected = max(1,cpus-1))
    
  })
  
  
  
#   ____________________________________________________________________________
#   thresholding settings                                                   ####
  
  
  output$n_edges_choice <- shiny::renderUI({
    shiny::req(r$normalized_counts)
    shiny::req(r$DEGs)
    shiny::req(r$DEGs[[input$input_deg_genes_net]])
    shiny::req(r$networks[[input$input_deg_genes_net]])
    shiny::req(r$networks[[input$input_deg_genes_net]]$mat)
    proposition = 1.5*length(r$DEGs[[input$input_deg_genes_net]])
    shiny::numericInput(ns("n_edges"), 
                        label = "Number of edges :", 
                        min = 1, value = proposition)
  })
  
  
  
  
  
  
  shiny::observeEvent(input$browser, {
    browser()
  })
  
  
#   ____________________________________________________________________________
#   bttn reactive                                                           ####

  
  shiny::observeEvent((input$launch_genie_btn), {
    shiny::req(r$normalized_counts, input$input_deg_genes_net, r$regulators, r$DEGs)
    
    if(r$splicing_aware) {
      targets <- get_locus(r$DEGs[[input$input_deg_genes_net]])
      data <- r$aggregated_normalized_counts
    }
    else {
      targets <- r$DEGs[[input$input_deg_genes_net]]
      data <- r$normalized_counts
    }
    
    
    if(length(intersect(targets, r$regulators)) < 2 ){
      shinyalert::shinyalert(
        "Not enough regulators provided",
        "GENIE3 requires a minimum of 2 regulators among the input genes to run.
        You coud maybe proceed to a less stringeant differential expression
        analysis to increase the number of input genes.",
        type = "error"
      )
    }
    
    shiny::req(length(intersect(targets, r$regulators)) >= 2)
    
    mat <- network_inference(normalized.count = data, targets = targets, 
                             conds = input$input_conditions_net,
                      regressors = intersect(targets, r$regulators),
                      nTrees = input$n_trees,
                      nCores = input$n_cores)
    
    
    r$networks[[input$input_deg_genes_net]]$mat <- mat
    
    
  })
  
  
  shiny::observeEvent((input$thr_btn), {
    shiny::req(r$normalized_counts)
    shiny::req(r$DEGs)
    shiny::req(r$DEGs[[input$input_deg_genes_net]])
    shiny::req(r$networks[[input$input_deg_genes_net]])
    shiny::req(r$networks[[input$input_deg_genes_net]]$mat)
    
    r$networks[[input$input_deg_genes_net]]$graph <- network_thresholding(
      r$networks[[input$input_deg_genes_net]]$mat, n_edges = input$n_edges)
    
    data <- network_data(r$networks[[input$input_deg_genes_net]]$graph, 
                         r$regulators)
    
    if(!is.null(r$gene_info)){
      data$nodes[,colnames(r$gene_info)] <- 
        r$gene_info[match(data$nodes$id, rownames(r$gene_info)),]
    }
      
    membership <- data$nodes$community
    names(membership) <- data$nodes$id
    
    r$networks[[input$input_deg_genes_net]]$membership <- membership
    
    r$networks[[input$input_deg_genes_net]]$nodes <- data$nodes
    r$networks[[input$input_deg_genes_net]]$edges <- data$edges
    
    r$networks[[input$input_deg_genes_net]]$conditions <- 
      input$input_conditions_net
    
    r$current_network <- input$input_deg_genes_net
    
  })
  
  
#   ____________________________________________________________________________
#   preview                                                                 ####

  
  output$net_preview <- visNetwork::renderVisNetwork({
    shiny::req(r$normalized_counts)
    shiny::req(r$DEGs)
    shiny::req(r$DEGs[[input$input_deg_genes_net]])
    shiny::req(r$networks[[input$input_deg_genes_net]])
    shiny::req(r$networks[[input$input_deg_genes_net]]$graph)
    shiny::req(r$networks[[input$input_deg_genes_net]]$nodes)
    shiny::req(r$networks[[input$input_deg_genes_net]]$edges)
    
    draw_network(nodes = r$networks[[input$input_deg_genes_net]]$nodes,
                 edges = r$networks[[input$input_deg_genes_net]]$edges)
  })
 
}
    
## To be copied in the UI
# mod_network_inference_ui("network_inference_ui_1")
    
## To be copied in the server
# callModule(mod_network_inference_server, "network_inference_ui_1")
 
