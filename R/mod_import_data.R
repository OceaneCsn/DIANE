
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
    shinyalert::useShinyalert(),
    shinybusy::add_busy_spinner(
      spin = "self-building-square",
      position = 'top-left',
      margins = c(70, 1100)
    ),
    
    ######################### Title and text
    
    shiny::h1("Upload expression data and experimental design"),
    shiny::hr(),
    
    
    #   ____________________________________________________________________________
    #   File upload                                                             ####
    
    boxPlus(
      title = "Expression file upload",
      width = 4,
      solidHeader = FALSE,
      status = "success",
      collapsible = TRUE,
      closable = FALSE,
      
      shiny::fluidRow(
        col_4(shinyWidgets::prettyCheckbox(
          ns("use_demo"),
          "Use demo data",
          value = TRUE,
          fill = TRUE,
          thick = TRUE,
          status = "success",
          animation = "smooth",
          icon = NULL,
          bigger = TRUE,
          width = "200%"
        )),
        col_8(shinyWidgets::pickerInput(
          inputId = ns('organism'),
          label = "Choose your organism if proposed :",
          choices = c("Arabidopsis thaliana", 
                      "Homo sapiens",
                      "Other")
        ))
      ),
      
      shiny::radioButtons(
        ns('sep'),
        
        'Separator : ',
        c(
          Comma = ',',
          Semicolon = ';',
          Tab = '\t'
        ),
        inline = TRUE
      ),
      
      shinyWidgets::dropdownButton(
        size = 'xs',
        label = "Input file requirements",
        shiny::includeMarkdown(
          system.file("extdata", "expressionFile.md", package = "DIANE")
        ),
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
      #   gene infos upload                                                           ####
      
      shinyWidgets::dropdownButton(
        size = 'xs',
        label = "Gene information file requirements",
        shiny::includeMarkdown(system.file("extdata", "infoFile.md", 
                                           package = "DIANE")),
        circle = TRUE,
        status = "primary",
        icon = shiny::icon("question"),
        width = "1200px",
        tooltip = shinyWidgets::tooltipOptions(title = "More details")
      ),
      shiny::radioButtons(
        ns('sep_gene_info'),
        'Separator : ',
        c(
          Tab = '\t'
        ),
        inline = TRUE
      ),
      shiny::fileInput(
        ns('gene_info_input'),
        'Choose CSV/TXT gene information file (optional)',
        accept = c(
          'text/csv',
          'text/comma-separated-values,text/plain',
          '.csv',
          '.txt'
        )
      ),
      
      
      valueBoxOutput(ns("data_dim")),
      valueBoxOutput(ns("conditions")),
      valueBoxOutput(ns("samples")),
      
      
      col_4(shiny::uiOutput(ns("variants_summary"))),
      col_4(shiny::uiOutput(ns("organism_summary"))),
      col_4(shiny::uiOutput(ns("gene_info_summary")))
    ),
    
    
    #   ____________________________________________________________________________
    #   Previews                                                                ####
    
    
    boxPlus(
      title = "Preview of the expression matrix",
      width = 4,
      solidHeader = FALSE,
      status = "success",
      collapsible = TRUE,
      closable = FALSE,
      shiny::plotOutput(ns("heatmap_preview"), height = 550),
      footer = "This might help you visualize the general aspect of the data and different sequencing depths 
      of your conditions."
    ),
    
    boxPlus(
      title = "Design file",
      width = 3,
      solidHeader = FALSE,
      status = "success",
      collapsible = TRUE,
      closable = FALSE,
      shiny::radioButtons(
        ns('sep_design'),
        
        'Separator : ',
        c(
          Comma = ',',
          Semicolon = ';',
          Tab = '\t'
        ),
        
        inline = TRUE
      ),
      shinyWidgets::dropdownButton(
        size = 'xs',
        label = "Design file requirements",
        shiny::includeMarkdown(system.file("extdata", "designFile.md", 
                                           package = "DIANE")),
        circle = TRUE,
        status = "primary",
        icon = shiny::icon("question"),
        width = "550px",
        tooltip = shinyWidgets::tooltipOptions(title = "More details")
      ),
      shiny::fileInput(
        ns('design'),
        'Choose CSV/TXT design file (optional)',
        accept = c(
          'text/csv',
          'text/comma-separated-values,text/plain',
          '.csv',
          '.txt'
        )
      ),
      DT::dataTableOutput(ns("design_preview")),
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
  
  # resets the global reactive variables that were maybe already created
  # when demo usage is toggled :

  shiny::observeEvent(input$use_demo,{
    r$raw_counts = NULL
    r$normalized_counts = NULL
    r$normalized_counts_pre_filter = NULL
    r$conditions = NULL
    r$design = NULL
    r$DEGs = list()
    r$tcc = NULL
    r$clusterings = list()
    r$current_comparison = NULL
    r$current_network = NULL
    r$top_tags = list()
    r$fit = NULL
    r$regulators = NULL
    r$use_demo = input$use_demo
    r$splicing_aware = NULL
    r$gene_info = NULL
    r$organism = NULL
    r$custom_go = NULL
  })
  
  #   ____________________________________________________________________________
  #   expression file                                                         ####

  
  raw_data <- shiny::reactive({
    
    if (input$use_demo) {
      r$use_demo = input$use_demo
      data("abiotic_stresses", package = "DIANE")
      d <- abiotic_stresses$raw_counts
    }
    else{
      req(input$raw_data)
      path = input$raw_data$datapath
      
      r$raw_counts = NULL
      r$normalized_counts = NULL
      r$normalized_counts_pre_filter = NULL
      r$conditions = NULL
      r$design = NULL
      r$DEGs = list()
      r$tcc = NULL
      r$clusterings = list()
      r$current_comparison = NULL
      r$current_network = NULL
      r$top_tags = list()
      r$fit = NULL
      r$regulators = NULL
      r$use_demo = input$use_demo
      r$splicing_aware = NULL
      r$gene_info = NULL
      r$organism = NULL
      r$custom_go = NULL
      
      d <-
        read.csv(
          path,
          sep = input$sep,
          header = TRUE,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      if ("Gene" %in% colnames(d)) {
        d <-
          read.csv(
            path,
            sep = input$sep,
            header = TRUE,
            stringsAsFactors = FALSE,
            row.names = "Gene",
            check.names = FALSE
          )
      }
      else{
        #bug here
        shinyalert::shinyalert(
          "Invalid input data",
          "Did you correctly set the separator?
                               Does your data contains a column named \"Gene\"?",
          type = "error"
        )
        stop()
      }
    }
    
    if (length(unique(colnames(d))) < length(colnames(d))) {
      shinyalert::shinyalert(
        "Invalid rownames",
        "Please specify unique rownames, in the form condition_replicateNumber",
        type = "error"
      )
      stop()
    }
    
    r$conditions <-
      stringr::str_split_fixed(colnames(d), "_", 2)[, 1]
    
    r$splicing_aware <- are_splice_variants(row.names(d))
    r$raw_counts <- d
    d
  })
  
  # TODO problem window id wrong design
  
  
#   ____________________________________________________________________________
#   splicing summary                                                        ####
  output$variants_summary <- shiny::renderUI({
    shiny::req(!is.null(r$splicing_aware))
    if (r$splicing_aware) {
      numberColor = "blue"
      number = "Alternatifve splicing aware"
      header = "gene identifiers"
    }
    else{
      numberColor = "blue"
      number = "No alternatifve splicing information"
      header = "in gene identifiers"
    }
    shinydashboardPlus::descriptionBlock(
      number = number,
      numberColor = numberColor,
      text = header,
      rightBorder = TRUE,
    )
  })
  
  
  #   ____________________________________________________________________________
  #   design loading                                                          ####
  
  design <- shiny::reactive({
    
    if (input$use_demo) {
      data("abiotic_stresses", package = "DIANE")
      d <- abiotic_stresses$design
    }
    else{
      req(r$conditions)
      req(input$design)
      path = input$design$datapath
      d <- read.csv(
        sep = input$sep_design,
        path,
        header = TRUE,
        stringsAsFactors = FALSE,
        row.names = "Condition"
      )
      if (sum(rownames(d) %in% r$conditions) < dim(d)[1]) {
        shinyalert::shinyalert(
          "Invalid design rownames...",
          paste("The conditions in your design file should be the experimental 
                conditions:", 
                paste(r$conditions, collapse = ', ')),
          type = "error"
        )
        stop()
      }
    }
    
    r$design <- d
    d
  })

  
  #   ____________________________________________________________________________
  #   organism                                                                ####
  
  organism <- shiny::reactive({
    if (input$use_demo) {
      d <- "Arabidopsis thaliana"
    }
    else{
      d <- input$organism
    }
    d
  })
  
  
  #   ____________________________________________________________________________
  #   genes info                                                              ####
  
  gene_info <- shiny::reactive({
    req(r$raw_counts)
    if (r$organism != "Other") {
      ids <- rownames(r$raw_counts)
      if(r$splicing_aware){
        ids <- get_locus(rownames(r$raw_counts))
      }
      d <- get_gene_information(ids, r$organism)
    }
    else{
        if(!is.null(input$gene_info_input)){
          path = input$gene_info_input$datapath
          d <- read.csv(
            sep = input$sep_gene_info,
            path,
            header = TRUE,
            stringsAsFactors = FALSE,
            row.names = "Gene",
            quote = ""
          )
        }
      else{
        d <- NULL
      }
    }
    d
  })
  ########### table view
  
  output$raw_data_preview <- DT::renderDataTable({
    head(raw_data(), n = 6)
  })
  
  ########## matrix preview
  output$heatmap_preview <- shiny::renderPlot({
    shiny::req(r$raw_counts)
    d <- raw_data()[rowSums(raw_data()) > 0,]
    
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
    shiny::req(r$conditions)
    valueBox(value = length((unique(r$conditions))),
             subtitle = "conditions",
             color = "teal")
  })
  
  output$samples <- renderValueBox({
    valueBox(value = length(colnames(raw_data())),
             subtitle = "samples",
             color = "olive")
  })
  
  output$gene_info_summary <- shiny::renderUI({
    ######## setting gene info here
    r$gene_info <- gene_info()
    
    if (is.null(r$gene_info)) {
      numberColor = "orange"
      number = "No additional gene data provided"
      header = ""
      numberIcon = "fa fa-times"
    }
    else{
      numberColor = "olive"
      number = "Additional gene data available"
      numberIcon = "fa fa-check"
      header = paste(colnames(r$gene_info), collapse = ', ')
    }
    shinydashboardPlus::descriptionBlock(
      number = number,
      numberColor = numberColor,
      numberIcon = numberIcon,
      text = header,
      rightBorder = FALSE
    )
  })
  
  output$organism_summary <- shiny::renderUI({
    ######## setting organism here
    r$organism <- organism()
    
    shiny::req(r$organism)
    
    shinydashboardPlus::descriptionBlock(
      number = r$organism,
      numberColor = "teal",
      text = "organism database",
      rightBorder = FALSE
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
