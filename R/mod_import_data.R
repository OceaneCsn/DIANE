

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
    
    shinydashboardPlus::boxPlus(
      title = "Expression file upload",
      width = 4,
      solidHeader = FALSE,
      status = "success",
      collapsible = TRUE,
      closable = FALSE,
      
      shiny::fluidRow(
        col_4(
          shinyWidgets::switchInput(
            ns("use_demo"),
            "Toggle to import your data",
            value = TRUE,
            onLabel = "Demo Arabidopsis data",
            offLabel = "Your dataset",
            onStatus = "success"
            
          )
        )
        ,
        col_8(shiny::uiOutput(ns("gene_ids")))
        
      ),
      
      shiny::uiOutput(ns("org_selection")),
      
      col_8(shiny::radioButtons(
        ns('sep'),
        'Separator : ',
        c(
          Comma = ',',
          Semicolon = ';',
          Tab = '\t'
        ),
        inline = TRUE
      )),
      
      shiny::fluidRow(col_4(shinyWidgets::dropdownButton(
        size = 'xs',
        label = "Input file requirements",
        shiny::includeMarkdown(
          system.file("extdata", "expressionFile.md", package = "DIANE")
        ),
        circle = TRUE,
        status = "primary",
        icon = shiny::icon("question"),
        width = "1200px",
        tooltip = shinyWidgets::tooltipOptions(title = "More details")
      ))),
      
      
      col_12(shiny::fileInput(
        ns('raw_data'),
        'Choose CSV/TXT expression file',
        accept = c(
          'text/csv',
          'text/comma-separated-values,text/plain',
          '.csv',
          '.txt'
        )
      )),
      
      
      
      #   ____________________________________________________________________________
      #   gene infos upload                                                           ####
      
      
      col_8(shiny::radioButtons(ns('sep_gene_info'),
                          'Separator : ',
                          c(Tab = '\t'),
                          inline = TRUE)),
      
      shiny::fluidRow(col_4(shinyWidgets::dropdownButton(
        size = 'xs',
        label = "Gene information file requirements",
        shiny::includeMarkdown(system.file("extdata", "infoFile.md",
                                           package = "DIANE")),
        circle = TRUE,
        status = "primary",
        icon = shiny::icon("question"),
        width = "1200px",
        tooltip = shinyWidgets::tooltipOptions(title = "More details")
      ))),
      
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
      
      shinydashboard::valueBoxOutput(ns("data_dim")),
      shinydashboard::valueBoxOutput(ns("conditions")),
      shinydashboard::valueBoxOutput(ns("samples")),
      
      
      col_4(shiny::uiOutput(ns(
        "variants_summary"
      ))),
      col_4(shiny::uiOutput(ns(
        "organism_summary"
      ))),
      col_4(shiny::uiOutput(ns(
        "gene_info_summary"
      ))),
      
      
      #   ____________________________________________________________________________
      #   seed settings                                                           ####
      
        
        shiny::uiOutput(ns("seed_field")),
      
        
        shinyWidgets::actionBttn(
          ns("change_seed"),
          label = "Change seed",
          style = "material-flat",
          color = "warning"
        ),
        
        
        shinyWidgets::actionBttn(
          ns("set_seed"),
          label = "Set seed",
          style = "material-flat",
          color = "success"
          
        ),
      col_4(shinyWidgets::dropdownButton(
        size = 'xs',
        label = "Input file requirements",
        shiny::includeMarkdown(
          system.file("extdata", "seed.md", package = "DIANE")
        ),
        circle = TRUE,
        status = "primary",
        icon = shiny::icon("question"),
        width = "1200px",
        tooltip = shinyWidgets::tooltipOptions(title = "More details")
      ))
    ),
    
    
    #   ____________________________________________________________________________
    #   Previews                                                                ####
    
    
    shinydashboardPlus::boxPlus(
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
    
    

    
    #   ____________________________________________________________________________
    #   design                                                                  ####
    
    shinydashboardPlus::boxPlus(
      title = "Design and gene information files",
      width = 4,
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
    
    shiny::br(),
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
  
  shiny::observeEvent(input$use_demo, {
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
  #   seed setting                                                            ####
  
  output$seed_field <- shiny::renderUI({
    shiny::req(r$seed)
    shiny::numericInput(
      ns("seed"),
      min = 0,
      max = 2 ^ 8,
      label = "Seed ensuring reproducibility (optional, can be left as default value) :",
      value = r$seed,
      width = "100%"
    )
  })
  
  
  shiny::observeEvent(input$change_seed, {
    r$seed = round(runif(n = 1, min = 0, max = 2 ^ 7))
    shiny::updateNumericInput(session,
                              ns("seed"),
                              value = r$seed)
  })
  
  shiny::observeEvent(input$set_seed, {
    r$seed <- input$seed
    print(paste("changed global seed to", r$seed))
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
    
    ############### checking organism compatibility
    shiny::req(r$organism)
    if (r$organism != "Other") {
      if (!check_IDs(rownames(d), r$organism)) {
        if (r$organism == "Arabidopsis thaliana")
          ex = "AT1G62510.1 or AT1G62510"
        
        if (r$organism == "Homo sapiens")
          ex = "ENSG00000000419"
        
        if (r$organism == "Mus musculus")
          ex = "ENSMUSG00000087910"
        
        if (r$organism == "Drosophilia melanogaster")
          ex = "FBgn0000036"
        
        if (r$organism == "Caenorhabditis elegans")
          ex = "WBGene00000042"
        
        if (r$organism == "Lupinus albus")
          ex = "Lalb_Chr00c02g0404151"
        
        if (r$organism == "Escherichia coli")
          ex = "acpS"
        
        shinyalert::shinyalert(
          "Invalid gene IDs",
          paste(
            "Some or all of the gene IDs in your Gene column do not match
          the expected pattern for the selected organism.
          For",
            r$organism,
            "they should be in the form",
            ex,
            "for example."
          ),
          type = "error"
        )
        #stop()
      }
      shiny::req(check_IDs(rownames(d), r$organism))
    }
    
    r$conditions <-
      stringr::str_split_fixed(colnames(d), "_", 2)[, 1]
    r$splicing_aware <- are_splice_variants(row.names(d))
    r$raw_counts <- d
    #r$gene_info <- gene_info()
    d
  })
  
  #   ____________________________________________________________________________
  #   splicing summary                                                        ####
  output$variants_summary <- shiny::renderUI({
    shiny::req(r$conditions)
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
          paste(
            "The Condition column in your design file should be the experimental
                conditions:",
            paste(r$conditions, collapse = ', ')
          ),
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
  
  
  org_choices <- shiny::reactive({
    choices = c("Arabidopsis thaliana", "Lupinus albus")
    if (requireNamespace("org.Mm.eg.db", quietly = TRUE))
      choices <- c(choices, "Mus musculus")
    
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE))
      choices <- c(choices, "Homo sapiens")
    
    if (requireNamespace("org.Ce.eg.db", quietly = TRUE))
      choices <- c(choices, "Caenorhabditis elegans")
    
    if (requireNamespace("org.Dm.eg.db", quietly = TRUE))
      choices <- c(choices, "Drosophilia melanogaster")
    
    if (requireNamespace("org.EcK12.eg.db", quietly = TRUE))
      choices <- c(choices, "Escherichia coli")
    
    c("Other", choices)
  })
  
  
  output$org_selection <- shiny::renderUI({
    shiny::req(!input$use_demo)
    shiny::selectInput(
      ns("org_select"),
      label = "Your organism :",
      choices = org_choices(),
      selected = "Other"
    )
  })
  
  shiny::observe({
    if (input$use_demo) {
      r$organism <- "Arabidopsis thaliana"
    }
    else{
      shiny::showModal(
        shiny::modalDialog(
          title = "Organism to study",
          shiny::htmlOutput(ns("org_install")),
          shinyWidgets::pickerInput(
            inputId = ns('organism'),
            label = "Choose your organism :",
            choices = c(org_choices()),
            selected = "Other"
          ),
          footer = list(shiny::actionButton(ns("org_chosen"), "OK"))
        )
      )
    }
  })
  
  output$org_install <- shiny::renderText({
    if (!golem::get_golem_options("server_version")) {
      "<b>The organisms listed below are the one detected on the system.</b> <br>
    To use new organisms, please close DIANE and install the corresponding
    package from R or Rstudio consoles.<br>

    <code> if (!requireNamespace(\"BiocManager\", quietly = TRUE))
      install.packages(\"BiocManager\") </code> <br>

    For Human : <code> BiocManager::install(\"org.Hs.eg.db\") </code> <br>
    For Mouse : <code> BiocManager::install(\"org.Mm.eg.db\") </code> <br>
    For Caenorhabditis elegans : <code> BiocManager::install(\"org.Ce.eg.db\") </code> <br>
    For E coli : <code> BiocManager::install(\"org.EcK12.eg.db\") </code> <br>
    For fruit fly : <code> BiocManager::install(\"org.Dm.eg.db\") </code> <br>

    Then, when you launch DIANE again, your organism should appear
    in the following selection menu.

    For now, only Arabidopsis, Human and Mouse are working.
    "
    }
    else{
      "For now, you can choose between all the organisms above"
    }
  })
  
  shiny::observeEvent(input$org_chosen, {
    r$organism <- input$organism
    shiny::removeModal()
    shiny::updateSelectInput(session, "org_select", selected = r$organism)
    
  })
  
  shiny::observe({
    shiny::req(!input$use_demo)
    r$organism <- input$org_select
  })
  
  #   ____________________________________________________________________________
  #   genes info                                                              ####
  
  gene_info <- shiny::reactive({
    req(r$raw_counts)
    req(r$conditions)
    req(r$organism)
    
    if (r$organism != "Other") {
      ids <- rownames(r$raw_counts)
      if (r$splicing_aware) {
        ids <- get_locus(rownames(r$raw_counts))
      }
      if (r$organism == "Lupinus albus") {
        d <-
          DIANE:::lupine$annotation[intersect(ids, rownames(DIANE:::lupine$annotation)), ]
      }
      else{
        d <- get_gene_information(ids, r$organism)
      }
      
    }
    else{
      if (!is.null(input$gene_info_input)) {
        path = input$gene_info_input$datapath
        d <- read.csv(
          sep = input$sep_gene_info,
          path,
          header = TRUE,
          stringsAsFactors = FALSE,
          row.names = "Gene"
        )

        if (length(unique(rownames(d))) < length(rownames(d))) {
          shinyalert::shinyalert(
            "Duplicated genes are not allowed in gene information file",
            type = "error"
          )
          stop()
        }
      }
      else{
        d <- NULL
      }
    }
    d
  })
  ########### table view
  
  output$raw_data_preview <- DT::renderDataTable({
    raw_data()
    shiny::req(r$raw_counts)
    head(r$raw_counts)
  })
  
  ########## matrix preview
  output$heatmap_preview <- shiny::renderPlot({
    shiny::req(r$raw_counts)
    d <- r$raw_counts[rowSums(r$raw_counts) > 0, ]
    
    draw_heatmap(d, title = "Expression data preview")
  })
  
  
  
  #   ____________________________________________________________________________
  #   ValueBoxes summaries                                                    ####
  
  output$gene_ids <- shiny::renderUI({
    shiny::req(r$organism)
    
    if (r$organism == "Other")
      txt <- "No gene ID requirement"
    else {
      data("regulators_per_organism", package = "DIANE")
      txt <- regulators_per_organism[[r$organism]]
    }
    
    shinydashboardPlus::descriptionBlock(
      number = "Expected gene IDs are in the form",
      numberColor = "teal",
      header =  sample(txt, size = 1),
      text = paste("for", r$organism),
      rightBorder = FALSE
    )
  })
  
  
  output$data_dim <- shinydashboard::renderValueBox({
    shiny::req(r$raw_counts)
    
    shinydashboard::valueBox(
      value = dim(r$raw_counts)[1],
      subtitle = "genes",
      color = "aqua",
      width = 4
    )
  })
  output$conditions <- shinydashboard::renderValueBox({
    shiny::req(r$conditions)
    
    shinydashboard::valueBox(value = length((unique(r$conditions))),
                             subtitle = "conditions",
                             color = "teal")
  })
  
  output$samples <- shinydashboard::renderValueBox({
    shiny::req(r$raw_counts)
    shinydashboard::valueBox(
      value = dim(r$raw_counts)[2],
      subtitle = "samples",
      color = "olive"
    )
  })
  
  output$gene_info_summary <- shiny::renderUI({
    shiny::req(r$raw_counts)
    shiny::req(r$organism)
    
    
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
