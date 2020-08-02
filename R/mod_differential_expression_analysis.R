#' differential_expression_analysis UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_differential_expression_analysis_ui <- function(id) {
  ns <- NS(id)
  tagList(
    shiny::h1("Differential expression analysis"),
    shiny::hr(),
    shinyalert::useShinyalert(),
    shinybusy::add_busy_spinner(
      spin = "self-building-square",
      position = 'top-left',
      margins = c(70, 1200)
    ),
    #   ____________________________________________________________________________
    #   Dispersion estimation                                                   ####
    
    col_3(
      boxPlus(
        title = "Settings",
        solidHeader = FALSE,
        status = "success",
        collapsible = TRUE,
        closable = FALSE,
        width = 12,
        
        shiny::h4("Estimation of disperion : "),
        
        
        shiny::fluidRow(
          col_8(
            shinyWidgets::actionBttn(
              ns("estimate_disp_btn"),
              label = "Launch estimation",
              color = "success",
              style = 'bordered'
            )
          ),
          col_2(
            shinyWidgets::dropdownButton(
              size = 'xs',
              shiny::includeMarkdown(
                system.file("extdata", "edgeR.md", package = "DIANE")
              ),
              circle = TRUE,
              status = "success",
              icon = shiny::icon("question"),
              width = "600px",
              tooltip = shinyWidgets::tooltipOptions(title = "More details")
            )
          )
        ),
        
        
        shiny::hr(),
        shiny::uiOutput(ns("disp_estimate_summary")),
        shiny::hr(),
        
        #   ____________________________________________________________________________
        #   DEG parameters                                                          ####
        
        
        shiny::h4("Conditions to compare for differential analysis : "),
        
        shiny::uiOutput(ns("condition_choices")),
        
        
        shiny::numericInput(
          ns("dea_fdr"),
          min = 0,
          max = 1,
          value = 0.05,
          label = "Adjusted pvalue (fdr)"
        ),
        shiny::numericInput(
          ns("dea_lfc"),
          min = 0,
          max = Inf,
          value = 1,
          label = "Minimum absolute log Fold Change :"
        ),
        
        shinyWidgets::actionBttn(
          ns("deg_test_btn"),
          label = "Detect differentially expressed genes",
          color = "success",
          style = 'bordered'
        ),
        
        shiny::hr(),
        shiny::uiOutput(ns("deg_test_summary")),
        shiny::hr(),
        shiny::uiOutput(ns("deg_number_summary")),
        
        shiny::hr(),
        shiny::br(),
        
        
        shiny::uiOutput(ns("dl_bttns"))
        
        
      )
    ),
    
    #   ____________________________________________________________________________
    #   Visualisation of the results                                            ####
    
    col_8(
      shinydashboard::tabBox(
        title = "Results",
        width = 12,
        shiny::tabPanel(title = "Results table",
                        DT::dataTableOutput(ns("deg_table"))),
        shiny::tabPanel(
          title = "MA - Vulcano plots",
          
          shinyWidgets::switchInput(
            inputId = ns("MA_vulcano_switch"),
            value = TRUE,
            onLabel = "MA",
            offLabel = "Volcano",
            onStatus = 'success'
          ),
          
          
          plotly::plotlyOutput(ns("ma_vulcano"), height = "700px")
          
        ),
        shiny::tabPanel(
          title = "Heatmap",
          shiny::uiOutput(ns("heatmap_conditions_choice")),
          shiny::plotOutput(ns("heatmap"), height = "700px")
        ),
        shiny::tabPanel(title = "Gene Ontology enrichment",
                        
                        col_4(shinyWidgets::actionBttn(
                          ns("go_enrich_btn"),
                          label = "Start GO enrichment analysis",
                          color = "success",
                          style = 'bordered'
                        )),
                        
                        
                        
        
                        
                        col_4(
                          shinyWidgets::radioGroupButtons(ns("draw_go"), 
                                       choices = c("Dot plot", "Enrichment map", "Data table"), 
                                       selected = "Dot plot",
                                       justified = TRUE,
                                       direction = "vertical",
                                       checkIcon = list(
                                         yes = icon("ok", 
                                                    lib = "glyphicon")))
                          
                        ),
                        
                        col_4(shinyWidgets::radioGroupButtons(ns("go_type"), 
                                                              choiceNames = c("Biological process", 
                                                                              "Cellular component", 
                                                                              "Molecular function"),
                                                              choiceValues = c("BP", "CC", "MF"),
                                                              selected = "BP",
                                                              justified = TRUE,
                                                              direction = "vertical",
                                                              checkIcon = list(
                                                                yes = icon("ok", 
                                                                           lib = "glyphicon"))),
                              shiny::uiOutput(ns("max_go_choice"))),
                        
                        shiny::uiOutput(ns("custom_data_go")),
                        
                        shiny::hr(),
                        
                        shiny::fluidRow(col_12(shiny::uiOutput(ns("go_results"))))
        )
        #shiny::verbatimTextOutput(ns("edgeR_summary")),
        
      )
      
    ),
    shiny::actionButton(ns("browser"), "backdoor")
  )
}


#   __________________________________________________________________________________________________________________________________
#   Server                                                                                                                        ####


#' differential_expression_analysis Server Function
#'
#' @noRd
mod_differential_expression_analysis_server <-
  function(input, output, session, r) {
    ns <- session$ns
    
    
    
    r_dea <- shiny::reactiveValues(
      top_tags = NULL,
      DEGs = NULL,
      ref = NULL,
      trt = NULL,
      lfc = NULL,
      fdr = NULL,
      gene_table = NULL
    )

    
    #   ____________________________________________________________________________
    #   Condition choices ui                                                    ####
    
    
    
    output$condition_choices <- shiny::renderUI({
      req(r$conditions)
      tagList(
        shinyWidgets::radioGroupButtons(
          inputId = ns("reference"),
          label = "Reference",
          choices = unique(r$conditions),
          justified = TRUE,
          checkIcon = list(yes = shiny::icon("ok",
                                             lib = "glyphicon"))
        ),
        
        shinyWidgets::radioGroupButtons(
          inputId = ns("perturbation"),
          label = "Perturbation",
          choices = unique(r$conditions),
          selected = unique(r$conditions)[2],
          justified = TRUE,
          checkIcon = list(yes = shiny::icon("ok",
                                             lib = "glyphicon"))
        )
      )
    })
    
    
    
#   ____________________________________________________________________________
#   custom go                                                               ####

    output$custom_data_go <- shiny::renderUI({
      shiny::req(r$organism == "Other")
      shiny::req(is.null(r$custom_go))

      tagList(
        col_2(
          shinyWidgets::dropdownButton(
            size = 'xs',
            shiny::includeMarkdown(
              system.file("extdata", "custom_go.md", package = "DIANE")
            ),
            circle = TRUE,
            status = "success",
            icon = shiny::icon("question"),
            width = "600px",
            tooltip = shinyWidgets::tooltipOptions(title = "More details")
          )
        ),
      col_10(shiny::h4("Your organism is not known to DIANE, but you can provide a matching between 
         gene IDs and GO IDs.")),
      
      
      col_6(shiny::radioButtons(
        ns('sep'),
        
        'Separator : ',
        c(
          Comma = ',',
          Semicolon = ';',
          Tab = '\t'
        ),
        inline = TRUE
      )),
      
      col_6(shiny::fileInput(
        ns('go_data'),
        'Choose CSV/TXT GO terms file',
        accept = c(
          'text/csv',
          'text/comma-separated-values,text/plain',
          '.csv',
          '.txt'
        )
      ))
      )
      
    })

    #   ____________________________________________________________________________
    #   Buttons reactives                                                       ####
    
    
    
    shiny::observeEvent((input$estimate_disp_btn), {
      shiny::req(r$tcc)
      r_dea$fit <- estimateDispersion(r$tcc)
      r$fit <- estimateDispersion(r$tcc)
    })
    
    shiny::observeEvent((input$deg_test_btn), {
      shiny::req(r$fit,
                 input$dea_fdr,
                 input$reference,
                 input$perturbation)
      
      if (input$reference == input$perturbation) {
        shinyalert::shinyalert("You tried to compare the same conditions! 
                               You may need some coffee...",
                               type = "error")
      }
      shiny::req(!input$reference == input$perturbation)
      
      
      r_dea$tags <-
        estimateDEGs(
          r$fit,
          reference = input$reference,
          perturbation = input$perturbation
        )
      
      r_dea$top_tags <-
        r_dea$tags$table[r_dea$tags$table$FDR < input$dea_fdr, ]
      r_dea$top_tags <-
        r_dea$top_tags[abs(r_dea$top_tags$logFC) > input$dea_lfc, ]
      r_dea$DEGs <- r_dea$top_tags$genes
      r_dea$ref <- input$reference
      r_dea$trt <- input$perturbation
      r$DEGs[[paste(r_dea$ref, r_dea$trt)]] <- r_dea$DEGs
      r$top_tags[[paste(r_dea$ref, r_dea$trt)]] <- r_dea$top_tags
      r_dea$go <- NULL
      
      r_dea$lfc <- input$dea_lfc
      r_dea$fdr <- input$dea_fdr
    })
    

    #   ____________________________________________________________________________
    #   Summaries                                                               ####
    
    output$disp_estimate_summary <- shiny::renderUI({
      if (is.null(r$fit)) {
        numberColor = "orange"
        number = "Needed"
        header = ""
        numberIcon = "fa fa-times"
      }
      else{
        numberColor = "olive"
        number = "Done"
        numberIcon = "fa fa-check"
        header = "See EdgeR summary tab for more details"
      }
      shinydashboardPlus::descriptionBlock(
        number = number,
        numberColor = numberColor,
        numberIcon = numberIcon,
        header = header,
        rightBorder = FALSE
      )
    })
    
    
    output$deg_test_summary <- shiny::renderUI({
      if (is.null(r$fit)) {
        numberColor = "red"
        number = "Dispersion estimation needed"
        header = ""
        numberIcon = "fa fa-times"
      }
      else{
        if (is.null(r_dea$top_tags)) {
          numberColor = "orange"
          number = "Test can be performed"
          header = ""
          numberIcon = "fa fa-times"
        }
        else{
          numberColor = "olive"
          number = "Done"
          numberIcon = "fa fa-check"
          header = "See plots and tables for more details"
        }
      }
      shinydashboardPlus::descriptionBlock(
        number = number,
        numberColor = numberColor,
        numberIcon = numberIcon,
        header = header,
        rightBorder = FALSE
      )
    })
    
    output$deg_number_summary <- shiny::renderUI({
      shiny::req(r$top_tags, r_dea$ref, r_dea$trt)
      shiny::req(r$top_tags[[paste(r_dea$ref, r_dea$trt)]])
      
      tagList(
        shiny::fluidRow(
          shinydashboardPlus::descriptionBlock(
            number = sum(r$top_tags[[paste(r_dea$ref, r_dea$trt)]]$logFC > 0),
            numberColor = "olive",
            numberIcon = "fa fa-caret-up",
            header = "up regulated",
            text = "genes",
            rightBorder = TRUE
          ),
          shinydashboardPlus::descriptionBlock(
            number = sum(r$top_tags[[paste(r_dea$ref, r_dea$trt)]]$logFC < 0),
            numberColor = "red",
            numberIcon = "fa fa-caret-down",
            header = "down-regulated",
            text = "genes",
            rightBorder = FALSE
          )
        )
      )
    })
    
    
    #   ____________________________________________________________________________
    #   Dl button                                                               ####
    
    output$dl_bttns <- shiny::renderUI({
      shiny::req(r$top_tags, r_dea$ref, r_dea$trt)
      shiny::req(r$top_tags[[paste(r_dea$ref, r_dea$trt)]])
      tagList(
        shiny::fluidRow(col_12(
          shinyWidgets::downloadBttn(
            outputId = ns("download_table_csv"),
            label = "Download result table as .csv",
            style = "bordered",
            color = "success"
          )
        )
        ),
        shiny::hr(),
        shiny::downloadButton(
          ns("report"), "Generate html report")
      )
    })
    
    output$download_table_csv <- shiny::downloadHandler(
      filename = function() {
        paste(paste0("DEGs_", r_dea$ref, "-", r_dea$trt, ".csv"))
      },
      content = function(file) {
        write.csv(r_dea$top_tags, file = file, quote = FALSE)
      }
    )
    
    #   ____________________________________________________________________________
    #   report                                                                  ####
    
    output$report <- shiny::downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "DEA_report.html",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "DEA_report.Rmd")
        tempImage <- file.path(tempdir(), "favicon.ico")
        file.copy("./R/DEA_report.Rmd", tempReport, overwrite = TRUE)
        file.copy("./inst/app/www/favicon.ico", tempImage, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(r_dea = r_dea)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      }
    )
  
    
    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
    ### download GO                                                             ####
    
    
    output$download_go_table <- shiny::downloadHandler(
      filename = function() {
        paste(paste0("enriched_GOterms", input$reference, '_VS_', 
                    input$perturbation, '_', 
                    input$go_type, ".csv"))
      },
      content = function(file) {
        write.csv(r_dea$go, file = file, quote = FALSE)
      }
    )
    
    
    #   ____________________________________________________________________________
    #   Result plots                                                            ####
    
    output$deg_table <- DT::renderDataTable({
      shiny::req(r$top_tags, r_dea$ref, r_dea$trt)
      shiny::req(r$top_tags[[paste(r_dea$ref, r_dea$trt)]])

      top <- r_dea$top_tags 
      top$Regulation <- ifelse(top$logFC > 0, "Up", "Down")
      
      columns <- c("logFC", "logCPM", "FDR", "Regulation")
      if (!is.null(r$gene_info)) {
        columns <- c(colnames(r$gene_info), columns)
        
        if (r$splicing_aware) ids <- get_locus(rownames(top), unique = FALSE)
        else ids <- rownames(top)
        top[,colnames(r$gene_info)] <- r$gene_info[match(ids, rownames(r$gene_info)),]
      }
      
      r_dea$gene_table <- top[, columns]
      
      DT::formatStyle(
        DT::datatable(top[, columns]),
        columns = c("Regulation"),
        target = c("cell", "row"),
        backgroundColor = DT::styleEqual(c("Up", "Down"), c("#72F02466", c("#FF000035")))
      )
    })
    
    
    output$heatmap_conditions_choice <- shiny::renderUI({
      shiny::req(r$conditions)
      shiny::req(r$top_tags, r_dea$ref, r_dea$trt)
      shiny::req(r$top_tags[[paste(r_dea$ref, r_dea$trt)]])
      shinyWidgets::checkboxGroupButtons(
        inputId = ns("conds_heatmap"),
        label = "Reference",
        choices = unique(r$conditions),
        selected = c(r_dea$ref, r_dea$trt),
        justified = TRUE,
        checkIcon = list(yes = shiny::icon("ok",
                                           lib = "glyphicon"))
      )
    })
    output$heatmap <- shiny::renderPlot({
      shiny::req(input$conds_heatmap, r_dea$DEGs, r$normalized_counts)
      shiny::req(r$top_tags, r_dea$ref, r_dea$trt)
      draw_heatmap(
        data = r$normalized_counts,
        subset = r_dea$DEGs,
        log = TRUE,
        conditions = input$conds_heatmap,
        profiles = TRUE,
        title = paste0(
          "LogCount of differentially expressed genes between : ",
          r_dea$top_tags$comparison
        )
      )
    })
    
    output$ma_vulcano <- plotly::renderPlotly({
      shiny::req(r$top_tags)
      shiny::req(r$top_tags[[paste(r_dea$ref, r_dea$trt)]])
      plotly::ggplotly(draw_DEGs(
        tags = r_dea$tags,
        fdr = input$dea_fdr,
        lfc = input$dea_lfc,
        MA = input$MA_vulcano_switch
      ))
    })
    
    

#   ____________________________________________________________________________
#   GO enrich                                                               ####
    
    shiny::observeEvent((input$go_enrich_btn), {
      shiny::req(r$normalized_counts)
      shiny::req(r_dea)
      shiny::req(r_dea$top_tags)
      
      
      if (r$organism == "Other") {
        
        if(is.null(r$custom_go)){
          if(!is.null(input$go_data)){
            pathName = input$go_data$datapath
            d <- read.csv(
              sep = input$sep,
              file = pathName,
              header = TRUE,
              stringsAsFactors = FALSE
            )
            print(ncol(d))
            r$custom_go <- d
          }
          else{
            shinyalert::shinyalert("Please input Gene to GO term file. ", 
                                   "For now, only Arabidopsis thaliana and 
        Homo sapiens are supported, but you can input your own gene - GO terms matching.",
                                   type = "error")
          }
        }
        shiny::req(r$custom_go)
          if (ncol(r$custom_go) != 2) {
            r$custom_go <- NULL
            shinyalert::shinyalert(
              "Invalid file",
              "It must contain two columns as described.
            Did you correctly set the separator?",
              type = "error"
            )
          }
          
          shiny::req(ncol(r$custom_go) == 2)

          GOs <- r$custom_go
          genes <- r_dea$top_tags$genes
          universe <- intersect(rownames(r$normalized_counts), GOs[,1])
          
          r_dea$go <- enrich_go_custom(genes, universe, GOs)
      }else{
        
        genes <- r_dea$top_tags$genes
        background <- rownames(r$normalized_counts)
        
        if(r$splicing_aware){
          genes <- get_locus(genes)
          background <- get_locus(background)
        }
        
        if(r$organism == "Arabidopsis thaliana"){
          genes <- convert_from_agi(genes)
          background <- convert_from_agi(background)
          org = org.At.tair.db::org.At.tair.db
        }
        
        if(r$organism == "Homo sapiens"){
          genes <- convert_from_ensembl(genes)
          background <- convert_from_ensembl(background)
          org = org.Hs.eg.db::org.Hs.eg.db
        }
        
        # TODO add check if it is entrez with regular expression here
        shiny::req(length(genes) > 0, length(background) > 0)
        r_dea$go <- enrich_go(genes, background, org = org, GO_type = input$go_type)
      }
      
    })
    
#   ____________________________________________________________________________
#   go results                                                              ####
    

    output$go_table <- DT::renderDataTable({
      shiny::req(r_dea$go)
      r_dea$go[,c("Description", "GeneRatio", "BgRatio", "p.adjust")]
    })
    
    output$max_go_choice <- shiny::renderUI({
      shiny::req(r_dea$go)
      shiny::req(input$draw_go =="Dot plot")
     shiny::numericInput(ns("n_go_terms"), 
                                label = "Top number of GO terms to plot :", 
                                min = 1, value = dim(r_dea$go)[1])
    })
    
    output$go_plot <- plotly::renderPlotly({
      shiny::req(r_dea$go)
      max = ifelse(is.na(input$n_go_terms), dim(r_dea$go)[1],input$n_go_terms )
      draw_enrich_go(r_dea$go, max_go = max)
    })
    
    output$go_map_plot <- shiny::renderPlot({
      shiny::req(r_dea$go)
      draw_enrich_go_map(r_dea$go)
    })
    
    output$go_results <- shiny::renderUI({
      
      shiny::req(r_dea$go)
      
      if(nrow(r_dea$go) == 0){
        shinyalert::shinyalert("No enriched GO terms were found",
                               "It can happen if input gene list is not big enough",
                               type = "error")
      }
      
      shiny::req(nrow(r_dea$go) > 0)
      
      if (input$draw_go == "Data table"){
        tagList(
          DT::dataTableOutput(ns("go_table")),
          
            shinyWidgets::downloadBttn(
              outputId = ns("download_go_table"),
              label = "Download enriched GO term as a csv table",
              style = "bordered",
              color = "success"
            )
          
        )
      }
      else{
        if (input$draw_go == "Enrichment map"){
          shiny::plotOutput(ns("go_map_plot"), height = "800px")
        }
        else
          plotly::plotlyOutput(ns("go_plot"), height = "800px")
      }
    })
    
    
    shiny::observeEvent(input$browser, {
      browser()
    })
  }

## To be copied in the UI
# mod_differential_expression_analysis_ui("differential_expression_analysis_ui_1")

## To be copied in the server
# callModule(mod_differential_expression_analysis_server, "differential_expression_analysis_ui_1")
