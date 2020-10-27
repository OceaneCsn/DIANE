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
    
    col_4(
      shinydashboardPlus::boxPlus(
        title = "Settings",
        solidHeader = FALSE,
        status = "success",
        collapsible = TRUE,
        closable = FALSE,
        width = 12,
        
        shiny::h4("Estimation of disperion"),
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
        ),
        
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
          style = "material-flat"
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
                        shiny::uiOutput(ns("table_ui"))),
        shiny::tabPanel(
          title = "MA - Vulcano plots",
          
          shinyWidgets::switchInput(
            inputId = ns("MA_vulcano_switch"),
            value = TRUE,
            onLabel = "MA",
            offLabel = "Volcano",
            onStatus = 'success'
          ),
          
          
          shiny::plotOutput(ns("ma_vulcano"), height = "700px")
          
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
                          style = "material-flat"
                        )),
                        col_4(
                          shinyWidgets::radioGroupButtons(ns("draw_go"), 
                                       choices = c("Dot plot", "Enrichment map", "Data table"), 
                                       selected = "Dot plot",
                                       justified = TRUE,
                                       direction = "vertical",
                                       checkIcon = list(
                                         yes = shiny::icon("ok", 
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
                                                                yes = shiny::icon("ok", 
                                                                           lib = "glyphicon"))),
                              shiny::uiOutput(ns("max_go_choice"))),
                        
                        shiny::uiOutput(ns("custom_data_go")),
                        
                        shiny::hr(),
                        
                        shiny::fluidRow(col_12(shiny::uiOutput(ns("go_results"))))
        ),
        shiny::tabPanel(
          title = "Compare genes lists (Venn)",
          shiny::h5("Once more than one differential expression analysis were performed, 
                    you can visualise and compare the different genes lists in a Venn
                    diagram."),
          shiny::uiOutput(ns("venn_lists_choice")),
          shinyWidgets::awesomeRadio(ns("up_down_radio"), label = "Differentially expressed genes to compare :", 
                                     choices = setNames(object = c("All", "Up", "Down"), 
                                                        c("All", "Up-regulated", "Down-regulated")),
                                     inline = TRUE, status = "success"),
          shiny::plotOutput(ns("venn"), height = "700px"),
          shiny::uiOutput(ns("dl_bttn_venn")),
          
          shiny::uiOutput(ns("venn_spec_comp_choice")),
          shiny::uiOutput(ns("venn_spec_comp_bttn"))
        )
      )
      
    )
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
        
        col_6(
        shinyWidgets::radioGroupButtons(
          inputId = ns("reference"),
          label = "Reference",
          choices = unique(r$conditions),
          justified = TRUE, direction = "vertical",
          checkIcon = list(yes = shiny::icon("ok",
                                             lib = "glyphicon"))
        )),
        
        col_6(
        shinyWidgets::radioGroupButtons(
          inputId = ns("perturbation"),
          label = "Perturbation",
          choices = unique(r$conditions),
          selected = unique(r$conditions)[2],
          justified = TRUE,direction = "vertical",
          checkIcon = list(yes = shiny::icon("ok",
                                             lib = "glyphicon"))
        ))
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
    

    shiny::observeEvent((input$deg_test_btn), {
      
      shiny::req(r$tcc)
      if(is.null(r$fit)){
        r_dea$fit <- estimateDispersion(r$tcc)
        r$fit <- r_dea$fit
        
        if(golem::get_golem_options("server_version"))
          loggit::loggit(custom_log_lvl = TRUE,
                       log_lvl = r$session_id,
                       log_msg = "DEA")
      }
      
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
      shiny::req(is.null(r$normalized_counts))

      numberColor = "red"
      number = "Normalisation needed"
      header = ""
      numberIcon = "fa fa-times"
      
      shinydashboardPlus::descriptionBlock(
        number = number,
        numberColor = numberColor,
        numberIcon = numberIcon,
        header = header,
        rightBorder = FALSE
      )
    })
    
    
    output$deg_test_summary <- shiny::renderUI({
      if (is.null(r$normalized_counts)) {
        numberColor = "red"
        number = "Normalisation needed"
        header = ""
        numberIcon = "fa fa-times"
      }
      else{
        if (is.null(r_dea$top_tags)) {
          numberColor = "orange"
          number = "Tests can be performed"
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
            style = "material-flat",            
            color = "success"
          )
        )
        ),
        shiny::hr(),
        shinyWidgets::downloadBttn(
          ns("report"), "Generate html report",
          style = "material-flat", color = "default")
      )
    })
    
    to_dl <- shiny::reactive({
      shiny::req(r_dea$gene_table)
      df <- r_dea$gene_table
      df$Gene_ID <- rownames(r_dea$gene_table)
      df[,!stringr::
           str_detect(colnames(df), "description")]
    })
    
    output$download_table_csv <- shiny::downloadHandler(
      filename = function() {
        paste(paste0("DEGs_", r_dea$ref, "-", r_dea$trt, ".csv"))
      },
      content = function(file) {
        write.table(to_dl(), file = file, row.names = FALSE, sep = ';',
          quote = FALSE)
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
        file.copy(system.file("extdata", "DEA_report.Rmd", package = "DIANE"),
                  tempReport, overwrite = TRUE)
        file.copy(system.file("extdata", "favicon.ico", package = "DIANE"),
                  tempImage, overwrite = TRUE)
        
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
    
    
    
    output$table_ui <- shiny::renderUI({
      if(is.null(r$normalized_counts)) {
        shinydashboardPlus::descriptionBlock(
          number = "Please normalize and filter raw data in normalization tab",
          numberColor = "orange",
          rightBorder = FALSE
        )
      }
      else DT::dataTableOutput(ns("deg_table"))
    })
    
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
    
    #   ____________________________________________________________________________
    #   Venn                                                                    ####
    
    
    output$venn_lists_choice <- shiny::renderUI({
      
      shiny::req(length(r$DEGs) > 1)
      
      shinyWidgets::checkboxGroupButtons(
        inputId = ns("venn_genes"),
        label = "Please select between 2 and 4 lists of genes to show in the Venn diagram :",
        choices = names(r$DEGs),
        justified = TRUE,
        checkIcon = list(yes = shiny::icon("ok",
                                           lib = "glyphicon"))
      )
    })
    
    
    
    venn_list <- shiny::reactive({
      shiny::req(length(input$venn_genes) >= 2 & length(input$venn_genes) <= 4)
      
      if(input$up_down_radio == "All"){
        venn_list <- r$DEGs[input$venn_genes]
      }
      else{
        venn_list <- list()
        for(comp in input$venn_genes){
            if(input$up_down_radio == "Up"){
            venn_list[[comp]] <- r$top_tags[[comp]][
              r$top_tags[[comp]]$logFC > 0 , "genes"]
          }
          else{
            venn_list[[comp]] <- r$top_tags[[comp]][
              r$top_tags[[comp]]$logFC < 0 , "genes"]
          }
        }
      }
      venn_list
    })
    
    output$venn <- shiny::renderPlot({
      shiny::req(venn_list)
      draw_venn(venn_list())
    })
    
    output$dl_bttn_venn <- shiny::renderUI({
      shiny::req(venn_list)
      shiny::req(length(input$venn_genes) >= 2 & length(input$venn_genes) <= 4)
      tagList(
        shiny::fluidRow(col_12(
          shinyWidgets::downloadBttn(
            outputId = ns("download_intersect_venn"),
            label = "Download the intersection of all sets (central part in the Venn diagram)",
            style = "material-flat",
            color = "success"
          )
        )
      )
      )
    })
    
    output$download_intersect_venn <- shiny::downloadHandler(
      filename = function() {
        paste(paste0("Venn_intersection_", paste(input$venn_genes, collapse = "-"), ".csv"))
      },
      content = function(file) {
        write.table(Reduce(intersect, venn_list()), file = file, row.names = FALSE, sep = ';',
                    quote = FALSE, col.names = FALSE)
      }
    )
    
    output$venn_spec_comp_choice <- shiny::renderUI({
      shiny::req(venn_list())
      tagList(
        shiny::selectInput(ns("venn_spec_comp"), label = "Genes specific to a comparison :",
                           choices = input$venn_genes, selected = input$venn_genes[1])
      )
    })
    
    output$venn_spec_comp_bttn <- shiny::renderUI({
      shiny::req(venn_list())
      tagList(
        shinyWidgets::downloadBttn(
          outputId = ns("download_specific_venn"),
          label = paste("Download genes specific to the", input$venn_spec_comp, "list"),
          style = "material-flat",
          color = "success"
        )
      )
    })
    
    
    output$download_specific_venn <- shiny::downloadHandler(
      filename = function() {
        paste(paste0("Venn_specific_to", input$venn_spec_comp, ".csv"))
      },
      content = function(file) {
        write.table(setdiff(venn_list()[[input$venn_spec_comp]], 
                            Reduce(intersect, venn_list())), 
                    file = file, row.names = FALSE, sep = ';',
                    quote = FALSE, col.names = FALSE)
      }
    )
    
    
    #   ____________________________________________________________________________
    #   heatmap                                                               ####
    
    
    output$heatmap_conditions_choice <- shiny::renderUI({
      shiny::req(r$conditions)
      shiny::req(r$top_tags, r_dea$ref, r_dea$trt)
      shiny::req(r$top_tags[[paste(r_dea$ref, r_dea$trt)]])
      shinyWidgets::checkboxGroupButtons(
        inputId = ns("conds_heatmap"),
        label = "Conditions :",
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
          paste0(r_dea$ref, " and ", r_dea$trt)
        )
      )
    })
    
    output$ma_vulcano <- shiny::renderPlot({
      shiny::req(r$top_tags, r_dea$DEGs)
      shiny::req(r$top_tags[[paste(r_dea$ref, r_dea$trt)]])
      draw_DEGs(
        tags = r_dea$tags,
        fdr = input$dea_fdr,
        lfc = input$dea_lfc,
        MA = input$MA_vulcano_switch
      )
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
            r$custom_go <- d
          }
          else{
            shinyalert::shinyalert("Please input Gene to GO term file. ", 
                                   "For now, only Arabidopsis thaliana, mus musculus, and 
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
          
          
      ################# known organisms
          
      }else{
        if (r$organism == "Lupinus albus"){
          genes <- r_dea$top_tags$genes
          background <- rownames(r$normalized_counts)
          

          if (r$splicing_aware) {
            genes <- get_locus(genes)
            background <- get_locus(background)
          }

          GOs <- DIANE:::lupine$go_list
          
          universe <- intersect(background, GOs[,1])
          r_dea$go <- enrich_go_custom(genes, universe, GOs)
          
        }
        
        
        else{
          genes <- r_dea$top_tags$genes
          background <- rownames(r$normalized_counts)
          
          if (r$splicing_aware) {
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
          
          if(r$organism == "Mus musculus"){
            genes <- convert_from_ensembl_mus(genes)
            background <- convert_from_ensembl_mus(background)
            org = org.Mm.eg.db::org.Mm.eg.db
          }
          
          if(r$organism == "Drosophilia melanogaster"){
            genes <- convert_from_ensembl_dm(genes)
            background <- convert_from_ensembl_dm(background)
            org = org.Dm.eg.db::org.Dm.eg.db
          }
          
          if(r$organism == "Caenorhabditis elegans"){
            genes <- convert_from_ensembl_ce(genes)
            background <- convert_from_ensembl_ce(background)
            org = org.Ce.eg.db::org.Ce.eg.db
          }
          
          if(r$organism == "Escherichia coli"){
            genes <- convert_from_ensembl_eck12(genes)
            background <- convert_from_ensembl_eck12(background)
            org = org.EcK12.eg.db::org.EcK12.eg.db
          }
          
          # TODO add check if it is entrez with regular expression here
          shiny::req(length(genes) > 0, length(background) > 0)
          r_dea$go <- enrich_go(genes, background, org = org, GO_type = input$go_type)
        }
        
      }
      if(golem::get_golem_options("server_version"))
        loggit::loggit(custom_log_lvl = TRUE,
                     log_lvl = r$session_id,
                     log_msg = "GO enrichment DEA")
      
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
              style = "material-flat",
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
    
    
  }

## To be copied in the UI
# mod_differential_expression_analysis_ui("differential_expression_analysis_ui_1")

## To be copied in the server
# callModule(mod_differential_expression_analysis_server, "differential_expression_analysis_ui_1")
