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
      shinydashboardPlus::box(
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
            shiny::includeMarkdown(system.file("extdata", "edgeR.md", package = "DIANE")),
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
        
        shiny::htmlOutput(ns("condition_choices_visualisation_2")),
        
        shiny::numericInput(
          ns("dea_fdr"),
          min = 0,
          max = 1,
          value = 0.05,
          label = "Adjusted pvalue ( FDR )"
        ),
        shiny::numericInput(
          ns("dea_lfc"),
          min = 0,
          max = Inf,
          value = 1,
          label = "Absolute Log Fold Change ( Log2 ( Perturbation / Reference ) ) :"
        ),
        
        shiny::uiOutput(ns("multiple_DE_parameters_ui")),
        
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
        
        #   ____________________________________________________________________________
        #   Go enrichment                                                           ####
        
        shiny::tabPanel(
          title = "Gene Ontology enrichment",
          
          shinyWidgets::radioGroupButtons(
            ns("up_down_go_radio"),label = "Genes to study :",
            choices = c("All", "Up-regulated", "Down-regulated"),
            selected = "All",
            direction = "horizontal",
            checkIcon = list(yes = shiny::icon("ok",
                                               lib = "glyphicon"))
          ),
          
          col_4(
            shinyWidgets::actionBttn(
              ns("go_enrich_btn"),
              label = "Start GO enrichment analysis",
              color = "success",
              style = "material-flat"
            )
          ),
          col_4(
            shinyWidgets::radioGroupButtons(
              ns("draw_go"),
              choices = c("Dot plot", "Enrichment map", "Data table"),
              selected = "Dot plot",
              justified = TRUE,
              direction = "vertical",
              checkIcon = list(yes = shiny::icon("ok",
                                                 lib = "glyphicon"))
            )
            
          ),
          
          col_4(
            shinyWidgets::radioGroupButtons(
              ns("go_type"),
              choiceNames = c(
                "Biological process",
                "Cellular component",
                "Molecular function"
              ),
              choiceValues = c("BP", "CC", "MF"),
              selected = "BP",
              justified = TRUE,
              direction = "vertical",
              checkIcon = list(yes = shiny::icon("ok",
                                                 lib = "glyphicon"))
            ),
            shiny::uiOutput(ns("max_go_choice"))
          ),
          
          shiny::uiOutput(ns("custom_data_go")),
          
          shiny::hr(),
          
          shiny::fluidRow(col_12(shiny::uiOutput(ns(
            "go_results"
          ))))
        ),
        shiny::tabPanel(
          title = "Compare genes lists (Venn)",
          shiny::h5(
            "Once more than one differential expression analysis were performed,
                    you can visualise and compare the different genes lists in a Venn
                    diagram."
          ),
          shiny::uiOutput(ns("venn_lists_choice_2")),
          shiny::plotOutput(ns("venn"), height = "700px"),
          shiny::uiOutput(ns("venn_spec_comp_choice_2")),
          shiny::uiOutput(ns("venn_spec_comp_bttn_2"))
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
          shinyWidgets::checkboxGroupButtons(
            inputId = ns("reference"),
            label = "Reference",
            choices = unique(r$conditions),
            selected = unique(r$conditions)[1],
            justified = TRUE,
            direction = "vertical",
            checkIcon = list(yes = shiny::icon("ok",
                                               lib = "glyphicon"))
          )
        ),
        
        col_6(
          shinyWidgets::checkboxGroupButtons(
            inputId = ns("perturbation"),
            label = "Perturbation",
            choices = unique(r$conditions),
            selected = unique(r$conditions)[2],
            justified = TRUE,
            direction = "vertical",
            checkIcon = list(yes = shiny::icon("ok",
                                               lib = "glyphicon"))
          )
        )
      )
    })
    
    output$condition_choices_visualisation_2 <- shiny::renderText({
      # req(input$reference, input$perturbation)
      if(is.null(input$reference)){
        reference_text = "<span style=\"color: #A52014\">NOTHING SELECTED</span>"
      } else {
        reference_text <- ifelse(test = length(input$reference) == 1,
                                 yes = input$reference,
                                 no = paste0(input$reference, collapse = " + "))
      }
      if(is.null(input$perturbation)){
        perturbation_text = "<span style=\"color: #A52014\">NOTHING SELECTED</span>"
      } else {
        perturbation_text <- ifelse(test = length(input$perturbation) == 1,
                                    yes = input$perturbation,
                                    no = paste0(input$perturbation, collapse = " + "))
      }
      comparison_type <- ifelse(length(input$reference) > 1 | length(input$perturbation) > 1,
                                yes = "Multiple comparison",
                                no = "Simple comparison")
      as.character(
        paste0(
          "<div style=\"text-align: center; font-family: 'Arial';\"> <h4 style=\"text-align: center\"><b>",comparison_type,"</b></h4>",
          shiny::column(6,
                        HTML(paste0("<b>Reference</b><br>",reference_text,""))
          ),
          shiny::column(6,
                        HTML(paste0("<b>Perturbation</b><br>",perturbation_text, ""))
          )," </div>"
        )
      )
    })
    
    
    output$multiple_DE_parameters_ui <- shiny::renderUI({
      shiny::req(!is.null(r$normalized_counts))
      if (length(input$reference) > 1 ||
          length(input$perturbation) > 1) {
        shiny::tagList(shiny::HTML(
          paste0(
            "<span style=\"display: inline-block; margin-right: 5px; max-width: 100%;
            # margin-bottom: 5px; \"><h5 style=\"display: inline-block; font-weight: 700;\"> Multiple DE comparison method  </h5>",
            as.character(
              shinyWidgets::dropdownButton(
                size = "xs",
                shiny::HTML(
                  "<p>There is two methods for multiple condition comparison.<br><b>Mean of conditions</b> :
                  the default method. Take the mean of all selected coniditions in order to perform differential expression analysis.<br>
                  <b>Mean of condition and same orientation DE</b> :
                  Fist perform differentiall expression using mean of conditions.
                  Then, perform every single possible differential expression analysis using the selected conditions,
                  and keep only genes that are differentially expressed with the same orientation in all comparison.</p>"
                ),
                # shiny::includeMarkdown(system.file("extdata", "edgeR.md", package = "DIANE")),
                circle = TRUE,
                status = "success",
                icon = shiny::icon("question"),
                width = "400px",
                tooltip = shinyWidgets::tooltipOptions(title = "More details"),
                inline = TRUE
              )
            ),
            "</span>"
          )
        ),
        shinyWidgets::pickerInput(
          inputId = ns("multiple_DE_parameters_selection"),
          # label = "Multiple DE comparison method",
          choices = c("Mean of conditions" = FALSE, "Mean of condition and same orientation DE" = TRUE)
        ))
      } else {
        print("Noooo")
        NULL
      }
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
            shiny::includeMarkdown(system.file("extdata", "custom_go.md", package = "DIANE")),
            circle = TRUE,
            status = "success",
            icon = shiny::icon("question"),
            width = "600px",
            tooltip = shinyWidgets::tooltipOptions(title = "More details")
          )
        ),
        col_10(
          shiny::h4(
            "Your organism is not known to DIANE, but you can provide a matching between
         gene IDs and GO IDs."
          )
        ),
        
        
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
        
        col_6(
          shiny::fileInput(
            ns('go_data'),
            'Choose CSV/TXT GO terms file',
            accept = c(
              'text/csv',
              'text/comma-separated-values,text/plain',
              '.csv',
              '.txt'
            )
          )
        )
      )
      
    })
    
    #   ____________________________________________________________________________
    #   Buttons reactives                                                       ####
    
    
    shiny::observeEvent((input$deg_test_btn), {
      shiny::req(r$tcc)
      if (is.null(r$fit)) {
        r_dea$fit <- estimateDispersion(r$tcc)
        r$fit <- r_dea$fit
        
        if (golem::get_golem_options("server_version"))
          loggit::loggit(
            custom_log_lvl = TRUE,
            log_lvl = r$session_id,
            log_msg = "DEA"
          )
      }
      
      shiny::req(r$fit,
                 input$dea_fdr)
      
      if (is.null(input$reference) | is.null(input$perturbation)) {
        shinyalert::shinyalert("You must select at least one reference and one perturbation.",
                               type = "error")
      }
      shiny::req(input$reference, input$perturbation)
      
      if (any(input$reference %in% input$perturbation) |
          any(input$perturbation %in% input$reference)) {
        shinyalert::shinyalert("You tried to compare the same conditions!
                               You may need some coffee...",
                               type = "error")
      }
      shiny::req(!any(input$reference %in% input$perturbation) | !any(input$perturbation %in% input$reference))
      
      
      r_dea$tags <-
        estimateDEGs(r$fit,
                     reference = input$reference,
                     perturbation = input$perturbation
                     systematic_orientation = ifelse(is.null(input$multiple_DE_parameters_selection),
                                                     yes = FALSE,
                                                     no = input$multiple_DE_parameters_selection
                     )
      
      r_dea$top_tags <-
        r_dea$tags$table[r_dea$tags$table$FDR < input$dea_fdr,]
      r_dea$top_tags <-
        r_dea$top_tags[abs(r_dea$top_tags$logFC) > input$dea_lfc,]
      r_dea$DEGs <- r_dea$top_tags$genes
      
      if(length(input$reference)>1){
        r_dea$ref <- paste0("(", paste0(input$reference, collapse = " "), ")")
      } else {
        r_dea$ref <- input$reference
      }
      if(length(input$perturbation)>1){
        r_dea$trt <- paste0("(", paste0(input$perturbation, collapse = " "), ")")
      } else {
        r_dea$trt <- input$perturbation
      }
      
      r$DEGs[[paste(r_dea$ref, r_dea$trt)]] <- r_dea$DEGs
      r$top_tags[[paste(r_dea$ref, r_dea$trt)]] <- r_dea$top_tags
      r_dea$go <- NULL
      
      r_dea$lfc <- input$dea_lfc
      r_dea$fdr <- input$dea_fdr
      
      # --- Creating data for table display and download --- #
      
      top <- r_dea$top_tags
      top$Regulation <- ifelse(top$logFC > 0, "Up", "Down")
      
      columns <- c("logFC", "logCPM", "FDR", "Regulation")
      if (!is.null(r$gene_info)) {
        columns <- c(colnames(r$gene_info), columns)
        
        if (r$splicing_aware)
          ids <- get_locus(rownames(top), unique = FALSE)
        else
          ids <- rownames(top)
        top[, colnames(r$gene_info)] <-
          r$gene_info[match(ids, rownames(r$gene_info)), ]
      }
      
      r_dea$gene_table <- top[, columns]
      
    })
    
    
    #   ____________________________________________________________________________
    #   Summaries                                                               ####
    
    output$disp_estimate_summary <- shiny::renderUI({
      shiny::req(is.null(r$normalized_counts))
      
      numberColor = "red"
      number = "Normalisation needed"
      header = ""
      numberIcon = shiny::icon('times')
      
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
        numberIcon = shiny::icon('times')
      }
      else{
        if (is.null(r_dea$top_tags)) {
          numberColor = "orange"
          number = "Tests can be performed"
          header = ""
          numberIcon = shiny::icon('times')
        }
        else{
          numberColor = "olive"
          number = "Done"
          numberIcon =icon('check')
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
            numberIcon = shiny::icon('caret-up'),
            header = "up regulated",
            text = "genes",
            rightBorder = TRUE
          ),
          shinydashboardPlus::descriptionBlock(
            number = sum(r$top_tags[[paste(r_dea$ref, r_dea$trt)]]$logFC < 0),
            numberColor = "red",
            numberIcon = shiny::icon('caret-down'),
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
            label = "Download result table as .tsv",
            style = "material-flat",
            color = "success"
          )
        )),
        shiny::hr(),
        shinyWidgets::downloadBttn(
          ns("report"),
          "Generate html report",
          style = "material-flat",
          color = "default"
        )
      )
    })
    
    output$download_table_csv <- shiny::downloadHandler(
      filename = function() {
        paste(paste0("DEGs_", r_dea$ref, "-", r_dea$trt, ".tsv"))
      },
      content = function(file) {
        
        df <- r_dea$gene_table
        df$Gene_ID <- rownames(r_dea$gene_table)
        # if(stringr::str_detect(colnames(df), "label"))
        #   df$label <- stringr::str_replace(df$label, ';', '-')
        write.table(#df[, !stringr::str_detect(colnames(df), "description")]
          df,
          file = file,
          row.names = FALSE,
          sep = '\t',
          quote = FALSE
        )
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
        file.copy(
          system.file("extdata", "DEA_report.Rmd", package = "DIANE"),
          tempReport,
          overwrite = TRUE
        )
        file.copy(system.file("extdata", "favicon.ico", package = "DIANE"),
                  tempImage,
                  overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(r_dea = r_dea, r = r)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(
          tempReport,
          output_file = file,
          params = params,
          envir = new.env(parent = globalenv())
        )
      }
    )
    
    
    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
    ### download GO                                                             ####
    
    
    output$download_go_table <- shiny::downloadHandler(
      filename = function() {
        paste(
          paste0(
            "enriched_GOterms",
            input$reference,
            '_VS_',
            input$perturbation,
            '_',
            input$go_type,
            ".csv"
          )
        )
      },
      content = function(file) {
        write.csv(r_dea$go, file = file, quote = FALSE)
      }
    )
    
    
    #   ____________________________________________________________________________
    #   Result plots                                                            ####
    
    
    
    output$table_ui <- shiny::renderUI({
      if (is.null(r$normalized_counts)) {
        shinydashboardPlus::descriptionBlock(
          number = "Please normalize and filter raw data in normalization tab",
          numberColor = "orange",
          rightBorder = FALSE
        )
      }
      else
        DT::dataTableOutput(ns("deg_table"))
    })
    
    output$deg_table <- DT::renderDataTable({
      shiny::req(r$top_tags, r_dea$ref, r_dea$trt, r_dea$gene_table)
      shiny::req(r$top_tags[[paste(r_dea$ref, r_dea$trt)]])
      
      r_dea$gene_table
      
      DT::formatStyle(
        DT::datatable(r_dea$gene_table),
        columns = c("Regulation"),
        target = c("cell", "row"),
        backgroundColor = DT::styleEqual(c("Up", "Down"), c("#72F02466", c("#FF000035")))
      )
    })
    
    #   ____________________________________________________________________________
    #   Venn                                                                    ####
    
    ###TODO : changer le nom du fichier tÃ©lÃ©chargeable
    ###TODO : EmpÃªcher l'utilisateur de faire des intersections qui n'ont pas de sens
    ###TODO : rendre la partie intersection plus intuitive DONE
    ###TODO : utiliser de beaux boutons graphiques pour les up/down/all. DONE (mais pas super)
    ###TODO : Trouver autre chose que FALSE DONE
    ###TODO : echelle de l'image ! (trouvÃ©! changer simplement le "res"...)
    ###FIXME : Boutons du choix de la mÃ©thode de normalisation (awesomeRadio) qui foire sur la page normalisation... Si j'ajoute un bouton de mÃªme type quelque part ici Ã§a remarche. En regardant, il manque une propriÃ©tÃ© (un petit padding) si j'ai pas un awesomeRadio (mÃªme inutile) dans cette partie du programme. Je ne comprends pas.
    
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

    
    output$venn_lists_choice_2 <- shiny::renderUI({
      shiny::req(length(r$DEGs) > 1)
      shiny::fluidRow(
        ###All the buttons containing the list of genes.
        shiny::column(3,
                      shinyWidgets::pickerInput(
                        inputId = ns("venn_list_1"),
                        label = "Gene list 1",
                        choices = c("None" = FALSE, names(r$DEGs))
                      )),
        shiny::column(3,
                      shinyWidgets::pickerInput(
                        inputId = ns("venn_list_2"),
                        label = "Gene list 2",
                        choices = c("None" = FALSE, names(r$DEGs))
                      )),
        shiny::column(3,
                      shinyWidgets::pickerInput(
                        inputId = ns("venn_list_3"),
                        label = "Gene list 3",
                        choices = c("None" = FALSE, names(r$DEGs))
                      )),
        shiny::column(3,
                      shinyWidgets::pickerInput(
                        inputId = ns("venn_list_4"),
                        label = "Gene list 4",
                        choices = c("None" = FALSE, names(r$DEGs))
                      )), 
        ###All the buttons containg the "up / down" chocices.
        shiny::column(3,
                      shinyWidgets::checkboxGroupButtons(
                        inputId = ns("up_down_button_venn_1"),
                        label = "Gene subset",
                        choices = c("Up", "Down"),
                        selected = c("Up", "Down"),
                        justified = TRUE,
                        size = "sm",
                        checkIcon = list(yes = icon("ok",
                                                    lib = "glyphicon"))
                      )
        ),
        shiny::column(3,
                      shinyWidgets::checkboxGroupButtons(
                        inputId = ns("up_down_button_venn_2"),
                        label = "Gene subset",
                        choices = c("Up", "Down"),
                        selected = c("Up", "Down"),
                        justified = TRUE,
                        size = "sm",
                        checkIcon = list(yes = icon("ok",
                                                    lib = "glyphicon"))
                      )
        ),
        shiny::column(3,
                      shinyWidgets::checkboxGroupButtons(
                        inputId = ns("up_down_button_venn_3"),
                        label = "Gene subset",
                        choices = c("Up", "Down"),
                        selected = c("Up", "Down"),
                        justified = TRUE,
                        size = "sm",
                        checkIcon = list(yes = icon("ok",
                                                    lib = "glyphicon"))
                      )
        ),
        shiny::column(3,
                      shinyWidgets::checkboxGroupButtons(
                        inputId = ns("up_down_button_venn_4"),
                        label = "Gene subset",
                        choices = c("Up", "Down"),
                        selected = c("Up", "Down"),
                        justified = TRUE,
                        size = "sm",
                        checkIcon = list(yes = icon("ok",
                                                    lib = "glyphicon"))
                      )
        ),
      )
    })
    
    
    ###List of input gene list for venn diagram. Based on what user input.
    venn_list <- shiny::reactive({
      shiny::req(sum( ###Check that at least two list have a value != FALSE
        c(
          input$venn_list_1,
          input$venn_list_2,
          input$venn_list_3,
          input$venn_list_4
        ) != FALSE
      ) >= 2)
      venn_list <- list()
      
      for (comp in 1:4) {
        ###We test the 4 input DE list fields.
        if (!isFALSE(input[[paste0("venn_list_", comp)]])) {
          ###If the gene list is set to false, we just go to the next
          selected_comparison <-
            input[[paste0("venn_list_", comp)]] ###Extraction of the value.
          if(all(input[[paste0("up_down_button_venn_", comp)]] == "")){
            venn_list[[selected_comparison]] <-
              r$top_tags[[selected_comparison]]$genes
          } else if (all(input[[paste0("up_down_button_venn_", comp)]] == "Up")) {
            #Only up is selected
            venn_list[[paste0(selected_comparison, " up")]] <-
              r$top_tags[[selected_comparison]][r$top_tags[[selected_comparison]]$logFC > 0 , "genes"]
          } else if (all(input[[paste0("up_down_button_venn_", comp)]] == "Down")) {
            #only down is selected
            venn_list[[paste0(selected_comparison, " down")]] <-
              r$top_tags[[selected_comparison]][r$top_tags[[selected_comparison]]$logFC < 0 , "genes"]
          } else {
            #Up and down are selected.
            venn_list[[selected_comparison]] <-
              r$top_tags[[selected_comparison]]$genes
          }
        }
      }
      # print(head(venn_list))
      venn_list
    })
    
    
    ###Venn diagram plot. The res parameter as a direct impact on text size.
    output$venn <- shiny::renderPlot({
      shiny::req(venn_list)
      validate(
        need(length(names(venn_list())) > 1, "Please specify between two and four genes list.")
      )
      draw_venn(venn_list())
    }, res = 100)
    
    
    output$venn_spec_comp_choice <- shiny::renderUI({
      shiny::req(venn_list())
      tagList(
        shiny::selectInput(
          ns("venn_spec_comp"),
          label = "Genes specific to a comparison :",
          choices = names(venn_list())[!FALSE],
        )
      )
    })
    
    ###Part with download intersection.
    output$venn_spec_comp_choice_2 <- shiny::renderUI({
      shiny::req(venn_list())
      shiny::req(length(names(venn_list())) > 1)
      tagList(
        shiny::h3("Download subsets of genes"),
        tags$table(style = "width: 100%; text-align: center;",
                   tags$tr(
                     tags$td(
                       style = "width: 30%; ",
                       shiny::selectInput(
                         ns("venn_genes_intersection"),
                         label = "Genes present in the intersection of :",
                         choices = names(venn_list())[!FALSE],
                         multiple = TRUE
                       )
                     ),
                     tags$td(style = "width: 40%; padding: 0 5px 0 5px;",
                             h4(
                               " Which are also absent from the following lists "
                             )),
                     tags$td(
                       style = "width: 30%",
                       shiny::selectInput(
                         ns("venn_genes_union_absent"),
                         label = "Genes absent in lists :",
                         choices = names(venn_list())[!FALSE],
                         multiple = TRUE
                       )
                     ),
                   )),
      )
    })
    
    output$download_specific_venn_2 <- shiny::downloadHandler(
      filename = function() {
        if (!is.null(input$venn_genes_union_absent)) {
          ###Nom pas encore trÃ¨s sexy :'(
          stringr::str_replace_all(paste(
            paste0(
              "Venn_specific_to ",
              paste0(input$venn_genes_intersection, collapse = "_and_")
              ,
              "_NOT_",
              paste0(input$venn_genes_union_absent, collapse = "_and_"),
              ".csv"
            ),
            collapse = "_"
          ),
          pattern = " ",
          replacement = "_")
        } else {
          stringr::str_replace_all(paste(paste0(
            "Venn_specific_to ",
            paste0(input$venn_genes_intersection, collapse = "_and_")
          )),
          pattern = " ",
          replacement = "_")
        }
      },
      
      content = function(file) {
        #if (!any(input$venn_genes_intersection == input$venn_genes_union_absent)) {
        write.table(
          ###Intersection of the list on the left - union of the list of the rigth.
          setdiff(Reduce(intersect, venn_list()[input$venn_genes_intersection]),
                  Reduce(union, venn_list()[input$venn_genes_union_absent])),
          file = file,
          row.names = FALSE,
          sep = ';',
          quote = FALSE,
          col.names = FALSE
        )
        #} #else {
        #shinyalert::shinyalert( ###Should just disable the button instead of this, here.
        #  "You cannot select identical conditions in the gene,
        #                       list to include and in the gene list to remove",
        #  type = "error"
        #)
        # }
      }
    )
    
    output$venn_spec_comp_bttn_2 <- shiny::renderUI({
      shiny::req(venn_list())
      shiny::req(length(names(venn_list())) > 1)
      shiny::validate(
        shiny::need(
          !any(
            input$venn_genes_intersection %in% input$venn_genes_union_absent
          ),
          "You cannot select identical conditions in the gene list to include and in the gene list to remove"
        ),
        shiny::need(
          input$venn_genes_intersection != "",
          "You must select at least one gene list to include."
        )
      )
      tagList(
        shinyWidgets::downloadBttn(
          outputId = ns("download_specific_venn_2"),
          label = paste(
            "Download genes specific to the",
            paste0(input$venn_genes_intersection, collapse = "/"),
            "list"
          ),
          style = "material-flat",
          color = "success"
        )
      )
    })
    
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
        if (is.null(r$custom_go)) {
          if (!is.null(input$go_data)) {
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
            shinyalert::shinyalert(
              "Please input Gene to GO term file. ",
              "Only some main model organisms are supported,
                                   but you can input your own gene - GO terms matching.",
              type = "error"
            )
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
        if(input$up_down_go_radio == "All")
          genes <- r_dea$top_tags$genes
        if(input$up_down_go_radio == "Up-regulated")
          genes <- r_dea$top_tags[r_dea$top_tags$logFC > 0,]$genes
        if(input$up_down_go_radio == "Down-regulated")
          genes <- r_dea$top_tags[r_dea$top_tags$logFC < 0,]$genes
        
        universe <-
          intersect(rownames(r$normalized_counts), GOs[, 1])
        
        
        if (length(universe) == 0) {
          r$custom_go <- NULL
          shinyalert::shinyalert(
            "Invalid first column",
            "The first column did not match any gene ID from
              differential expression analysis",
            type = "error"
          )
        }
        
        shiny::req(length(universe) > length(genes))
        r_dea$go <-
          enrich_go_custom(genes, universe, GOs, GO_type = input$go_type)
        
        
        ################# known organisms
        
      } else{
        
        if(input$up_down_go_radio == "All")
          genes <- r_dea$top_tags$genes
        if(input$up_down_go_radio == "Up-regulated")
          genes <- r_dea$top_tags[r_dea$top_tags$logFC > 0,]$genes
        if(input$up_down_go_radio == "Down-regulated")
          genes <- r_dea$top_tags[r_dea$top_tags$logFC < 0,]$genes
        
        background <- rownames(r$normalized_counts)
        
        if (r$splicing_aware) {
          genes <- get_locus(genes)
          background <- get_locus(background)
        }
        
        if (r$organism == "Lupinus albus") {
          GOs <- DIANE:::lupine$go_list
          universe <- intersect(background, GOs[, 1])
          r_dea$go <- enrich_go_custom(genes, universe, GOs,
                                       GO_type = input$go_type)
        }
        else if (stringr::str_detect(r$organism, "Oryza")) {
          data("go_matchings", package = "DIANE")
          
          GOs <- go_matchings[[r$organism]]
          universe <- intersect(background, GOs[, 1])
          r_dea$go <- enrich_go_custom(genes, universe, GOs,
                                       GO_type = input$go_type)
        }
        else{
          if (r$organism == "Arabidopsis thaliana") {
            genes <- convert_from_agi(genes)
            background <- convert_from_agi(background)
            org = org.At.tair.db::org.At.tair.db
          }
          
          if (r$organism == "Homo sapiens") {
            genes <- convert_from_ensembl(genes)
            background <- convert_from_ensembl(background)
            org = org.Hs.eg.db::org.Hs.eg.db
          }
          
          if (r$organism == "Mus musculus") {
            genes <- convert_from_ensembl_mus(genes)
            background <- convert_from_ensembl_mus(background)
            org = org.Mm.eg.db::org.Mm.eg.db
          }
          
          if (r$organism == "Drosophilia melanogaster") {
            genes <- convert_from_ensembl_dm(genes)
            background <- convert_from_ensembl_dm(background)
            org = org.Dm.eg.db::org.Dm.eg.db
          }
          
          if (r$organism == "Caenorhabditis elegans") {
            genes <- convert_from_ensembl_ce(genes)
            background <- convert_from_ensembl_ce(background)
            org = org.Ce.eg.db::org.Ce.eg.db
          }
          
          if (r$organism == "Escherichia coli") {
            genes <- convert_from_ensembl_eck12(genes)
            background <- convert_from_ensembl_eck12(background)
            org = org.EcK12.eg.db::org.EcK12.eg.db
          }
          
          # TODO add check if it is entrez with regular expression here
          shiny::req(length(genes) > 0, length(background) > 0)
          r_dea$go <-
            enrich_go(genes,
                      background,
                      org = org,
                      GO_type = input$go_type)
        }
        
      }
      if (golem::get_golem_options("server_version"))
        loggit::loggit(
          custom_log_lvl = TRUE,
          log_lvl = r$session_id,
          log_msg = "GO enrichment DEA"
        )
      
    })
    
    #   ____________________________________________________________________________
    #   go results                                                              ####
    
    
    output$go_table <- DT::renderDataTable({
      shiny::req(r_dea$go)
      r_dea$go[, c("Description", "GeneRatio", "BgRatio", "p.adjust")]
    })
    
    output$max_go_choice <- shiny::renderUI({
      shiny::req(r_dea$go)
      shiny::req(input$draw_go == "Dot plot")
      shiny::numericInput(
        ns("n_go_terms"),
        label = "Top number of GO terms to plot :",
        min = 1,
        value = dim(r_dea$go)[1]
      )
    })
    
    output$go_plot <- plotly::renderPlotly({
      shiny::req(r_dea$go)
      max = ifelse(is.na(input$n_go_terms),
                   dim(r_dea$go)[1],
                   input$n_go_terms)
      draw_enrich_go(r_dea$go, max_go = max)
    })
    
    output$go_map_plot <- shiny::renderPlot({
      shiny::req(r_dea$go)
      draw_enrich_go_map(r_dea$go)
    })
    
    output$go_results <- shiny::renderUI({
      shiny::req(r_dea$go)
      
      if (nrow(r_dea$go) == 0) {
        shinyalert::shinyalert(
          "No enriched GO terms were found",
          "It can happen if input gene list is not big enough",
          type = "error"
        )
      }
      
      shiny::req(nrow(r_dea$go) > 0)
      
      if (input$draw_go == "Data table") {
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
        if (input$draw_go == "Enrichment map") {
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
