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
            offLabel = "Vulcano"
          ),
          
          
          shiny::plotOutput(ns("ma_vulcano"), height = "700px")
          
        ),
        shiny::tabPanel(
          title = "Heatmap",
          shiny::uiOutput(ns("heatmap_conditions_choice")),
          shiny::plotOutput(ns("heatmap"), height = "700px")
        )
        #shiny::tabPanel(title = "EdgeR summary",
        #shiny::verbatimTextOutput(ns("edgeR_summary")),
        
      )
      
    ),
    shiny::actionButton(ns("browser"), "backdoor"),
  )
}
# TODO place spinners correctly

# TODO cleaner demo data

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
      trt = NULL
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
    #   Buttons reactives                                                       ####
    
    
    shiny::observeEvent((input$estimate_disp_btn), {
      shiny::req(r$tcc)
      #r_dea$fit <- estimateDispersion(r$tcc)
      r$fit <- estimateDispersion(r$tcc)
    })
    
    shiny::observeEvent((input$deg_test_btn), {
      shiny::req(r$fit,
                 input$dea_fdr,
                 input$reference,
                 input$perturbation)
      
      print(input$reference == input$perturbation)
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
      
    })
    
    # output$edgeR_summary <- shiny::renderPrint({
    #   shiny::req(r_dea$fit)
    #   print(r_dea$fit$coefficients)
    # })
    #
    #
    #   ____________________________________________________________________________
    #   Summaries                                                               ####
    
    output$disp_estimate_summary <- shiny::renderUI({
      if (is.null(r$fit)) {
        number_color = "orange"
        number = "Needed"
        header = ""
        number_icon = "fa fa-times"
      }
      else{
        number_color = "olive"
        number = "Done"
        number_icon = "fa fa-check"
        header = "See EdgeR summary tab for more details"
      }
      shinydashboardPlus::descriptionBlock(
        number = number,
        number_color = number_color,
        number_icon = number_icon,
        header = header,
        right_border = FALSE
      )
    })
    
    
    output$deg_test_summary <- shiny::renderUI({
      if (is.null(r$fit)) {
        number_color = "red"
        number = "Dispersion estimation needed"
        header = ""
        number_icon = "fa fa-times"
      }
      else{
        if (is.null(r_dea$top_tags)) {
          number_color = "orange"
          number = "Test can be performed"
          header = ""
          number_icon = "fa fa-times"
        }
        else{
          number_color = "olive"
          number = "Done"
          number_icon = "fa fa-check"
          header = "See plots and tables for more details"
        }
      }
      shinydashboardPlus::descriptionBlock(
        number = number,
        number_color = number_color,
        number_icon = number_icon,
        header = header,
        right_border = FALSE
      )
    })
    
    output$deg_number_summary <- shiny::renderUI({
      shiny::req(r$top_tags, r_dea$ref, r_dea$trt)
      shiny::req(r$top_tags[[paste(r_dea$ref, r_dea$trt)]])
      
      tagList(
        shiny::fluidRow(
          shinydashboardPlus::descriptionBlock(
            number = sum(r$top_tags[[paste(r_dea$ref, r_dea$trt)]]$logFC > 0),
            number_color = "olive",
            number_icon = "fa fa-caret-up",
            header = "up regulated",
            text = "genes",
            right_border = TRUE
          ),
          shinydashboardPlus::descriptionBlock(
            number = sum(r$top_tags[[paste(r_dea$ref, r_dea$trt)]]$logFC < 0),
            number_color = "red",
            number_icon = "fa fa-caret-down",
            header = "down-regulated",
            text = "genes",
            right_border = FALSE
          )
        )
      )
    })
    
    
    #   ____________________________________________________________________________
    #   Dl button                                                               ####
    
    output$dl_bttns <- shiny::renderUI({
      shiny::req(r$top_tags, r_dea$ref, r_dea$trt)
      shiny::req(r$top_tags[[paste(r_dea$ref, r_dea$trt)]])
      shiny::fluidRow(
        shinyWidgets::downloadBttn(
          outputId = ns("download_table_csv"),
          label = "Download result table as .csv",
          style = "bordered",
          color = "default"
        )
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
    #   Result plots                                                            ####
    
    output$deg_table <- DT::renderDataTable({
      shiny::req(r$top_tags, r_dea$ref, r_dea$trt)
      shiny::req(r$top_tags[[paste(r_dea$ref, r_dea$trt)]])
      
      top <- r_dea$top_tags
      top$Regulation <- ifelse(top$logFC > 0, "Up", "Down")
      DT::formatStyle(
        DT::datatable(top[, c("logFC", "logCPM", "FDR", "Regulation")],
                      options(list(
                        pageLength = 20,
                        lengthMenu = c(10, 20)
                      ))),
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
    
    output$ma_vulcano <- shiny::renderPlot({
      shiny::req(r$top_tags)
      shiny::req(r$top_tags[[paste(r_dea$ref, r_dea$trt)]])
      draw_DEGs(
        tags = r_dea$tags,
        fdr = input$dea_fdr,
        lfc = input$dea_lfc,
        MA = input$MA_vulcano_switch
      )
    })
    
    shiny::observeEvent(input$browser, {
      browser()
    })
    
    
  }

## To be copied in the UI
# mod_differential_expression_analysis_ui("differential_expression_analysis_ui_1")

## To be copied in the server
# callModule(mod_differential_expression_analysis_server, "differential_expression_analysis_ui_1")
