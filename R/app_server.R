#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
 
  
#   ____________________________________________________________________________
#   reactive values                                                         ####

  r <- shiny::reactiveValues(
    raw_counts = NULL,
    normalized_counts = NULL,
    normalized_counts_pre_filter = NULL,
    tcc = NULL,
    conditions = NULL,
    design = NULL,
    DEGs = list(),
    top_tags = list(),
    clusterings = list(),
    current_comparison = NULL
  )
  

#   ____________________________________________________________________________
#   Server modules                                                          ####

  shiny::callModule(mod_context_server, "context_ui_1")
  shiny::callModule(mod_import_data_server, "import_data_ui_1", r)
  shiny::callModule(mod_normalisation_server, "normalisation_ui_1", r)
  shiny::callModule(mod_differential_expression_analysis_server, 
                    "differential_expression_analysis_ui_1", r)
  
  # clustering modules
  shiny::callModule(mod_clustering_server, "clustering_ui_1", r)
  shiny::callModule(mod_cluster_exploration_server, 
                    "cluster_exploration_ui_1", r)
  
  shiny::callModule(mod_network_inference_server, "network_inference_ui_1")
  shiny::callModule(mod_network_analysis_server, "network_analysis_ui_1")
  shiny::callModule(mod_module_analysis_server, "module_analysis_ui_1")
  
}
