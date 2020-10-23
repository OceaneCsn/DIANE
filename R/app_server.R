#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
#' 
options(shiny.maxRequestSize=30*1024^2)
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
    current_comparison = NULL,
    current_network = NULL,
    regulators = NULL,
    use_demo = NULL,
    networks = list(),
    splicing_aware = NULL,
    gene_info = NULL,
    organism = NULL,
    custom_go = NULL,
    session_id = as.character(floor(runif(1)*1e20))
  )
  
  
  #   ____________________________________________________________________________
  #   logs                                                                    ####
  
  LOG_FILE = "./logs/loggit.log"
  SESSION_ID_FILE = "./logs/next_id.txt"
  
  
  if(golem::get_golem_options("server_version")){
    if (!file.exists(SESSION_ID_FILE)){
      file.create(SESSION_ID_FILE)
      write(1, SESSION_ID_FILE )
    }
  
  
    session_id <- readLines(SESSION_ID_FILE)
    close( file( SESSION_ID_FILE, open="w" ) )
    write(as.numeric(session_id) + 1, SESSION_ID_FILE )
    
    loggit::set_logfile(LOG_FILE)
    loggit::set_timestamp_format("%Y-%m-%d %H:%M:%S")
    
    
    loggit::loggit(custom_log_lvl = TRUE,
                   log_lvl = session_id,
                   log_msg = "connection")
    
    r$session_id <- session_id
  }
  
  #session_id <- as.character(floor(runif(1)*1e20))
  
  

#   ____________________________________________________________________________
#   Server modules                                                          ####

  shiny::callModule(mod_context_server, "context_ui_1")
  shiny::callModule(mod_import_data_server, "import_data_ui_1", r)
  shiny::callModule(mod_normalisation_server, "normalisation_ui_1", r)
  
  shiny::callModule(mod_module_levels_server, "module_levels_ui_1", r)
  shiny::callModule(mod_differential_expression_analysis_server, 
                    "differential_expression_analysis_ui_1", r)
  
  # clustering modules
  shiny::callModule(mod_clustering_server, "clustering_ui_1", r)
  shiny::callModule(mod_cluster_exploration_server, 
                    "cluster_exploration_ui_1", r)
  
  shiny::callModule(mod_network_inference_server, "network_inference_ui_1", r)
  shiny::callModule(mod_network_analysis_server, "network_analysis_ui_1", r)

  shiny::callModule(mod_datasets_server, "datasets_ui_1")
  
}
