# library(shinydashboard)
# library(shinythemes)
# try(library(shinydashboardPlus), silent = T)


logo_diane <- dashboardthemes::shinyDashboardLogoDIY(
  boldText = ""
  ,
  mainText = ""
  ,
  textSize = 16
  ,
  badgeText = "DIANE"
  ,
  badgeTextColor = "white"
  ,
  badgeTextSize = 7
  ,
  badgeBackColor = "#5FBF64"
  ,
  badgeBorderRadius = 5
  
)


#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
#'
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # List the first level UI elements here
    shinydashboard::dashboardPage(
      shinydashboard::dashboardHeader(title = logo_diane),
      
      
#   ____________________________________________________________________________
#   sidebar                                                                 ####

      shinydashboard::dashboardSidebar(
        width = "300px",
        shinydashboard::sidebarMenu(
          shinydashboard::menuItem(
            "Context",
            tabName = "context_tab",
            icon = shiny::icon("seedling")
          ),
          shinydashboard::menuItem(
            "Data import",
            tabName = "data_import_tab",
            icon = shiny::icon("table")
          ),
          shinydashboard::menuItem(
            "Normalisation",
            tabName = "normalisation_tab",
            icon = shiny::icon("table")
          ),
          shinydashboard::menuItem(
            "Gene expression levels",
            tabName = "levels_tab",
            icon = shiny::icon("chart-line")
          ),
          shinydashboard::menuItem(
            "Differential Expression",
            tabName = "dea_tab",
            icon = shiny::icon("greater-than-equal")
          ),
          shinydashboard::menuItem(
            "Expression based clustering",
            startExpanded = TRUE,
            icon = shiny::icon("circle-notch"),
            shinydashboard::menuSubItem(tabName = "clustering_tab",
                                        text = "Run a clustering"),
            shinydashboard::menuSubItem(tabName = "cluster_exploration_sub_tab",
                                        text = "Explore clusters")
          ),
          
          shinydashboard::menuItem(
            "Gene Regulatory Network",
            startExpanded = TRUE,
            icon = shiny::icon("project-diagram"),
            shinydashboard::menuSubItem(tabName = "network_inference_tab",
                                        text = "Network inference"),
            shinydashboard::menuSubItem(tabName = "network_analysis_tab",
                                        text = "Network analysis")
            )
          )
      ),


#   ____________________________________________________________________________
#   body                                                                    ####

      
      shinydashboard::dashboardBody(
        dashboardthemes::shinyDashboardThemes(theme = "grey_light"),
        
        
        shinydashboard::tabItems(
          shinydashboard::tabItem(tabName = "context_tab",
                                  mod_context_ui("context_ui_1")),
          shinydashboard::tabItem(tabName = "data_import_tab",
                                  mod_import_data_ui("import_data_ui_1")),
          shinydashboard::tabItem(tabName = "normalisation_tab",
                                  mod_normalisation_ui("normalisation_ui_1")),
          shinydashboard::tabItem(
            tabName = "levels_tab",
            mod_module_levels_ui("module_levels_ui_1")
          ),
          shinydashboard::tabItem(
            tabName = "dea_tab",
            mod_differential_expression_analysis_ui("differential_expression_analysis_ui_1")
          ),
          
          
          
          
#   ____________________________________________________________________________
#   clustering                                                              ####

          shinydashboard::tabItem(tabName = "clustering_tab",
                                  mod_clustering_ui("clustering_ui_1")),
          shinydashboard::tabItem(
            tabName = "cluster_exploration_sub_tab",
            mod_cluster_exploration_ui("cluster_exploration_ui_1")
          ),
          
         
          
#   ____________________________________________________________________________
#   network                                                                 ####

          shinydashboard::tabItem(tabName = "network_inference_tab",
                                  mod_network_inference_ui("network_inference_ui_1")),
          shinydashboard::tabItem(tabName = "network_analysis_tab",
                                  mod_network_analysis_ui("network_analysis_ui_1")
          )
          # shinydashboard::tabItem(tabName = "module_analysis_tab",
          #                         mod_module_analysis_ui("module_analysis_ui_1")
          # )
         
        )
        
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  golem::add_resource_path('www', app_sys('app/www'))
  #add_resource_path('md', app_sys('app/www/md'))
  
  tags$head(golem::favicon(),
            golem::bundle_resources(path = app_sys('app/www'),
                             app_title = 'DIANE'),
            # Add here other external resources
            #shinyalert::useShinyalert())
  )
}
