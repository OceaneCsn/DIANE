library(shinydashboard)
library(shinythemes)
try(library(shinydashboardPlus), silent = T)
#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # List the first level UI elements here
    shinydashboard::dashboardPage(
      skin = "black", 
      
      shinydashboard::dashboardHeader(title = "DIANE"),
      shinydashboard::dashboardSidebar(
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
            "Differential Expression Analysis",
            tabName = "dea_tab",
            icon = shiny::icon("table")
          ),
          shinydashboard::menuItem(
            "Expression based clustering",
            tabName = "coseq",
            icon = shiny::icon("greater-than-equal")
          ),
          shinydashboard::menuItem(
            "Network inference",
            tabName = "network_tab",
            icon = shiny::icon("circle-notch")
          )
        )
      ),
      
      shinydashboard::dashboardBody(

        shinydashboard::tabItems(
          shinydashboard::tabItem( tabName = "context_tab", mod_context_ui("context_ui_1")),
          shinydashboard::tabItem( tabName = "data_import_tab", mod_import_data_ui("import_data_ui_1")),
          shinydashboard::tabItem( tabName = "normalisation_tab", mod_normalisation_ui("normalisation_ui_1")),
          shinydashboard::tabItem( tabName = "dea_tab", mod_differential_expression_analysis_ui("differential_expression_analysis_ui_1"))

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
  add_resource_path('www', app_sys('app/www'))
  #add_resource_path('md', app_sys('app/www/md'))
  
  tags$head(favicon(),
            bundle_resources(path = app_sys('app/www'),
                             app_title = 'DIANE'),
            # Add here other external resources
            #shinyalert::useShinyalert()
  )
}
