library(shinydashboard)
library(shinythemes)
library(shinipsum)
library(shinydashboardPlus)
#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # List the first level UI elements here
    dashboardPage(
      skin = "black",
      
      dashboardHeader(title = "DIANE"),
      dashboardSidebar(
        sidebarMenu(
          # menuItem(
          #   "Context",
          #   tabName = "context_tab",
          #   icon = icon("seedling")
          # ),
          menuItem(
            "Data import",
            tabName = "data_import_tab",
            icon = icon("table")
          ),
          menuItem(
            "Differential Expression Analysis",
            tabName = "DEA_tab",
            icon = icon("table")
          ),
          menuItem(
            "Expression based clustering",
            tabName = "coseq",
            icon = icon("greater-than-equal")
          ),
          menuItem(
            "Network inference",
            tabName = "network_tab",
            icon = icon("circle-notch")
          )
        )
      ),
      
      dashboardBody(

        tabItems(
          tabItem( tabName = "context_tab", mod_context_ui("context_ui_1")),
          tabItem( tabName = "data_import_tab", mod_import_data_ui("import_data_ui_1"))

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
  
  tags$head(favicon(),
            bundle_resources(path = app_sys('app/www'),
                             app_title = 'DIANE'),
            # Add here other external resources
            # for example, you can add shinyalert::useShinyalert() )
            #shinyalert::useShinyalert()
  )
}
