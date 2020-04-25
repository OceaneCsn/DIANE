# No Remotes ----
# Attachments ----
to_install <- c("attempt", "config", "DT", "ggplot2", "glue", "golem", "htmltools", "plotly", "shinipsum", "shiny", "shinydashboard", "shinydashboardPlus", "shinythemes", "shinyWidgets", "stringr")
  for (i in to_install) {
    message(paste("looking for ", i))
    if (!requireNamespace(i)) {
      message(paste("     installing", i))
      install.packages(i)
    }
  }
