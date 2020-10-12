# Building a Prod-Ready, Robust Shiny Application.
# 
# README: each step of the dev files is optional, and you don't have to 
# fill every dev scripts before getting started. 
# 01_start.R should be filled at start. 
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
# 
# 
###################################
#### CURRENT FILE: DEV SCRIPT #####
###################################

# Engineering

## Dependencies ----
## Add one line by package you want to add as dependency
#usethis::use_package( "thinkr" )
usethis::use_package( "stringr" )
usethis::use_package( "ggplot2" )
usethis::use_package( "plotly" )
usethis::use_package( "shinydashboard" )
usethis::use_package( "shinythemes" )
usethis::use_package( "shinydashboardPlus", min_version = "0.7.5" )
usethis::use_package( "shinyWidgets" )
usethis::use_package( "stats" )
usethis::use_package( "pheatmap" )
usethis::use_package( "HTSFilter" )
usethis::use_package( "TCC" )
usethis::use_package( "shinybusy" )
usethis::use_package( "reshape2" )
usethis::use_package( "limma" )
usethis::use_package( "markdown" )
usethis::use_package( "edgeR" )
usethis::use_package( "dashboardthemes" )
usethis::use_package( "shinyalert" )
usethis::use_package( "RColorBrewer" )
usethis::use_package( "coseq" )
usethis::use_package( "utils" )
usethis::use_package( "shinyjs" )
usethis::use_package( "parallel" )
usethis::use_package( "GENIE3" )
usethis::use_package( "igraph" )
usethis::use_package( "visNetwork" )
usethis::use_package( "gridExtra" )
usethis::use_package( "org.At.tair.db" )
usethis::use_package( "clusterProfiler" )
usethis::use_package( "plotly" )
usethis::use_package( "org.Hs.eg.db" ) 
usethis::use_package( "foreach" ) 
usethis::use_package( "doParallel" ) 
usethis::use_package( "doRNG" ) 
usethis::use_package( "ade4" ) 
usethis::use_package( "ggraph" ) 
usethis::use_package( "GO.db" ) 
usethis::use_package( "ggridges" ) 
usethis::use_package( "randomForest" ) 
usethis::use_package( "rfPermute" ) 
usethis::use_package( "knitr" ) 
usethis::use_package( "rmarkdown" ) 
usethis::use_package( "ggVennDiagram" ) 
usethis::use_package( "loggit" ) 

#usethis::use_package( "MASS" )


## Add modules ----
## Create a module infrastructure in R/
golem::add_module( name = "import_data" ) # Name of the module
golem::add_module( name = "differential_expression_analysis" ) # Name of the module
golem::add_module( name = "context" )
golem::add_module( name = "normalisation" )
golem::add_module( name = "clustering" )
golem::add_module( name = "cluster_exploration" )
golem::add_module( name = "network_inference" )
golem::add_module( name = "network_analysis" )
golem::add_module( name = "module_analysis" )
golem::add_module( name = "module_levels" )
golem::add_module( name = "datasets" )


## Add helper functions ----
## Creates ftc_* and utils_*
golem::add_fct( "heatmap" )
golem::add_fct( "normalisation")
golem::add_fct( "dea")
golem::add_fct( "coseq")
golem::add_fct( "glm")
golem::add_fct( "network_inference")
golem::add_fct( "network_analysis")
golem::add_fct( "go")
golem::add_fct( "regressors_grouping")
golem::add_fct( "GENIE3_oob")
golem::add_fct( "edge_testing")
golem::add_fct( "geneIDs")


#golem::add_utils( "helpers" )

## External resources
## Creates .js and .css files at inst/app/www
#golem::add_js_file( "script" )
#golem::add_js_handler( "handlers" )
#golem::add_css_file( "custom" )

## Add internal datasets ----
## If you have data in your package
usethis::use_data(abiotic_stresses, version = 3, overwrite = T)
usethis::use_data(gene_annotations , version = 3, overwrite = T)
usethis::use_data(regulators_per_organism , version = 3, overwrite = T)
usethis::use_data(lupine, version = 3, overwrite = T, internal = TRUE)


#usethis::use_data_raw( name = "raw_data_demo", open = FALSE ) 

## Tests ----
## Add one line by test you want to create
#usethis::use_test( "app" )

# Documentation

golem::browser_button()

## Vignette ----
usethis::use_vignette("DIANE_Programming_Interface")
usethis::use_vignette("DIANE")

devtools::build_vignettes()




usethis::use_pkgdown()

pkgdown::build_site()


usethis::use_github_action("pkgdown")
usethis::use_github_action("test-coverage")
usethis::use_github_action_check_standard()


usethis::use_travis()
usethis::use_travis_badge()


## Code coverage ----
## (You'll need GitHub there)

#usethis::use_github()
#usethis::use_travis()
#usethis::use_appveyor()

# You're now set! ----
# go to dev/03_deploy.R

golem::add_shinyserver_file()

rstudioapi::navigateToFile("dev/03_deploy.R")

