# Set options here
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode

# Detach all loaded packages and clean your environment
#golem::detach_all_attached()
rm(list = ls(all.names = TRUE))

# Document and reload your package
golem::document_and_reload()

# Run the application
run_app()


# NOTE : fix point d'interrogation (coseq, poisson, ew!)

# TODO : ajouter des ? pour acp, et poisson regression

# Dire dans la vignette qu'en fait il faut pas normaliser les données de démo