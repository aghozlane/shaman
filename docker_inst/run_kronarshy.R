source("/srv/shiny-server/renv/activate.R")
shiny::runApp("/srv/shiny-server/kronarshy/", port=5438, host="0.0.0.0")