library(packrat)
packrat::init("/opt/packman/shaman/")
library(shiny)
runApp("/srv/shiny-server/kronarshy/", port=5438, host="0.0.0.0")