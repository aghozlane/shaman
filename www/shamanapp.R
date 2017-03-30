#!/usr/bin/Rscript --vanilla
library(shiny)
system("Rscript -e 'library(\"shiny\");runGitHub(\"pierreLec/KronaRShy\", port=5438)'",wait=FALSE)
runGitHub("aghozlane/shaman")
