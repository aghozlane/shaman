#!/usr/bin/Rscript --vanilla
system("Rscript -e 'library(\"shiny\");runApp(appDir= KronaRShy\", port=5438)'",wait=FALSE)
runGitHub("aghozlane/shaman")
