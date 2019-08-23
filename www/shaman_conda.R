#!/usr/bin/Rscript --vanilla
prefix=Sys.getenv(x = "PREFIX")
prefix=paste(prefix, "bin", sep=.Platform$file.sep)
system(paste("Rscript -e 'library(\"shiny\");runApp(appDir=", prefix,"KronaRShy\", port=5438)'", sep=.Platform$file.sep, wait=FALSE)
runApp(paste("appDir= ", prefix, "shaman", sep=.Platform$file.sep)
