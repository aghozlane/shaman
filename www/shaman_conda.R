#!/usr/bin/env Rscript --vanilla
prefix=Sys.getenv(x = "PREFIX")
prefix=paste(prefix, "bin", sep=.Platform$file.sep)
system(paste("python ", prefix, "shaman_bioblend/shaman_bioblend.py -w ", prefix,"shamanapp/www/masque/ -s -d", sep=.Platform$file.sep))
system(paste(prefix, "Rscript -e 'library(\"shiny\");runApp(appDir=", prefix,"KronaRShy\", port=5438)'", sep=.Platform$file.sep, wait=FALSE)
runApp(paste("appDir= ", prefix, "shamanapp", sep=.Platform$file.sep)
