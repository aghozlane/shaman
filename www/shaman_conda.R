#!/usr/bin/env Rscript --vanilla
prefix=Sys.getenv(x = "CONDA_PREFIX")
prefix=paste(prefix, "bin", sep=.Platform$file.sep)
system(paste("python ", prefix, "shaman_bioblend/shaman_bioblend.py -w ", prefix,"shamanapp/www/masque/ -s -d", sep=.Platform$file.sep))
workdir = getwd()
kronarshy=tempfile("kronarshy",fileext=".R", tmpdir =workdir )
write(paste("library(\"shiny\");runApp(appDir=\"", prefix,"KronaRShy\", port=5438)", sep=.Platform$file.sep),kronarshy)
system(paste("Rscript", kronarshy), wait=FALSE)
shamanapp=tempfile("shaman", fileext=".R", tmpdir =workdir)
write(paste("library(\"shiny\");runApp(appDir=\"", prefix,"shamanapp\", port=80, host=\"0.0.0.0\")", sep=.Platform$file.sep), shamanapp)
system(paste("Rscript", shamanapp), wait=FALSE)
