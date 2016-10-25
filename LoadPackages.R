if (!require("Rcpp")){
  install.packages("Rcpp")
}
if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}
if(!require(shinydashboard)){
  install.packages('shinydashboard')
  library(shinydashboard)
}

if(!require(rjson)){
  install.packages('rjson')
}

if(!require(devtools)){
  install.packages('devtools')
}

if (!require(psych)) {
  install.packages('psych')
  library(psych)
}

if (!require(ggplot2)) {
  install.packages('ggplot2')
  library(ggplot2)
}

if (!require(vegan)) {
  install.packages('vegan')
  library(vegan)
}

if (!require(dendextend)) {
  install.packages('dendextend')
  library(dendextend)
}

if (!require(circlize)) {
  install.packages('circlize')
  library(circlize)
}


if (!require(devtools)) {
  install.packages('devtools')
  library(devtools)
}
# if(!require(treeWeightD3)){
#   devtools::install_git('https://gitlab.pasteur.fr/plechat/treeWeightD3')
#   library(treeWeightD3)
# }
if (!require(BiocInstaller)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("BiocInstaller")
  library(BiocInstaller)
}

if (!require(d3heatmap)) {
  devtools::install_github('aghozlane/d3heatmap')
  library(d3heatmap)
}

# Let us use biomformat instead of biom
if (!require(biomformat)){
 devtools::install_github("aghozlane/biomformat")
library(biomformat)
}

if (!require(scatterD3)) {
  devtools::install_github('aghozlane/scatterD3')
  library(scatterD3)
}

if (!require(rNVD3)) {
  devtools::install_github('aghozlane/rNVD3')
  library(rNVD3)
}

if (!require(genefilter)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("genefilter")
  library(genefilter)
}

if (!require(googleVis)) {
  install.packages('googleVis')
  suppressPackageStartupMessages(library(googleVis))
}

if (!require(shinyjs)) {
  install.packages('shinyjs')
  library(shinyjs)
}

if(!require(d3vennR)){
  install_github("timelyportfolio/d3vennR")
  library(d3vennR)
}

if (!require(DT)) {
  install.packages('DT')
  library(DT)
}

if (!require(RColorBrewer)) {
  install.packages('RColorBrewer')
  library(RColorBrewer)
}

if (!require(gplots)) {
  install.packages('gplots')
  library(gplots)
}

if (!require(DESeq2)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  library(DESeq2)
}

if (!require(ade4)) {
  install.packages('ade4')
  library(ade4)
}

if (!require(scales)) {
  install.packages('scales')
  library(scales)
}
# if(!require(plotly)){
#   install.packages("plotly")
#   library(plotly)  
# }

# Allow to upload 50M files only shaman server
if(Sys.info()["nodename"] == "shaman"){
  options(shiny.maxRequestSize=50*1024^2)
}else{
  # No limit
  options(shiny.maxRequestSize=500000000000000*1024^2)
}
source("Rfunctions/Data_Management.R")
source("Rfunctions/Stat_Model.R")
source("Rfunctions/DiagPlot.R")
source("Rfunctions/VisuPlot.R")
source("Rfunctions/CompPlot.R")
source("Rfunctions/DiffTable.R")


