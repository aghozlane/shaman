if (!require("Rcpp")){
  install.packages("Rcpp")
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

#library(plotly)

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

if (!require(d3heatmap)) {
  install.packages('d3heatmap')
  library(d3heatmap)
}

if (!require(devtools)) {
  install.packages('devtools')
  library(devtools)
}

# Let us use biomformat instead of biom
if (!require(biomformat)){
  library(devtools)
 devtools::install_github("biomformat", "aghozlane")
library(biomformat)
}

if (!require(scatterD3)) {
  install.packages('scatterD3')
  library(scatterD3)
}

if (!require(rNVD3)) {
  library(devtools)
  install_github('rNVD3', 'aghozlane')
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

if(!require(plotly)){
  install.packages('plotly')
  library(plotly)
}

if(!require(d3vennR)){
  install_github("timelyportfolio/d3vennR")
  library(d3vennR)
}

if (!require(psych)) {
  install.packages('psych')
  library(psych)
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


if(!require(biomformat)){
  devtools::install_github("biomformat", "aghozlane")
}


# Allow to upload 50M files
options(shiny.maxRequestSize=50*1024^2) 
source("Rfunctions/Data_Management.R")
source("Rfunctions/Stat_Model.R")
source("Rfunctions/DiagPlot.R")
source("Rfunctions/VisuPlot.R")
source("Rfunctions/CompPlot.R")
source("Rfunctions/DiffTable.R")


