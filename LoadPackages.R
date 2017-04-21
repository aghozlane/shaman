options(download.file.method = 'wget')

if (!require("backports")){
  install.packages("backports")
  library(backports)
}

if (!require("backports")){
  install.packages("backports")
  
}

if (!require("readr")){
  install.packages("readr")
  library(readr)
}

if (!require("jsonlite")){
  install.packages("jsonlite")
  library(jsonlite)
}


if (!require("shinyFiles")){
  install.packages("shinyFiles")
  library(shinyFiles)
}


if (!require("Rcpp")){
  install.packages("Rcpp")
}
if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}
if(!require(shinydashboard)){
  devtools::install_github('aghozlane/shinydashboard')
  library(shinydashboard)
}

if(!require(rjson)){
  install.packages('rjson')
  library(rjson)
}
if(!require(ape)){
  install.packages('ape')
  library(ape)
}

if(!require(GUniFrac)){
  install.packages('GUniFrac')
  library(GUniFrac)
}

if(!require(devtools)){
  install.packages('devtools')
  library(devtools)
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
if(!require(treeWeightD3)){
   devtools::install_github('pierreLec/treeWeightD3')
   library(treeWeightD3)
}
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

if (!require(phytools)) {
  install.packages('phytools')
  library(phytools)
}

if(!require(philentropy)){
  devtools::install_github("HajkD/philentropy", build_vignettes = TRUE, dependencies = TRUE)
  library(philentropy)
}

if(!require(PhyloTreeMetaR)){
  devtools::install_github("pierreLec/PhyloTreeMetaR")
  library(PhyloTreeMetaR)
}


# if(!require(plotly)){
#   install.packages("plotly")
#   library(plotly)  
# }


