#options(download.file.method = 'wget', bitmapType='cairo')
options(bitmapType='cairo')

if (!require("Rcpp")){
  install.packages("Rcpp",  repos="https://cran.univ-paris1.fr/")
}
if(!require(shiny)){
  install.packages("shiny",  repos="https://cran.univ-paris1.fr/")
  library(shiny)
}

if(!require(rjson)){
  install.packages('rjson',  repos="https://cran.univ-paris1.fr/")
  library(rjson)
}
if(!require(ape)){
  install.packages('ape',  repos="https://cran.univ-paris1.fr/")
  library(ape)
}

if(!require(GUniFrac)){
  install.packages('GUniFrac',  repos="https://cran.univ-paris1.fr/")
  library(GUniFrac)
}

if (!require(psych)) {
  install.packages('psych',  repos="https://cran.univ-paris1.fr/")
  library(psych)
}

if (!require(ggplot2)) {
  install.packages('ggplot2',  repos="https://cran.univ-paris1.fr/")
  library(ggplot2)
}

if (!require(vegan)) {
  install.packages('vegan',  repos="https://cran.univ-paris1.fr/")
  library(vegan)
}

if (!require(dendextend)) {
  install.packages('dendextend',  repos="https://cran.univ-paris1.fr/")
  library(dendextend)
}

if (!require(circlize)) {
  install.packages('circlize',  repos="https://cran.univ-paris1.fr/")
  library(circlize)
}


# if (!require(BiocInstaller)){
#   source("https://bioconductor.org/biocLite.R")
#   biocLite("BiocInstaller")
#   library(BiocInstaller)
# }

if(!require(BiocManager)){
  install.packages("BiocManager")
}
if (!require(genefilter)) {
  BiocManager::install("genefilter", ask=FALSE)
  library(genefilter)
}

if(!require("SummarizedExperiment")){
  BiocManager::install("SummarizedExperiment", ask=FALSE)
  library(SummarizedExperiment)
}

if (!require(googleVis)) {
  install.packages('googleVis',  repos="https://cran.univ-paris1.fr/")
  suppressPackageStartupMessages(library(googleVis))
}

if (!require(shinyjs)) {
  install.packages('shinyjs',  repos="https://cran.univ-paris1.fr/")
  library(shinyjs)
}

if (!require(DT)) {
  install.packages('DT',  repos="https://cran.univ-paris1.fr/")
  library(DT)
}

if (!require(RColorBrewer)) {
  install.packages('RColorBrewer',  repos="https://cran.univ-paris1.fr/")
  library(RColorBrewer)
}

if (!require(gplots)) {
  install.packages('gplots',  repos="https://cran.univ-paris1.fr/")
  library(gplots)
}

if (!require(ade4)) {
  install.packages('ade4',  repos="https://cran.univ-paris1.fr/")
  library(ade4)
}

if (!require(scales)) {
  install.packages('scales',  repos="https://cran.univ-paris1.fr/")
  library(scales)
}

if (!require(phytools)) {
  install.packages('phytools',  repos="https://cran.univ-paris1.fr/")
  library(phytools)
}

if(!require(philentropy)){
  install.packages("philentropy",  repos="https://cran.univ-paris1.fr/")
  library(philentropy)
}

if (!require("shinyWidgets")){
  install.packages("shinyWidgets")
  #devtools::install_github("aghozlane/shinyWidgets")
  library(shinyWidgets)
}

if (!require("sendmailR")){
  install.packages("sendmailR",  repos="https://cran.univ-paris1.fr/")
  library(sendmailR)
}

if (!require("shinyBS")){
  install.packages("shinyBS",  repos="https://cran.univ-paris1.fr/")
  library(shinyBS)
}

library(tools)

if (!require("flexdashboard")){
  install.packages("flexdashboard",  repos="https://cran.univ-paris1.fr/")
  library(flexdashboard)
}

if (!require("backports")){
  install.packages("backports",  repos="https://cran.univ-paris1.fr/")
  library(backports)
}


if (!require("readr")){
  install.packages("readr",  repos="https://cran.univ-paris1.fr/")
  library(readr)
}

if (!require("jsonlite")){
  install.packages("jsonlite",  repos="https://cran.univ-paris1.fr/")
  library(jsonlite)
}


if (!require("shinyFiles")){
  install.packages("shinyFiles",  repos="https://cran.univ-paris1.fr/")
  library(shinyFiles)
}

if (!require("htmltools")){
  install.packages("htmltools",  repos="https://cran.univ-paris1.fr/")
  library(htmltools)
}

# if (!require("V8")){
#   install.packages("V8",  repos="https://cran.univ-paris1.fr/")
#   library(htmltools)
# }

# if(!require(plotly)){
#   install.packages("plotly")
#   library(plotly)  
# }

if (!require("rAmCharts")){
  install.packages("rAmCharts",  repos="https://cran.univ-paris1.fr/")
  library(rAmCharts)
}

if(!require("colourpicker")){
  install.packages("colourpicker",  repos="https://cran.univ-paris1.fr/")
  library(colourpicker)
}

if(!require("data.table")){
  install.packages("data.table",  repos="https://cran.univ-paris1.fr/")
  library(data.table)
}

if(!require("UpSetR")){
  install.packages("UpSetR",  repos="https://cran.univ-paris1.fr/")
  library(UpSetR)
}

if(!require("ggrepel")){
  install.packages("ggrepel",  repos="https://cran.univ-paris1.fr/")
  library(ggrepel)
}

# if(!require("networkD3")){
#   install.packages("networkD3")
#   library(networkD3)
# }

if(!require("igraph")){
  install.packages("igraph",  repos="https://cran.univ-paris1.fr/")
  library(igraph)
}

if(!require("visNetwork")){
  install.packages("visNetwork",  repos="https://cran.univ-paris1.fr/")
  library(visNetwork)
}

if (!require("shinytoastr")){
  install.packages("shinytoastr")
  library(shinytoastr)
}

if (!require("scatterD3")) {
  #devtools::install_github('aghozlane/scatterD3')
  install.packages("scatterD3", repos="https://cran.univ-paris1.fr/")
  library(scatterD3)
}

if (!require(devtools)) {
  install.packages('devtools',  repos="https://cran.univ-paris1.fr/")
  library(devtools)
}

if(!require("shinydashboardshaman")){
  devtools::install_github('aghozlane/shinydashboardshaman')
  library(shinydashboardshaman)
}

if (!require("d3heatmap")) {
  #devtools::install_github('aghozlane/d3heatmap')
  devtools::install_github("rstudio/d3heatmap")
  library(d3heatmap)
}

# Let us use biomformat instead of biom
#torename
if (!require("biomformatshaman")){
  devtools::install_github("aghozlane/biomformatshaman")
  library(biomformatshaman)
}

if (!require("rNVD3shaman")) {
  devtools::install_github('aghozlane/rNVD3shaman')
  library(rNVD3shaman)
}

if (!require("DESeq2shaman")) {
  devtools::install_github("aghozlane/DESeq2shaman")
  library(DESeq2shaman)
}

if(!require("PhyloTreeMetaR")){
  devtools::install_github("pierreLec/PhyloTreeMetaR")
  library(PhyloTreeMetaR)
}

if(!require("treeWeightD3")){
  devtools::install_github('pierreLec/treeWeightD3')
  library(treeWeightD3)
}

if(!require("d3vennR")){
  devtools::install_github("timelyportfolio/d3vennR")
  library(d3vennR)
}