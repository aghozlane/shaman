#options(download.file.method = 'wget', bitmapType='cairo')
options(bitmapType='cairo')

if (!require("Rcpp")){
  #RcppArmadillo_0.9.800.3.0
  install.packages("Rcpp",  repos="http://cran.irsn.fr/")
}
# Limited shiny 1.3.2
if(!require(shiny)){
  install.packages("shiny",  repos="http://cran.irsn.fr/")
}

if(!require(rjson)){
  install.packages('rjson',  repos="http://cran.irsn.fr/")
}

if(!require(ape)){
  install.packages('ape',  repos="http://cran.irsn.fr/")
}

if(!require(GUniFrac)){
  install.packages('GUniFrac',  repos="http://cran.irsn.fr/")
}

if (!require(psych)) {
  install.packages('psych',  repos="http://cran.irsn.fr/")
  library(psych)
}

if (!require(ggplot2)) {
  install.packages('ggplot2',  repos="http://cran.irsn.fr/")
}

if (!require(vegan)) {
  install.packages('vegan',  repos="http://cran.irsn.fr/")
}

if (!require(dendextend)) {
  install.packages('dendextend',  repos="http://cran.irsn.fr/")
}

if (!require(circlize)) {
  install.packages('circlize',  repos="http://cran.irsn.fr/")
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
}

if(!require("SummarizedExperiment")){
  BiocManager::install("SummarizedExperiment", ask=FALSE)
}

if (!require(googleVis)) {
  install.packages('googleVis',  repos="http://cran.irsn.fr/")
  #suppressPackageStartupMessages(library(googleVis))
}

if (!require(shinyjs)) {
  install.packages('shinyjs',  repos="http://cran.irsn.fr/")
}

if (!require(DT)) {
  install.packages('DT',  repos="http://cran.irsn.fr/")
}

if (!require(RColorBrewer)) {
  install.packages('RColorBrewer',  repos="http://cran.irsn.fr/")
}

if (!require(gplots)) {
  install.packages('gplots',  repos="http://cran.irsn.fr/")
}

if (!require(ade4)) {
  install.packages('ade4',  repos="http://cran.irsn.fr/")
}

if (!require(scales)) {
  install.packages('scales',  repos="http://cran.irsn.fr/")
}

if (!require(phytools)) {
  install.packages('phytools',  repos="http://cran.irsn.fr/")
}

if(!require(philentropy)){
  install.packages("philentropy",  repos="http://cran.irsn.fr/")
}

if (!require("shinyWidgets")){
  install.packages("shinyWidgets", repos="http://cran.irsn.fr/")
  #devtools::install_github("aghozlane/shinyWidgets")
}

if (!require("sendmailR")){
  install.packages("sendmailR",  repos="http://cran.irsn.fr/")
}

if (!require("shinyBS")){
  install.packages("shinyBS",  repos="http://cran.irsn.fr/")
}

library(tools)

if (!require("flexdashboard")){
  install.packages("flexdashboard",  repos="http://cran.irsn.fr/")
}

if (!require("backports")){
  install.packages("backports",  repos="http://cran.irsn.fr/")
}


if (!require("readr")){
  install.packages("readr",  repos="http://cran.irsn.fr/")
}

if (!require("jsonlite")){
  install.packages("jsonlite",  repos="http://cran.irsn.fr/")
}

if (!require("shinyFiles")){
  install.packages("shinyFiles",  repos="http://cran.irsn.fr/")
}

if (!require("htmltools")){
  install.packages("htmltools",  repos="http://cran.irsn.fr/")
}

if (!require("rAmCharts")){
  install.packages("rAmCharts",  repos="http://cran.irsn.fr/")
}

if(!require("colourpicker")){
  install.packages("colourpicker",  repos="http://cran.irsn.fr/")
}

if(!require("data.table")){
  install.packages("data.table",  repos="http://cran.irsn.fr/")
}

if(!require("UpSetR")){
  install.packages("UpSetR",  repos="http://cran.irsn.fr/")
}

if(!require("ggrepel")){
  install.packages("ggrepel",  repos="http://cran.irsn.fr/")
}

if(!require("igraph")){
  install.packages("igraph",  repos="http://cran.irsn.fr/")
}

if(!require("visNetwork")){
  install.packages("visNetwork",  repos="http://cran.irsn.fr/")
}

if (!require("shinytoastr")){
  install.packages("shinytoastr",  repos="http://cran.irsn.fr/")
}

if (!require("scatterD3")) {
  #devtools::install_github('aghozlane/scatterD3')
  install.packages("scatterD3", repos="http://cran.irsn.fr/")
}

if (!require("stringr")) {
  install.packages("stringr", repos="http://cran.irsn.fr/")
}

if(!require("FactoMineR")){
  install.packages("FactoMineR", repos="http://cran.irsn.fr/")
}

if(!require("factoextra")){
  install.packages("factoextra", repos="http://cran.irsn.fr/")
}

if (!require(devtools)) {
  install.packages('devtools',  repos="http://cran.irsn.fr/")
}

if(!require("shinydashboardshaman")){
  devtools::install_github('aghozlane/shinydashboardshaman')
}

if (!require("d3heatmap")) {
  #devtools::install_github('aghozlane/d3heatmap')
  devtools::install_github("rstudio/d3heatmap")
}

if (!require("biomformatshaman")){
  devtools::install_github("aghozlane/biomformatshaman")
}

if (!require("rNVD3shaman")) {
  devtools::install_github('aghozlane/rNVD3shaman')
}

if (!require("DESeq2shaman")) {
  devtools::install_github("aghozlane/DESeq2shaman")
}

if(!require("PhyloTreeMetaR")){
  devtools::install_github("pierreLec/PhyloTreeMetaR")
}

if(!require("treeWeightD3")){
  devtools::install_github('pierreLec/treeWeightD3')
}

if(!require("d3vennR")){
  devtools::install_github("timelyportfolio/d3vennR")
}

libshaman=c("shiny", "rjson", "ape", "GUniFrac", "psych", "ggplot2", "vegan", "dendextend", 
            "circlize", "genefilter", "SummarizedExperiment", "shinyjs", "DT", "RColorBrewer", 
            "gplots", "ade4", "scales", "phytools", "philentropy", "shinyWidgets", "sendmailR", 
            "shinyBS", "tools", "flexdashboard", "backports", "readr", "jsonlite", "shinyFiles", 
            "htmltools", "rAmCharts", "colourpicker", "data.table", "UpSetR", "ggrepel", "igraph", 
            "visNetwork", "shinytoastr", "scatterD3", "devtools", "shinydashboardshaman", "d3heatmap", 
            "biomformatshaman", "rNVD3shaman", "DESeq2shaman", "PhyloTreeMetaR", "treeWeightD3", "d3vennR",
            "googleVis", "stringr", "FactoMineR", "factoextra")
lapply(libshaman, require, character.only = TRUE)