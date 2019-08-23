#options(download.file.method = 'wget', bitmapType='cairo')
options(bitmapType='cairo')

if (!require("Rcpp")){
  install.packages("Rcpp",  repos="https://cran.univ-paris1.fr/")
}

if(!require(shiny)){
  install.packages("shiny",  repos="https://cran.univ-paris1.fr/")
}

if(!require(rjson)){
  install.packages('rjson',  repos="https://cran.univ-paris1.fr/")
}

if(!require(ape)){
  install.packages('ape',  repos="https://cran.univ-paris1.fr/")
}

if(!require(GUniFrac)){
  install.packages('GUniFrac',  repos="https://cran.univ-paris1.fr/")
}

if (!require(psych)) {
  install.packages('psych',  repos="https://cran.univ-paris1.fr/")
  library(psych)
}

if (!require(ggplot2)) {
  install.packages('ggplot2',  repos="https://cran.univ-paris1.fr/")
}

if (!require(vegan)) {
  install.packages('vegan',  repos="https://cran.univ-paris1.fr/")
}

if (!require(dendextend)) {
  install.packages('dendextend',  repos="https://cran.univ-paris1.fr/")
}

if (!require(circlize)) {
  install.packages('circlize',  repos="https://cran.univ-paris1.fr/")
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
  install.packages('googleVis',  repos="https://cran.univ-paris1.fr/")
  #suppressPackageStartupMessages(library(googleVis))
}

if (!require(shinyjs)) {
  install.packages('shinyjs',  repos="https://cran.univ-paris1.fr/")
}

if (!require(DT)) {
  install.packages('DT',  repos="https://cran.univ-paris1.fr/")
}

if (!require(RColorBrewer)) {
  install.packages('RColorBrewer',  repos="https://cran.univ-paris1.fr/")
}

if (!require(gplots)) {
  install.packages('gplots',  repos="https://cran.univ-paris1.fr/")
}

if (!require(ade4)) {
  install.packages('ade4',  repos="https://cran.univ-paris1.fr/")
}

if (!require(scales)) {
  install.packages('scales',  repos="https://cran.univ-paris1.fr/")
}

if (!require(phytools)) {
  install.packages('phytools',  repos="https://cran.univ-paris1.fr/")
}

if(!require(philentropy)){
  install.packages("philentropy",  repos="https://cran.univ-paris1.fr/")
}

if (!require("shinyWidgets")){
  install.packages("shinyWidgets")
  #devtools::install_github("aghozlane/shinyWidgets")
}

if (!require("sendmailR")){
  install.packages("sendmailR",  repos="https://cran.univ-paris1.fr/")
}

if (!require("shinyBS")){
  install.packages("shinyBS",  repos="https://cran.univ-paris1.fr/")
}

library(tools)

if (!require("flexdashboard")){
  install.packages("flexdashboard",  repos="https://cran.univ-paris1.fr/")
}

if (!require("backports")){
  install.packages("backports",  repos="https://cran.univ-paris1.fr/")
}


if (!require("readr")){
  install.packages("readr",  repos="https://cran.univ-paris1.fr/")
}

if (!require("jsonlite")){
  install.packages("jsonlite",  repos="https://cran.univ-paris1.fr/")
}

if (!require("shinyFiles")){
  install.packages("shinyFiles",  repos="https://cran.univ-paris1.fr/")
}

if (!require("htmltools")){
  install.packages("htmltools",  repos="https://cran.univ-paris1.fr/")
}

if (!require("rAmCharts")){
  install.packages("rAmCharts",  repos="https://cran.univ-paris1.fr/")
}

if(!require("colourpicker")){
  install.packages("colourpicker",  repos="https://cran.univ-paris1.fr/")
}

if(!require("data.table")){
  install.packages("data.table",  repos="https://cran.univ-paris1.fr/")
}

if(!require("UpSetR")){
  install.packages("UpSetR",  repos="https://cran.univ-paris1.fr/")
}

if(!require("ggrepel")){
  install.packages("ggrepel",  repos="https://cran.univ-paris1.fr/")
}

if(!require("igraph")){
  install.packages("igraph",  repos="https://cran.univ-paris1.fr/")
}

if(!require("visNetwork")){
  install.packages("visNetwork",  repos="https://cran.univ-paris1.fr/")
}

if (!require("shinytoastr")){
  install.packages("shinytoastr",  repos="https://cran.univ-paris1.fr/")
}

if (!require("scatterD3")) {
  #devtools::install_github('aghozlane/scatterD3')
  install.packages("scatterD3", repos="https://cran.univ-paris1.fr/")
}

if (!require(devtools)) {
  install.packages('devtools',  repos="https://cran.univ-paris1.fr/")
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
            "googleVis")
lapply(libshaman, require, character.only = TRUE)