options(download.file.method = 'wget', bitmapType='cairo')

if (!require("Rcpp")){
  install.packages("Rcpp",  repos="http://cran.univ-paris1.fr/")
}
if(!require(shiny)){
  install.packages("shiny",  repos="http://cran.univ-paris1.fr/")
  library(shiny)
}
if(!require(shinydashboard)){
  devtools::install_github('aghozlane/shinydashboard')
  library(shinydashboard)
}

if(!require(rjson)){
  install.packages('rjson',  repos="http://cran.univ-paris1.fr/")
  library(rjson)
}
if(!require(ape)){
  install.packages('ape',  repos="http://cran.univ-paris1.fr/")
  library(ape)
}

if(!require(GUniFrac)){
  install.packages('GUniFrac',  repos="http://cran.univ-paris1.fr/")
  library(GUniFrac)
}

if(!require(devtools)){
  install.packages('devtools',  repos="http://cran.univ-paris1.fr/")
  library(devtools)
}

if (!require(psych)) {
  install.packages('psych',  repos="http://cran.univ-paris1.fr/")
  library(psych)
}

if (!require(ggplot2)) {
  install.packages('ggplot2',  repos="http://cran.univ-paris1.fr/")
  library(ggplot2)
}

if (!require(vegan)) {
  install.packages('vegan',  repos="http://cran.univ-paris1.fr/")
  library(vegan)
}

if (!require(dendextend)) {
  install.packages('dendextend',  repos="http://cran.univ-paris1.fr/")
  library(dendextend)
}

if (!require(circlize)) {
  install.packages('circlize',  repos="http://cran.univ-paris1.fr/")
  library(circlize)
}


if (!require(devtools)) {
  install.packages('devtools',  repos="http://cran.univ-paris1.fr/")
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
  install.packages('googleVis',  repos="http://cran.univ-paris1.fr/")
  suppressPackageStartupMessages(library(googleVis))
}

if (!require(shinyjs)) {
  install.packages('shinyjs',  repos="http://cran.univ-paris1.fr/")
  library(shinyjs)
}

if(!require(d3vennR)){
  install_github("timelyportfolio/d3vennR")
  library(d3vennR)
}

if (!require(DT)) {
  install.packages('DT',  repos="http://cran.univ-paris1.fr/")
  library(DT)
}

if (!require(RColorBrewer)) {
  install.packages('RColorBrewer',  repos="http://cran.univ-paris1.fr/")
  library(RColorBrewer)
}

if (!require(gplots)) {
  install.packages('gplots',  repos="http://cran.univ-paris1.fr/")
  library(gplots)
}

if (!require(DESeq2)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  library(DESeq2)
}

if (!require(ade4)) {
  install.packages('ade4',  repos="http://cran.univ-paris1.fr/")
  library(ade4)
}

if (!require(scales)) {
  install.packages('scales',  repos="http://cran.univ-paris1.fr/")
  library(scales)
}

if (!require(phytools)) {
  install.packages('phytools',  repos="http://cran.univ-paris1.fr/")
  library(phytools)
}

if(!require(philentropy)){
  install.packages("philentropy",  repos="http://cran.univ-paris1.fr/")
  library(philentropy)
}

if(!require(PhyloTreeMetaR)){
  devtools::install_github("pierreLec/PhyloTreeMetaR")
  library(PhyloTreeMetaR)
}

if (!require("shinytoastr")){
  devtools::install_github("mangothecat/shinytoastr")
  library(shinytoastr)
}

if (!require("shinyWidgets")){
  devtools::install_github("aghozlane/shinyWidgets")
  library(shinyWidgets)
}

if (!require("sendmailR")){
  install.packages("sendmailR",  repos="http://cran.univ-paris1.fr/")
  library(sendmailR)
}

if (!require("shinyBS")){
  install.packages("shinyBS",  repos="http://cran.univ-paris1.fr/")
  library(shinyBS)
}

library(tools)

if (!require("flexdashboard")){
  install.packages("flexdashboard",  repos="http://cran.univ-paris1.fr/")
  library(flexdashboard)
}

if (!require("backports")){
  install.packages("backports",  repos="http://cran.univ-paris1.fr/")
  library(backports)
}


if (!require("readr")){
  install.packages("readr",  repos="http://cran.univ-paris1.fr/")
  library(readr)
}

if (!require("jsonlite")){
  install.packages("jsonlite",  repos="http://cran.univ-paris1.fr/")
  library(jsonlite)
}


if (!require("shinyFiles")){
  install.packages("shinyFiles",  repos="http://cran.univ-paris1.fr/")
  library(shinyFiles)
}

if (!require("htmltools")){
  install.packages("htmltools",  repos="http://cran.univ-paris1.fr/")
  library(htmltools)
}

# if (!require("V8")){
#   install.packages("V8",  repos="http://cran.univ-paris1.fr/")
#   library(htmltools)
# }

# if(!require(plotly)){
#   install.packages("plotly")
#   library(plotly)  
# }


