# Allow to upload 50M files only shaman server
#if(Sys.info()["nodename"] == "ShinyPro"){
#  options(shiny.maxRequestSize=1000000000*1024^2)
#}else{
# Limit with the raw data submission to 2Gb
options(shiny.maxRequestSize=2000000000)
#}


source('LoadPackages.R')
source("css/owncss.R")
source("Rfunctions/Data_Management.R")
source("Rfunctions/Stat_Model.R")
source("Rfunctions/DiagPlot.R")
source("Rfunctions/VisuPlot.R")
source("Rfunctions/CompPlot.R")
source("Rfunctions/DiffTable.R")
source('Rfunctions/directoryInput.R')
source('Rfunctions/internal_masque.R')
