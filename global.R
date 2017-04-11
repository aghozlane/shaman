# Allow to upload 50M files only shaman server
if(Sys.info()["nodename"] == "ShinyPro"){
  options(shiny.maxRequestSize=50*1024^2)
}else{
  # No limit
  options(shiny.maxRequestSize=500000000000000*1024^2)
}


source('LoadPackages.R')
source("css/owncss.R")
source("Rfunctions/Data_Management.R")
source("Rfunctions/Stat_Model.R")
source("Rfunctions/DiagPlot.R")
source("Rfunctions/VisuPlot.R")
source("Rfunctions/CompPlot.R")
source("Rfunctions/DiffTable.R")
source('Rfunctions/directoryInput.R')
