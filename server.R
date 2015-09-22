library(shiny)
library(psych)
library(ggplot2)
#library(gdata)
#install_github('rCharts', 'ramnathv')
source("internal.R")

renderDataTable <- DT::renderDataTable
dataTableOutput <- DT::dataTableOutput

shinyServer(function(input, output) {
  
  
#####################################################
##
##                    LOAD FILES
##
#####################################################


  ## Counts file
  dataInputCounts <-reactive({ 
    
    inFile <- input$fileCounts
  
    if (is.null(inFile)) return(NULL)
    
    ## Get the extension
#     tmp = strsplit(inFile$name, ".",fixed=T)[[1]]
#     ext = tmp[length(tmp)]
#     
#     ## header
#     header = FALSE
#     if(input$header==1) header=TRUE
    
    ## Read data
#     if(ext=="csv") data = read.csv(inFile$datapath,sep=",",header=header)
#     if(ext=="xls") data = read.csv(inFile$datapath,sep="\t",header=header)
    
    data = read.csv(inFile$datapath,sep=",",header=TRUE)

    ## Rownames
    rownames(data)=data[,1];data=data[,-1]
    
    return(as.data.frame(data))
  })
  


  ## Taxo File
  dataInputTaxo <-reactive({ 
    
    inFile <- input$fileTaxo
    
    if (is.null(inFile)) return(NULL)
    

    data = read.csv(inFile$datapath,sep=",",header=TRUE)
    
    ## Rownames
    rownames(data)=data[,1];data=data[,-1]
    
    return(as.data.frame(data))
  })


  ## BIOM File
  dataInputBiom <-reactive({ 
    
    inFile <- input$fileBiom
    
    if (is.null(inFile)) return(NULL)
    
    
    data = read.csv(inFile$datapath,sep=",",header=TRUE)
    
    ## Rownames
    rownames(data)=data[,1];data=data[,-1]
    
    return(as.data.frame(data))
  })



  ## Input data
  dataInput <-reactive({ 
    
    data = NULL
    
    if(input$FileFormat=="fileCounts")
    {
      Counts = dataInputCounts()
      Taxo = dataInputTaxo()
      data = GetDataFromCT(Counts,Taxo)
    }
    
    if(input$FileFormat=="fileBiom")
    {
      tmpBIOM = dataInputBiom()
      if(!is.null(tmpBIOM)) data = GetDataFromBIOM(Counts,Taxo)
    }
    
    return(data)
  })
  
  
  
  
#####################################################
##
##                DYNAMIC MENU
##
#####################################################
  
  
  
  output$dymMenu <- renderMenu({
    
    input$LoadData
    data=isolate(dataInput())

    if(nrow(data$Counts)>0 && nrow(data$Taxo)>0)
    {
      sidebarMenu(
        menuItem("Statistical analysis",
                 menuSubItem("Run differential analysis",tabName="RunDiff"),
                 menuSubItem("Diagnostic plots",tabName="DiagPlot"),
                 menuSubItem("Tables",tabName="TableDiff"),
                 icon = icon("bar-chart-o"), tabName = "AnaStat"),
        menuItem("Data visualisation", icon = icon("area-chart"), tabName = "Visu"),
        menuItem("Krona plot", icon = icon("pie-chart"), tabName = "Krona")

      )
    }
  })
  
  

  output$LoadButton <- renderUI({
    
    data=isolate(dataInput())
    button = NULL
    
    if(nrow(data$Counts)>0 && nrow(data$Taxo)>0)
    {
      button = actionButton("LoadData",icon=icon("upload"),"Load your data")
    }
    
    return(button)
    
  })
  
    
})