library(shiny)
library(psych)
library(ggplot2)
#library(gdata)
#install_github('rCharts', 'ramnathv')
source("internal.R")

renderDataTable <- DT::renderDataTable
dataTableOutput <- DT::dataTableOutput

shinyServer(function(input, output,session) {
  
  
#####################################################
##
##                    LOAD FILES
##
#####################################################
  
  namesfile = "www/BaseContrast.txt"
  #file.create(namesfile,showWarnings=FALSE)
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
    
    data = read.csv(inFile$datapath,sep="\t",header=TRUE)

    ## Rownames
    rownames(data)=data[,1];data=data[,-1]
    
    return(as.data.frame(data))
  })
  


  ## Taxo File
  dataInputTaxo <-reactive({ 
    
    inFile <- input$fileTaxo
    
    if (is.null(inFile)) return(NULL)
    

    data = read.csv(inFile$datapath,sep="\t",header=TRUE)
    
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
      if(!is.null(tmpBIOM)) data = GetDataFromBIOM(tmpBIOM)
    }
    
    return(data)
  })
  
  
  
  
#####################################################
##
##                DYNAMIC MENU
##
#####################################################
  
  
  
  output$dymMenu <- renderMenu({
    
    data = dataInput()
    
#     if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0)
#     {
      sidebarMenu(
        menuItem("Statistical analysis",
                 menuSubItem("Run differential analysis",tabName="RunDiff"),
                 menuSubItem("Diagnostic plots",tabName="DiagPlot"),
                 menuSubItem("Tables",tabName="TableDiff"),
                 icon = icon("bar-chart-o"), tabName = "AnaStat"),
        menuItem("Data visualisation", icon = icon("area-chart"), tabName = "Visu"),
        menuItem("Krona plot", icon = icon("pie-chart"), tabName = "Krona")

      )
#     }
  })
  
  

  #####################################################
  ##
  ##                DATA TABLE
  ##
  #####################################################
  
  ## Counts Table
  output$DataCounts <- renderDataTable(
    dataInput()$counts, 
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE
    ))
  
  ## Taxonomy table
  output$DataTaxo <- renderDataTable(
    dataInput()$taxo, 
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE
    ))

  
  ## Tab box for visualisation
  output$TabBoxData <- renderUI({
    
    data=dataInput()
    
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0)
    {
      tabBox(width = NULL, selected = "Counts table",
             tabPanel("Counts table",dataTableOutput("DataCounts")),
             tabPanel("Taxonomy",dataTableOutput("DataTaxo"))  
      )
    }
    
  })


  #####################################################
  ##
  ##                TARGET FILE
  ##
  #####################################################
  

  ## Load target file
  dataInputTarget <-reactive({ 
    
    inFile <- input$fileTarget
    
    if (is.null(inFile)) return(NULL)
    
    
    data = read.csv(inFile$datapath,sep="\t",header=TRUE)
    
    return(as.data.frame(data))
  })


  # Infobox design
  output$RowTarget <- renderInfoBox({
    
    target = dataInputTarget()
    
    InterVar = input$InterestVar
    Interaction = input$Interaction2
    alltmp = c(InterVar,Interaction)
    
    if(!is.null(target)) 
    {
      #### Ajout fontion check target
      infoBox(h6(strong("Target format")), subtitle = h6("Your target file is OK"), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
    }
    else infoBox(h6(strong("Warning")), subtitle = h6("Label of the target file must correspond to counts table column names") ,color = "orange",width=NULL,fill=TRUE, icon = icon("warning"))
  })


  ## Interest Variables
  output$SelectInterestVar <- renderUI({
        
    target=dataInputTarget()
    
    if(!is.null(target)) 
    {
      namesTarget = colnames(target)[2:ncol(target)]
      selectInput("InterestVar",h6(strong("Select the variables")),namesTarget,selected=namesTarget,multiple=TRUE)
    }
  
  })


  ## Interactions
  output$SelectInteraction2 <- renderUI({
        
    target = dataInputTarget()

    if(!is.null(target)) 
    {
      Interac = GetInteraction2(target)
      selectInput("Interaction2",h6(strong("Add interactions")),Interac,selected=NULL,multiple=TRUE)
    }
    
  })


  ## Reference radio buttons
  output$RefSelect <- renderUI({
    
    target = dataInputTarget()
    RB=list()
    if(!is.null(target)) 
    {
      InterVar = input$InterestVar
      if(length(InterVar)>0)
      {
        names = paste0("Ref",InterVar)
          for(i in 1:length(names))
          {
            val = unique(target[,InterVar[i]])
            RB[[i]] = selectInput(names[i],paste("Reference for",InterVar[i]),as.vector(val))
          }
      }
    }
    return(RB)
    
  })



#####################################################
##
##            DEFINE CONTRAST
##
#####################################################

  
  output$contrastMat <- renderUI({
    
    
    #dds=RunDESeq2()
    #names = resultsNames(dds)
    names=c('test1',"test2","test3")
    Contrast=list()
    
    for(i in 1:length(names)){Contrast[[i]] = textInput(names[i],names[i],0)}
  
    return(Contrast)
    
  })


  output$ContrastOverview <- renderUI({
    
    #dds=RunDESeq2()
    #names = resultsNames(dds)
    names=c('test1',"test2","test3")
    cont = input$ContrastList
    ContrastBase = read.table(namesfile,header=TRUE)
    
    res = PrintContrasts(names,ContrastBase[,cont])
    return(res)
  })



  BaseContrast <- function(input,namesfile)
  {  
    #dds=RunDESeq2()
    #names = resultsNames(dds)
    
    oldContrast = read.table(namesfile,header=TRUE)
    names=c('test1',"test2","test3")
    v_tmp = c()
    
    for(i in 1:length(names))
    {  
      Tinput = paste("input$",names[i],sep="")
      print(Tinput)
      print(input$test1)
      
      expr=parse(text=Tinput)
      val = eval(expr) 
      print(val)
      v_tmp[i] = as.numeric(val)
    }
    print(v_tmp)
    colnamesTmp = c(colnames(oldContrast),input$ContrastName)
    mat = cbind(oldContrast,v_tmp)
    write.table(mat,namesfile,row.names=FALSE,col.names = colnamesTmp)
  }


  AddCont <-eventReactive(input$AddContrast,{
    
    
    BaseContrast(input,namesfile)
    tmp = read.table(namesfile,header=TRUE)
    Contrast = colnames(tmp)
    updateSelectInput(session, "ContrastList","Contrasts",Contrast)
    
  
  })

  ## 
  observeEvent(input$AddContrast,{  
    
    AddCont()
    
  })


  
  ## Get the results from MEMHDX
  RemoveCont <-eventReactive(input$RemoveContrast,{
    
    tmp = read.table(namesfile,header=TRUE)
    print(input$ContrastList)
    matKept = tmp[,-which(colnames(tmp)%in%input$ContrastList)]
    ContrastKept = colnames(matKept)
    print(ContrastKept)
    write.table(matKept,namesfile,row.names=FALSE,col.names = ContrastKept)
    
    updateSelectInput(session, "ContrastList","Contrasts",ContrastKept)
  })
  
  
  
  ## Run MEMHDX via RunProcess button
  observeEvent(input$RemoveContrast,{  
    
    RemoveCont()
    
  })

#####################################################
##
##                OPTIONS DIFF ANALYSIS
##
#####################################################



#   ## Select variable for reference
#   output$RadioSelectVarRef <- renderUI({
#     
#     target = dataInputTarget()
#     Var = input$InterestVar
#     
#     if(!is.null(target)) 
#     {
#       radioButtons("SelectVarRef",h6(strong("Select the reference for each variable")),Var)
#     }
#     
#   })
# 
#   
#   ## Select the reference
#   output$SelectVarRef <- renderUI({
#     
#     target = dataInputTarget()
#     Var = input$InterestVar
#     
#     if(!is.null(target)) 
#     {
#       selectInput("SelectRef",h6(strong("Select the reference for each variable")),Var)
#     }
#     
#     
# 
#     if(!is.null(target)) 
#     {
#       mod = target[,input$SelectVarRef]
#       selectInput("SelectRef",h6(strong("")),mod)
#     }
#     
#   })

#   # Infobox design
#   output$DesignBox <- renderInfoBox({
#     
#     target = dataInputTarget()
#     
#     InterVar = input$InterestVar
#     Interaction = input$Interaction2
#     alltmp = c(InterVar,Interaction)
#     
#     if(!is.null(target)) 
#     {
#       design = paste0(alltmp, collapse= "+")
#       infoBox(h6(strong("Design")), subtitle = as.expression(design), icon = icon("info"),color = "green",width=NULL)
#     }
#   })


})