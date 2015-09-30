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
  
  ## Create base for contrast
  rand = floor(runif(1,0,1e9))
  namesfile = paste("www/base/BaseContrast_",rand,".txt",sep="")
  file.create(namesfile,showWarnings=FALSE)

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
    
    
    ## Add NA
    data=as.matrix(data)
    indNa = which(data=="")
    data[indNa]=NA
    
    
    
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

  
  ## Tab box for data visualisation
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
    rownames(target) <- as.character(target[, 1])
    return((data))
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


  # Infobox design
  output$RowTarget <- renderInfoBox({
    
    target = dataInputTarget()
    
    if(!is.null(target)) 
    {
      #### Ajout fontion check target
      infoBox(h6(strong("Target file")), subtitle = h6("Your target file is OK"), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
    }
    else infoBox(h6(strong("Target file")), subtitle = h6("Label of the target file must correspond to counts table column names") ,color = "orange",width=NULL,fill=TRUE, icon = icon("warning"))
  })




  ## taget table
  output$DataTarget <- renderDataTable(
  dataInputTarget(),
  options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                 pageLength = 10,scrollX=TRUE
  ))


  ## Box for target visualisation
  output$BoxTarget <- renderUI({
    
    target = dataInputTarget()
    
    if(!is.null(target) &&  nrow(target)>0)
    {
      box(title="Target file overview",width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = TRUE,
          dataTableOutput("DataTarget")
      )  
    }
    
  })



#####################################################
##
##            DEFINE CONTRAST
##
#####################################################
  
  output$contrastMat <- renderUI({
    
    resDiff = ResDiffAnal()
    dds = resDiff$dds
    names = resultsNames(dds)

    Contrast=list()
    
    for(i in 1:length(names)){Contrast[[i]] = textInput(names[i],names[i],0)}
  
    return(Contrast)
    
  })


  output$ContrastOverview <- renderPrint({
    
    resDiff = ResDiffAnal()
    dds = resDiff$dds
    names = resultsNames(dds)
    
    cont = input$ContrastList
    filesize = file.info(namesfile)[,"size"]
    
    if(filesize!=0)
    { 
      ContrastBase = read.table(namesfile,header=TRUE)
      ind = which(colnames(ContrastBase)%in%cont)
      div(HTML(PrintContrasts(names,sapply(ContrastBase[,ind],as.numeric),cont)))
    }
  })


  ## Add contrast function
  AddCont <-eventReactive(input$AddContrast,{
    
    resDiff = ResDiffAnal()
    dds = resDiff$dds
    names = resultsNames(dds)
    
    BaseContrast(input,names,namesfile)
    tmp = read.table(namesfile,header=TRUE)
    Contrast = colnames(as.matrix(tmp))
    updateSelectInput(session, "ContrastList","Contrasts",Contrast)
    
  })

  ## Add contrast 
  observeEvent(input$AddContrast,{  
    
    AddCont()
    
  })

 
  
  ## Remove contrast function
  RemoveCont <-eventReactive(input$RemoveContrast,{
    
    ## get the size of the contrast base file
    filesize = file.info(namesfile)[,"size"]
    if(filesize!=0)
    { 
      tmp = read.table(namesfile,header=TRUE)
      ind = which(colnames(tmp)%in%input$ContrastList)
      matKept = as.matrix(tmp[,-ind])
      ContrastKept = colnames(tmp)[-ind]
      
      if(ncol(matKept)>0) write.table(matKept,namesfile,row.names=FALSE,col.names = ContrastKept)
      else file.create(namesfile,showWarnings=FALSE)
      updateSelectInput(session, "ContrastList","Contrasts",ContrastKept)
    }
  })
  
  
  
  ## Remove contrast
  observeEvent(input$RemoveContrast,{  
    
    RemoveCont()
    
  })


  # Infobox Contrast
  output$InfoContrast <- renderInfoBox({
    
    test = FALSE
    input$AddContrast
    
    filesize = isolate(file.info(namesfile)[,"size"])
    if(filesize!=0) 
    {
      tmp = read.table(namesfile,header=TRUE)
      if(any(as.vector(tmp)!=0)) test = TRUE
    }
    
    if(test) 
    {
      infoBox(h6(strong("Contrasts")), subtitle = h6("Contrasts OK"), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
    }
    else infoBox(h6(strong("Contrasts")), subtitle = h6("At least one contrast (non null) must be defined") ,color = "orange",width=NULL,fill=TRUE, icon = icon("warning"))
  })



#####################################################
##
##                DESEQ2 run
##
#####################################################



  # Infobox Contrast
  output$InfoDESeq <- renderInfoBox({
    
    input$RunDESeq
    box = NULL
    target = isolate(dataInputTarget())
    taxo = input$TaxoSelect

    if(!is.null(target) && taxo!="...") 
    {
      infoBox(h6(strong("Statistical analysis")), subtitle = h6("Differential analysis is done !"), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
    }
    else infoBox(h6(strong("Statistical analysis")), subtitle = h6("Not done !"), icon = icon("warning"),color = "orange",width=NULL,fill=TRUE) 
    
  })
  

  ## Get the results from DESeq2
  ResDiffAnal <-eventReactive(input$RunDESeq,{
    
    data = dataInput() 

    target = dataInputTarget()
    design = GetDesign(input)
    counts = GetCountsMerge(data,input$TaxoSelect)
   
    Get_dds_object(input,counts,target,design)

    
  })
  
  
  
  ## Run DESeq2 via RunDESeq button
  observeEvent(input$RunDESeq,{  
    
    ResDiffAnal()
    
  })

  
#####################################################
##
##                Taxonomy
##
#####################################################
  
  
  # Infobox Contrast
  output$SelectTaxo <- renderUI({
    
    data = dataInput()
    if(!is.null(data$taxo) && nrow(data$taxo)>0)
    { 
      tmp = colnames(data$taxo)
      selectInput("TaxoSelect",h6(strong("Select the taxonomy")),c("...",tmp))
    }
    else selectInput("TaxoSelect",h6(strong("Select the taxonomy")),c("..."))

  })
  


  
  # Infobox taxo
  output$InfoTaxo <- renderInfoBox({
  
    taxo = input$TaxoSelect
    print(taxo)
    if(taxo!="...") 
    {
      infoBox(h6(strong("Taxonomy")), subtitle = h6(taxo), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
    }
    else infoBox(h6(strong("Taxonomy")), subtitle = h6("Select the taxonomy for the analysis") ,color = "orange",width=NULL,fill=TRUE, icon = icon("warning"))
  })



#####################################################
##
##                Diagnostic plots
##
#####################################################


  
  output$VarIntBarPlot <- renderUI({
    
    int = input$InterestVar
    if(length(int)>=2) intSel = int[c(1,2)]
    else intSel = int[1]
    
    selectizeInput("VarInt",h6(strong("Select the variables of interest (max 2)")),int, selected = intSel,multiple = TRUE,options = list(maxItems = 2))
    
  })
  
  
  output$PlotDiag <- renderPlot({
    
    resDiff = ResDiffAnal()
    Plot_diag(input,resDiff)
  })


  SizeFactor_table <-reactive({ 
    res = ResDiffAnal()
    return(t(data.frame(Factor=sizeFactors(res$dds))))
    
  })

  output$SizeFactTable <- renderDataTable(
    SizeFactor_table(),
    options = list(scrollX=TRUE,searching = FALSE
  ))

#####################################################
##
##                EXPORT DIAG GRAPH
##
#####################################################

  ## PDF  
  output$exportPDFdiag <- downloadHandler(
    filename <- function() { paste(input$DiagPlot,'meta16S.pdf',sep="_")},
    content <- function(file) {
      pdf(file, width = 6, height = 4)
      print(Plot_diag(input,ResDiffAnal()))
      dev.off()
    }
  )

  
  ## PNG
  output$exportPNGdiag <- downloadHandler(
    filename <- function() { paste(input$DiagPlot,'meta16S.png',sep="_") },
    content <- function(file) {
      png(file, width = 600, height = 400)
      print(Plot_diag(input,ResDiffAnal()))
      dev.off()
    }
  )




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