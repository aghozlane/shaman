if (!require("Rcpp")){
  install.packages("Rcpp")
}

if(!require(shinydashboard)){
  install.packages('shinydashboard')
  library(shinydashboard)
}

if(!require(rjson)){
  install.packages('rjson')
}

if(!require(devtools)){
  install.packages('devtools')
}

#library(plotly)

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

if (!require(d3heatmap)) {
  install.packages('d3heatmap')
  library(d3heatmap)
}

if (!require(biom)) {
  install.packages('biom')
  library(biom)
}

if (!require(devtools)) {
  install.packages('devtools')
  library(devtools)
}

if (!require(scatterD3)) {
  install.packages('scatterD3')
  library(scatterD3)
}

if (!require(rNVD3)) {
  library(devtools)
  install_github('rNVD3', 'ramnathv')
  library(rNVD3)
}

if (!require(genefilter)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("genefilter")
  library(genefilter)
}

if (!require(googleVis)) {
  install.packages('googleVis')
  library(googleVis)
}

library(shinyjs)

# Allow to upload 50M files
options(shiny.maxRequestSize=50*1024^2) 
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
  namesfile = tempfile(pattern = "BaseContrast", tmpdir = tempdir(), fileext = "")
  file.create(namesfile,showWarnings=FALSE)

  ## Popup messages
  observe(if(input$AddRegScatter) info("By adding the regression line, you will lose interactivity."))

  ## Counts file
  dataInputCounts <-reactive({ 
    
    inFile <- input$fileCounts
  
    if (is.null(inFile)) return(NULL)

    data = read.csv(inFile$datapath,sep="\t",header=TRUE)

    ## Rownames
    if(!TRUE%in%duplicated(data[,1])) rownames(data)=data[,1];data=data[,-1]

    return(as.data.frame(data))
  })
  


  ## Taxo File
  dataInputTaxo <-reactive({ 
    
    inFile <- input$fileTaxo
    
    if (is.null(inFile)) return(NULL)
    
    if(input$TypeTaxo=="Table") 
    {
      data = read.csv(inFile$datapath,sep="\t",header=TRUE)
      
      ## Rownames
      if(!TRUE%in%duplicated(data[,1])){ 
        DataNames=data[,1]
        colNames=colnames(data)[-1]
        data=as.matrix(data[,-1])
        rownames(data)=DataNames
        colnames(data) = colNames
        }
    }
    
    if(input$TypeTaxo=="RDP") 
    {
      data = read_rdp(inFile$datapath,input$RDP_th)
    }
    
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
    data = read_biom(inFile$datapath)
    
    return(data)
  })



  ## Input data
  dataInput <-reactive({ 
    
    data = NULL
    check = NULL
    percent = NULL
    
    if(input$FileFormat=="fileCounts")
    {
      Counts = dataInputCounts()
      Taxo = dataInputTaxo()
      if(!is.null(Counts) && !is.null(Taxo))
      { 
        tmp = GetDataFromCT(Counts,Taxo)
        data = list(counts=tmp$counts,taxo=tmp$taxo)
        check = list(CheckCounts=tmp$CheckCounts,CheckTaxo=tmp$CheckTaxo,CheckPercent=tmp$CheckPercent)
        percent = tmp$Percent
      }    
    }
    
    if(input$FileFormat=="fileBiom")
    {
      tmpBIOM = dataInputBiom()
      if(!is.null(tmpBIOM))
      {
        tmp = GetDataFromBIOM(tmpBIOM)
        data = list(counts=tmp$counts,taxo=tmp$taxo)
        check = list(CheckCounts=tmp$CheckCounts,CheckTaxo=tmp$CheckTaxo,CheckPercent=tmp$CheckPercent)
        percent = tmp$Percent
      }
    }
    
    return(list(data=data,check=check,percent=percent))
  })
  
  

  ## Merge counts data
  dataMergeCounts <-reactive({ 
    input$RunDESeq
    
    counts = NULL
    CheckTarget = FALSE
    normFactors = NULL
    CT_noNorm = NULL
    #labeled= NULL
    data = dataInput()$data
    target = isolate(dataInputTarget()$target)
    taxo = isolate(input$TaxoSelect)
    
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0 && !is.null(taxo) && taxo!="..." && !is.null(target)) 
    {
      design = GetDesign(input)
      tmp = isolate(GetCountsMerge(input,data,taxo,target,design))
      counts = tmp$counts
      CheckTarget = tmp$CheckTarget
      #target = tmp$target
      #labeled = tmp$labeled
      normFactors = tmp$normFactors
      CT_noNorm = tmp$CT_noNorm
    }
    return(list(counts=counts,CheckTarget=CheckTarget,normFactors=normFactors,CT_noNorm=CT_noNorm))
    #return(list(counts=counts,target=target,labeled=labeled,normFactors=normFactors,CT_noNorm=CT_noNorm))
  })


  # Infobox Error counts
  output$InfoErrorCounts <- renderInfoBox({
    
    tmp = dataInput()
    data = tmp$data
    check = tmp$check
    cond = (!is.null(data$counts) && nrow(data$counts)>0)
    res =infoBox(h6(strong("Count table")), subtitle = h6("Load the count table") ,color = "light-blue",width=NULL,fill=TRUE, icon = icon("upload"))
    
    if(cond)
    {
      if(!is.null(check$CheckCounts$Warning)) res = infoBox(h6(strong("Count table")), subtitle = h6(check$CheckCounts$Warning), icon = icon("warning"),color = "orange",width=NULL,fill=TRUE)
      if(!is.null(check$CheckCounts$Error)) res = infoBox(h6(strong("Count table")), subtitle = h6(check$CheckCounts$Error), icon = icon("thumbs-o-down"),color = "red",width=NULL,fill=TRUE)
      if(is.null(check$CheckCounts$Error) && is.null(check$CheckCounts$Warning)) res = infoBox(h6(strong("Count table")), subtitle = h6(paste("Format of the count table seems to be OK")), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
    }
    
    return(res)
  })

  # Infobox Error counts
  output$InfoErrorTaxo <- renderInfoBox({
    
    tmp = dataInput()
    data = tmp$data
    check = tmp$check
    cond = (!is.null(data$taxo) && nrow(data$taxo)>0)
    res =infoBox(h6(strong("Taxonomy table")), subtitle = h6("Load the taxonomy table") ,color = "light-blue",width=NULL,fill=TRUE, icon = icon("upload"))
    
    if(cond)
    {
      if(!is.null(check$CheckTaxo$Warning)) res = infoBox(h6(strong("Taxonomy table")), subtitle = h6(check$CheckTaxo$Warning), icon = icon("warning"),color = "orange",width=NULL,fill=TRUE)
      if(!is.null(check$CheckTaxo$Error)) res = infoBox(h6(strong("Taxonomy table")), subtitle = h6(check$CheckTaxo$Error), icon = icon("thumbs-o-down"),color = "red",width=NULL,fill=TRUE)
      if(is.null(check$CheckTaxo$Error) && is.null(check$CheckTaxo$Warning)) res = infoBox(h6(strong("Taxonomy table")), subtitle = h6(paste("Format of the taxonomy table seems to be OK")), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
    }
    
    return(res)
  })


  # Infobox Error counts
  output$valueErrorPercent <- renderInfoBox({
    
    tmp = dataInput()
    data = tmp$data
    check = tmp$check
    cond = (!is.null(data$counts) && nrow(data$counts)>0 && !is.null(data$taxo) && nrow(data$taxo)>0)
    res = valueBox(paste0(0, "%"),h6(strong("Annotated features")), color = "light-blue",width=NULL,icon = icon("list"))
    
    if(cond)
    {
      percent = round(100*tmp$percent,2)
      if(percent==0) res = valueBox(paste0(percent, "%"),h6(strong("Annotated features")), color = "red",width=NULL,icon = icon("list"))  
      if(percent!=0) res = valueBox(paste0(percent, "%"),h6(strong("Annotated features")), color = "green",width=NULL,icon = icon("list"))  
      
    }
    
    return(res)
  })


  
#####################################################
##
##                DYNAMIC MENU
##
#####################################################
  
  
  
  output$dymMenu <- renderMenu({
    
    tmp = dataInput()
    data = tmp$data
    check = tmp$check
    
    ## Check error in the counts and taxonomy table 
    CheckOK = (is.null(check$CheckCounts$Error) && is.null(check$CheckTaxo$Error)  && is.null(check$CheckPercent))

    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0 && CheckOK)
    {
      sidebarMenu(
        menuItem("Statistical analysis",
                 menuSubItem("Run differential analysis",tabName="RunDiff"),
                 menuSubItem("Diagnostic plots",tabName="DiagPlotTab"),
                 menuSubItem("Tables",tabName="TableDiff"),
                 icon = icon("bar-chart-o"), tabName = "AnaStat"
        ),
        menuItem("Visualization",icon = icon("area-chart"), tabName = "Visu"),
        menuItem("Perspective plots", icon = icon("pie-chart"), tabName = "Krona")

      )
     }
    else{
      sidebarMenu()
    }
    
  })
  
  

  #####################################################
  ##
  ##                DATA TABLE
  ##
  #####################################################
  
  ## Counts Table
  output$DataCounts <- renderDataTable(
    dataInput()$data$counts, 
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE
    ))
  
  ## Taxonomy table
  output$DataTaxo <- renderDataTable(
    dataInput()$data$taxo, 
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE
    ))

  
  ## Tab box for data visualisation
  output$TabBoxData <- renderUI({
    
    data=dataInput()$data
    
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0)
    {
      tabBox(width = NULL, selected = "Count table",
             tabPanel("Count table",dataTableOutput("DataCounts")),
             tabPanel("Taxonomy",dataTableOutput("DataTaxo")),
             tabPanel("Summary",h5(strong("Percentage of annotation")),htmlOutput("SummaryView"),
                      br(),h5(strong("Number of features by level:")),plotOutput("SummaryViewBarplot",width = 1200,height=500))
      )
    }
    
  })

  output$SummaryView <- renderGvis({
    data = dataInput()$data
    taxo = data$taxo
    counts = data$counts
    res = NULL
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0)
    {
      taxo = rbind(taxo,rep(NA,ncol(taxo)))
      tmpPercent = round(apply(is.na(taxo),2,table)["FALSE",]/nrow(counts)*100,2)

      #print(tmpPercent)
      df <- data.frame(Label = colnames(taxo),Value = tmpPercent)
    
      res = gvisGauge(df,options=list(min=0, max=100, greenFrom=80,
                                    greenTo=100, yellowFrom=60, yellowTo=80,
                                    redFrom=0, redTo=60, width=1200, height=300))
    }
    return(res)
  })
  
  
  output$SummaryViewBarplot <- renderPlot({
    data = dataInput()$data
    taxo = data$taxo
    counts = data$counts
    res = NULL
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0)
    {
      colors=rep(c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78", "#2ca02c","#98df8a","#d62728","#ff9896","#9467bd","#c5b0d5","#8c564b",
                   "#c49c94","#e377c2","#f7b6d2","#7f7f7f", "#c7c7c7","#bcbd22","#dbdb8d","#17becf","#9edae5"),ceiling(ncol(taxo)/20))
      tmp = apply(taxo,2,unique)
      nbfeatures = as.numeric(lapply(tmp,length)) -as.numeric(lapply(lapply(tmp,is.na),any))
      df <- data.frame(Label = colnames(taxo),Count = nbfeatures)
      df$Label = factor(df$Label,levels =colnames(taxo) )
      res = ggplot(df,aes(x=Label,y=Count,fill=Label))+geom_bar(stat="identity")
      res = res + theme_bw() + xlab("Taxonomy") + scale_fill_manual(values=colors) + guides(fill=FALSE)
    }
    return(res)
  })
  
  #####################################################
  ##
  ##                TARGET FILE
  ##
  #####################################################
  

  ## Load target file
  dataInputTarget <-reactive({ 
    
    inFile <- input$fileTarget
    counts = dataInput()$data$counts
    
    
    if (is.null(inFile)) return(NULL)
    
    
    data = read.csv(inFile$datapath,sep="\t",header=TRUE)
    rownames(data) <- as.character(data[, 1])
    ind = which(rownames(data)%in%colnames(counts))
    #ord = order(rownames(data))
    #data = data[ord,]
    ### A SUPPRIMER 
    #rownames(data) <- colnames(counts)
    
    # Percent annotated
#     print(ind)
#     print(colnames(counts))
#     print(rownames(data))
    labeled = length(ind)/length(colnames(counts))*100.0
    
    return(list(target=data[ind,], labeled=labeled))
  })



  ## Interest Variables
  output$SelectInterestVar <- renderUI({
        
    target=dataInputTarget()$target
    
    if(!is.null(target)) 
    {
      namesTarget = colnames(target)[2:ncol(target)]
      selectInput("InterestVar",h6(strong("Select the variables")),namesTarget,selected=namesTarget,multiple=TRUE)
    }
  
  })


  ## Interactions
  output$SelectInteraction2 <- renderUI({
        
    target = dataInputTarget()$target

    if(!is.null(target)) 
    {
      Interac = GetInteraction2(target)
      selectInput("Interaction2",h6(strong("Add interactions")),Interac,selected=NULL,multiple=TRUE)
    }
    
  })


  ## Reference radio buttons
  output$RefSelect <- renderUI({
    
    target = dataInputTarget()$target
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
    #target = dataInputTarget()
    #labeled = dataMergeCounts()$labeled
    labeled = dataInputTarget()$labeled
    if(!is.null(labeled)) 
    {
      #### Ajout fontion check target
      #infoBox(h6(strong("Target file")), subtitle = h6("Your target file is OK"), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
      labeled = round(labeled,2)
      if(labeled>0) res = valueBox(paste0(labeled, "%"),h6(strong("Labeled features")), color = "green",width=NULL,icon = icon("list"))
      else res = valueBox(paste0(labeled, "%"),h6(strong("Labeled features")), color = "red",width=NULL,icon = icon("list"))  
    }
    #else infoBox(h6(strong("Target file")), subtitle = h6("Label of the target file must correspond to count table column names") ,color = "light-blue",width=NULL,fill=TRUE, icon = icon("warning"))
    else res = valueBox(paste0(0, "%"),h6(strong("Labeled features")), color = "light-blue",width=NULL,icon = icon("list"))
    return(res)
  }
  )




  ## target table
  output$DataTarget <- renderDataTable(
  dataInputTarget()$target,
  options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                 pageLength = 10,scrollX=TRUE
  ))

  ## Counts table for the selected taxonomy level
  output$CountsMerge <- renderDataTable(
    round(counts(ResDiffAnal()$dds,normalized=TRUE)),
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE
    ))


  ## Box for merged counts
  output$BoxCountsMerge <- renderUI({
    
    counts = dataMergeCounts()$counts
    taxo = input$TaxoSelect
    
    if(!is.null(counts) && taxo != "...")
    {
      box(title=paste("Count table (",taxo,")",sep=""),width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = TRUE,
          dataTableOutput("CountsMerge"),
          downloadButton('ExportCounts', 'Export normalised counts'),
          downloadButton('ExportRelative', 'Export relative abundance')
      )  
    }
    
  })

  ## Export in .csv
  output$ExportCounts <- downloadHandler(
    filename = function() { 'NormCounts.csv' },
    content = function(file){write.csv(dataMergeCounts()$counts, file)}
  )

  ## Export in .csv
  output$ExportRelative <- downloadHandler(
    filename = function() { 'RelativeAb.csv' },
    content = function(file){write.csv(dataMergeCounts()$counts/colSums(dataMergeCounts()$counts), file)}
  )


  ## Box for target visualisation
  output$BoxTarget <- renderUI({
    
    target = dataInputTarget()$target

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
    
    for(i in 1:length(names)){Contrast[[i]] = numericInput(names[i],names[i],0,step=1,min=-1,max=1)}
  
    return(Contrast)
    
  })


  output$ContrastOverview <- renderPrint({
    
    resDiff = ResDiffAnal()
    dds = resDiff$dds
    names = resultsNames(dds)
    
    cont = input$ContrastList
    filesize = file.info(namesfile)[,"size"]
    if(is.na(filesize)){filesize=0}
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
    filesize = file.info(namesfile)[,"size"]
    if(is.na(filesize)){filesize=0}
    if(filesize!=0) tmp = read.table(namesfile,header=TRUE)
    Contrast = colnames(as.matrix(tmp))
    updateSelectInput(session, "ContrastList","Contrasts",Contrast)
    updateSelectInput(session, "ContrastList_table","Contrasts",Contrast)
    updateSelectInput(session, "ContrastList_table_Visu","For which contrasts",Contrast)
    updateSelectInput(session, "ContrastList_table_FC","Contrasts (Min = 2)",Contrast)
  })

  ## Add contrast 
  observeEvent(input$AddContrast,{  
    
    AddCont()
    
  },priority=1)


  ## Add contrast function
  AddContEasy <-eventReactive(input$AddContrastEasy,{
    
    resDiff = ResDiffAnal()
    dds = resDiff$dds
    names = resultsNames(dds)
#     
#     BaseContrast(input,names,namesfile)
#     tmp = read.table(namesfile,header=TRUE)
#     Contrast = colnames(as.matrix(tmp))
    updateSelectInput(session, "ContrastList","Contrasts",Contrast)
    updateSelectInput(session, "ContrastList_table","Contrasts",Contrast)
    updateSelectInput(session, "ContrastList_table_Visu","For which contrasts",Contrast)
    updateSelectInput(session, "ContrastList_table_FC","Contrasts (Min = 2)",Contrast)
  })
  
  ## Add contrast 
  observeEvent(input$AddContrastEasy,{  
    
    AddContEasy()
    
  })


  AddContFromFile <-eventReactive(input$fileContrast,{ 
    
    res = ReadContrastFile()
    createdCont = NULL
    filesize = file.info(namesfile)[,"size"]
    if(is.na(filesize)){filesize=0}
    if(filesize!=0){ createdCont = read.table(namesfile,header=TRUE) }
      
    if(!is.null(res))
    { 
      if(!is.null(createdCont)) res = cbind(res,createdCont)
      updateSelectInput(session, "ContrastList","Contrasts",colnames(res))
      updateSelectInput(session, "ContrastList_table","Contrasts",colnames(res))
      updateSelectInput(session, "ContrastList_table_Visu","For which contrasts",colnames(res))
      updateSelectInput(session, "ContrastList_table_FC","Contrasts (Min = 2)",colnames(res))
      write.table(res,namesfile,row.names=FALSE)
    }
  })
 
 observeEvent(input$fileContrast,{ 
    
   AddContFromFile()
  },priority=1)


  ## Remove contrast function
  RemoveCont <-eventReactive(input$RemoveContrast,{
    
    ## get the size of the contrast base file
    filesize = file.info(namesfile)[,"size"]
    if(is.na(filesize)){filesize=0}
    if(filesize!=0)
    { 
      tmp = read.table(namesfile,header=TRUE)
      ind = which(colnames(tmp)%in%input$ContrastList)
      matKept = as.matrix(tmp[,-ind])
      ContrastKept = colnames(tmp)[-ind]
      
      if(ncol(matKept)>0) write.table(matKept,namesfile,row.names=FALSE,col.names = ContrastKept)
      else file.create(namesfile,showWarnings=FALSE)
      updateSelectInput(session, "ContrastList","Contrasts",ContrastKept)
      updateSelectInput(session, "ContrastList_table","Contrasts",ContrastKept)
      updateSelectInput(session, "ContrastList_table_Visu","For which contrasts",ContrastKept)
      updateSelectInput(session, "ContrastList_table_FC","Contrasts (Min = 2)",ContrastKept)
    }
  })
  
  
  
  ## Remove contrast
  observeEvent(input$RemoveContrast,{  
    
    RemoveCont()
    
  },priority=1)

  ## Remove all contrasts
  RemoveAllCont <-eventReactive(input$RunDESeq,{
    
      file.create(namesfile,showWarnings=FALSE)
      updateSelectInput(session, "ContrastList","Contrasts",NULL)
      updateSelectInput(session, "ContrastList_table","Contrasts",NULL)
      updateSelectInput(session, "ContrastList_table_Visu","For which contrasts",NULL)
      updateSelectInput(session, "ContrastList_table_FC","Contrasts (Min = 2)",NULL)
  })

  ## Remove all contrast
  observeEvent(input$RunDESeq,{  
    
    RemoveAllCont()
    
  })

# Infobox Contrast
output$InfoContrast <- renderInfoBox({
  input$RunDESeq
  input$AddContrast
  input$RemoveContrast
  input$fileContrast
  resDiff = ResDiffAnal()
  res=NULL
  if(!is.null(resDiff)){
  
    res = infoBox("Contrasts", subtitle = h6("At least one contrast (non null) must be defined"), icon = icon("warning"),color = "light-blue",width=NULL,fill=TRUE)
    test = FALSE
    filesize = isolate(file.info(namesfile)[,"size"])
    
    if(is.na(filesize)){filesize=0}
    if(filesize!=0) 
    {
      tmp = read.table(namesfile,header=TRUE)
      if(any(as.vector(tmp)!=0)) test = TRUE
    }
    
    if(test) res = infoBox("Contrasts", subtitle = h6("Contrasts OK"), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
  }
    return(res)
  })


  
  output$contrastBox <- renderUI({
    
    resDiff = ResDiffAnal()
    int = input$Interaction2
    target = resDiff$target
    
    InterVar = input$InterestVar
    ModInterestAll = unique(target[,InterVar])
    test = c()
    for(i in 1:length(InterVar)){ test =c(test,input$Select1_contrast%in%target[,InterVar[i]]) }
    
    ModInterestCond = unique(target[,which(test)])
    #alltmp = c(InterVar,Interaction)
    
    if(!is.null(resDiff))
    { 
      box(title="Contrasts",width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = FALSE,
          fluidRow(
            column(width=3,selectInput("Select1_contrast",label=h6(strong("Compare")),ModInterestAll)),
            column(width=3,selectInput("Select2_contrast",label=h6(strong("To")),ModInterestCond)),
            if(length(int)>=1) column(width=3,selectInput("Select3_contrast",label=h6(strong("For")),c("All","WT","Delta"))),
            column(width=3,br(),br(),actionButton("AddContrastEasy","Add",icon = icon("plus")))
          )
      )
    }
    
  })


  output$contrastBoxAdvanced <- renderUI({
    
    resDiff = ResDiffAnal()
    
    if(!is.null(resDiff))
    { 
      box(title="Contrasts (advanced user)",width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = FALSE,
          fluidRow(
            column(width=12,
                   fileInput('fileContrast', h6(strong('Select a file of contrasts')),width="60%")
            ),
            hr(),
            column(width=12,h6(strong("Define contrasts by yourself"))),
            column(width=6,textInput("ContrastName",label = NULL,value = "Contrast name")),
            column(width=6,actionButton("AddContrast","Add contrast",icon = icon("plus")))
          ),
          fluidRow(column(width=12,uiOutput("contrastMat")))
      )
    }
    
  })


  output$contrastDefined <- renderUI({
    resDiff = ResDiffAnal()
    
    if(!is.null(resDiff))
    { 
      box(title="Defined contrasts",width = NULL, status = "primary", solidHeader = TRUE,collapsible = FALSE,collapsed = FALSE,
          fluidRow(
            column(width=11,
                   selectInput("ContrastList","Contrasts","",multiple=TRUE,size=4,selectize=FALSE,width = '100%'),
                   fluidRow(
                      column(width=6,actionButton("RemoveContrast","Remove",icon = icon("remove"))),
                      column(width=6,downloadButton("exportContrast", "Export"))
                   ),
                   htmlOutput("ContrastOverview")
            )
            
          )
      )
    }
  })


  ReadContrastFile <-reactive({ 
    
    inFile <- input$fileContrast
    
    if (is.null(inFile)) return(NULL)
    
    res = read.csv(inFile$datapath,sep=" ",header=TRUE)
    
    return(res)
  })
  

  output$exportContrast <- downloadHandler(
    filename <- function() {"Contrasts.txt"},
    content <- function(file) {
      file.copy(namesfile,file)
    }
  )




#####################################################
##
##                DESEQ2 run
##
#####################################################



  # Infobox Contrast
  output$InfoDESeq <- renderInfoBox({
      input$RunDESeq
      
      resDiff = ResDiffAnal()
      if(!is.null(resDiff)){
        infoBox(h6(strong("Statistical analysis")), subtitle = h6("Differential analysis is done !"), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
      }

    
  })
  

  ## Get the results from DESeq2
  ResDiffAnal <-eventReactive(input$RunDESeq,{
    
    counts = dataMergeCounts()$counts
    CT_noNorm = dataMergeCounts()$CT_noNorm
    normFactors = dataMergeCounts()$normFactors
    ## HEEEERRREE
    target = dataInputTarget()$target
    #print(target)
    design = GetDesign(input)
   
    Get_dds_object(input,counts,target,design,normFactors,CT_noNorm)

    
  })


  ## Run DESeq2 via RunDESeq button
  observeEvent(input$RunDESeq,{
    withProgress(message="Analysis in progress...",ResDiffAnal())
  })

  
#####################################################
##
##                Taxonomy
##
#####################################################
  
  
  # Infobox Contrast
  output$SelectTaxo <- renderUI({
    
    data = dataInput()$data
    if(!is.null(data$taxo) && nrow(data$taxo)>0)
    { 
      tmp = colnames(data$taxo)
      selectInput("TaxoSelect",h6(strong("Select the taxonomy level")),c("...",tmp,"OTU"))
    }
    else selectInput("TaxoSelect",h6(strong("Select the taxonomy level")),c("..."))

  })
  


  
  # Infobox taxo
  output$InfoTaxo <- renderInfoBox({
  
    data = dataInput()$data
    taxo = input$TaxoSelect
    
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0 && !is.null(taxo) && taxo!="...") 
    {
      counts = dataMergeCounts()$counts
      nfeature = nrow(counts)
      infoBox(h6(strong("Taxonomy")), subtitle = h6(paste(taxo, ", nb features: ",nfeature,sep="")), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
    }
    else infoBox(h6(strong("Taxonomy")), subtitle = h6("Select the taxonomy for the analysis") ,color = "light-blue",width=NULL,fill=TRUE, icon = icon("warning"))
  })



#####################################################
##
##                Diagnostic plots
##
#####################################################


  
  output$VarIntDiag <- renderUI({
    
    target=dataInputTarget()$target
    
    if(!is.null(target)) 
    {
      namesTarget = colnames(target)[2:ncol(target)]
      selectizeInput("VarInt",h6(strong("Select the variables of interest (max 2)")),namesTarget, selected = namesTarget[1],multiple = TRUE,options = list(maxItems = 2))
    }
    
  })
  
  
  output$PlotDiag <- renderPlot({
    input$RunDESeq
    
    resDiff = isolate(ResDiffAnal())
    Plot_diag(input,resDiff)
  },height = reactive(input$heightDiag))



  output$PlotpcoaEigen <- renderPlot({
    
    resDiff = ResDiffAnal()
    Plot_diag_pcoaEigen(input,resDiff)
  },height = 400)

  output$PlotEigen <- renderPlot({
    
    resDiff = ResDiffAnal()
    Plot_diag_Eigen(input,resDiff)
  },height =400)

  SizeFactor_table <-reactive({ 
    res = ResDiffAnal()
    return(t(data.frame(Factor=res$normFactors)))
    
  })

  output$SizeFactTable <- renderDataTable(
    SizeFactor_table(),
    options = list(scrollX=TRUE,searching = FALSE
  ))


  ## Select Modality DiagPlot

  output$ModMat <- renderUI({
    
    VarInt = input$VarInt
    target = dataInputTarget()$target
    
    Mod = list()
    
    for(i in 1:length(VarInt)){
      value = as.character(unique(as.factor(target[,VarInt[i]])))
      Mod[[i]] = selectizeInput(paste("Mod",VarInt[i],sep=""),VarInt[i],value,selected=value, multiple = TRUE)
    }
  
    return(Mod)
    
  })
  
  
  ## Select Modality VisuPlot
  
  output$ModVisu <- renderUI({
    Mod = NULL
    
    VisuVarInt = input$VisuVarInt
    target = dataInputTarget()$target
    if(!is.null(VisuVarInt))
    {
      Mod = list()
      
      for(i in 1:length(VisuVarInt)){
        value = as.character(unique(as.factor(target[,VisuVarInt[i]])))
        Mod[[i]] = selectizeInput(paste("ModVisu",VisuVarInt[i],sep=""),VisuVarInt[i],value,selected=value, multiple = TRUE)
      }
    }
    return(Mod)
    
  })
  
  

#####################################################
##
##                EXPORT DIAG GRAPH
##
#####################################################

  #### Export Diag
  output$exportdiag <- downloadHandler(
    filename <- function() { paste(input$DiagPlot,paste('SHAMAN',input$Exp_format,sep="."),sep="_") },
    content <- function(file) {
      if(input$Exp_format=="png") png(file, width = input$widthDiagExport, height = input$heightDiagExport)
      if(input$Exp_format=="pdf") pdf(file, width = input$widthDiagExport/96, height = input$heightDiagExport/96)
      if(input$Exp_format=="eps") postscript(file, width = input$widthDiagExport/96, height = input$heightDiagExport/96)
      if(input$Exp_format=="svg") svg(file, width = input$widthDiagExport/96, height = input$heightDiagExport/96)
      
      print(Plot_diag(input,ResDiffAnal()))
      dev.off()
    }
  )

#####################################################
##
##                EXPORT VISU GRAPH
##
#####################################################

  ## PDF  
## PDF  
output$exportPDFVisu <- downloadHandler(
  filename <- function() { paste("test",'SHAMAN.ps',sep="_")},
  content <- function(file) {
    resDiff = ResDiffAnal()
    filesize = file.info(namesfile)[,"size"]
    if(is.na(filesize)){filesize=0}
    if(filesize!=0)
    {
      BaseContrast = read.table(namesfile,header=TRUE)
      #ggsave(filename = filename, Plot_Visu_Heatmap_FC(input,BaseContrast,resDiff),width = input$widthHeat, height = input$heightHeat)
      postscript(file, width = input$widthHeat, height = input$heightHeat)
      if(input$HeatMapType=="Counts")  Plot_Visu_Heatmap(input,resDiff)
      if(input$HeatMapType=="Log2FC")     Plot_Visu_Heatmap_FC(input,BaseContrast,resDiff)
    dev.off()
    }
  }
)
  
  
  ## PNG
  output$exportPNGVisu <- downloadHandler(
    filename <- function() { paste("test",'SHAMAN.png',sep="_") },
    content <- function(file) {
      png(file, width = 600, height = 400)
      Plot_Visu(input,ResDiffAnal())
      dev.off()
    }
  )

#### Export Visu
output$exportVisu <- downloadHandler(
  filename <- function() { paste(input$PlotVisuSelect,paste('SHAMAN',input$Exp_format_Visu,sep="."),sep="_") },
  content <- function(file) {

      taxo = input$TaxoSelect

      if(input$Exp_format_Visu=="png") png(file, width = input$widthVisuExport, height = input$heightVisuExport)
      if(input$Exp_format_Visu=="pdf") pdf(file, width = input$widthVisuExport/96, height = input$heightVisuExport/96)
      if(input$Exp_format_Visu=="eps") postscript(file, width = input$widthVisuExport/96, height = input$heightVisuExport/96,paper="special")
      if(input$Exp_format_Visu=="svg") svg(file, width = input$widthVisuExport/96, height = input$heightVisuExport/96)

      if(input$PlotVisuSelect=="Barplot") print(Plot_Visu_Barplot(input,ResDiffAnal())$gg)
      if(input$PlotVisuSelect=="Heatmap"){
        if(input$HeatMapType=="Counts") print(Plot_Visu_Heatmap(input,ResDiffAnal(),export=TRUE))
        if(input$HeatMapType=="Log2FC") {      
          BaseContrast = read.table(namesfile,header=TRUE)
          filesize = file.info(namesfile)[,"size"]
          if(is.na(filesize)){filesize=0}
          if(filesize!=0) print(Plot_Visu_Heatmap_FC(input,BaseContrast,ResDiffAnal(),export=TRUE))
        }
      } 
    
      if(input$PlotVisuSelect=="Boxplot") print(Plot_Visu_Boxplot(input,ResDiffAnal(),alpha=ifelse(input$Exp_format_Visu=="eps",1,0.7)))
      if(input$PlotVisuSelect=="Scatterplot") print(Plot_Visu_Scatterplot(input,ResDiffAnal(),export=TRUE,lmEst = FALSE))
      if(input$PlotVisuSelect=="Diversity") print(Plot_Visu_Diversity(input,ResDiffAnal(),type="point"))
      if(input$PlotVisuSelect=="Rarefaction") print( Plot_Visu_Rarefaction(input,ResDiffAnal(),ranges$x,ranges$y,ylab=taxo))
      dev.off()
    
  }
)



#####################################################
##
##                DIFF TABLES
##
#####################################################


#   output$ContrastListTable <- renderUI({
#     
#     filesize = file.info(namesfile)[,"size"]
#     if(filesize!=0)
#     { 
#       tmp = read.table(namesfile,header=TRUE)
#       cont = colnames(tmp)
#     }  
#     selectInput("ContrastList_table",h6(strong("Contrast list")),cont, multiple = FALSE)
#     
#   })


  output$ContrastOverviewTable <- renderPrint({
    
    resDiff = ResDiffAnal()
    dds = resDiff$dds
    names = resultsNames(dds)
    
    cont = input$ContrastList_table
    filesize = file.info(namesfile)[,"size"]
    if(is.na(filesize)){filesize=0}
    if(filesize!=0)
    { 
      ContrastBase = read.table(namesfile,header=TRUE)
      ind = which(colnames(ContrastBase)%in%cont)
      div(HTML(PrintContrasts(names,sapply(ContrastBase[,ind],as.numeric),cont)))
    }
  })




  ## Get the diff table
  dataDiff <-reactive({ 
    input$AddContrast
    input$RemoveContrast
    input$fileContrast
    input$RunDESeq
    
    resDiff = ResDiffAnal()
    filesize = file.info(namesfile)[,"size"]
    res=NULL 
    if(is.na(filesize)){filesize=0}
    if(filesize!=0)
    { 
      BaseContrast = read.table(namesfile,header=TRUE)
      res = TableDiff_print(input,BaseContrast,resDiff, info = NULL) 
    } 
    
    return(res)
  })



  ## Complete diff table
  output$DataDiffcomplete <- renderDataTable(
    dataDiff()$complete,
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE
    ))
  
  ## Up diff table
  output$DataDiffup <- renderDataTable(
    dataDiff()$up,
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE
    ))
  
  ## Down diff table
  output$DataDiffdown <- renderDataTable(
    dataDiff()$down,
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE
    ))


  ## TabBox for diff table
  output$TabBoxDataDiff <- renderUI({
    
    data = dataDiff()
    
    if(!is.null(data))
    {
      
      tabBox(width = NULL, selected = "Complete",
             tabPanel("Complete",dataTableOutput("DataDiffcomplete")),
             tabPanel("Up",dataTableOutput("DataDiffup")),
             tabPanel("Down",dataTableOutput("DataDiffdown"))  
      )
    }
    
  })


#####################
###
###
###################


#### Export diff table
output$exportDiffTable <- downloadHandler(
  filename = function() {paste(input$WhichExportTable,input$ContrastList_table,'SHAMAN.csv',sep="_")},
  content = function(file){
    switch(input$WhichExportTable,
           "Complete" = write.csv(dataDiff()$complete, file),
           "Up" =  write.csv(dataDiff()$up, file),
           "Down" =  write.csv(dataDiff()$down, file)
    )
  }
)



output$ExportTableButton <- renderUI({
  
  res = NULL
  table = dataDiff()$complete
  if(nrow(table)>0) res = downloadButton('exportDiffTable', 'Export table')
  
  return(res)
})








  ## Run button

output$RunButton <- renderUI({
  
  res = NULL
  target = dataInputTarget()$target
  taxo = input$TaxoSelect
  if(!is.null(target) && taxo!="...") res = actionButton("RunDESeq",strong("Run analysis"),icon = icon("caret-right"))
  
  return(res)
})




#####################################################
##
##                VISUALISATION
##
#####################################################


  output$PlotVisuBar <- renderChart({
    resDiff = ResDiffAnal()
    res = NULL
    if(!is.null(resDiff$dds) && length(input$VisuVarInt)>=1) tmp = Plot_Visu_Barplot(input,resDiff)
    if(!is.null(tmp)) res = tmp$plotd3
    return(res)
    })

# output$PlotVisu <- renderPlotly({
#   resDiff = ResDiffAnal()
#   if(!is.null(resDiff$dds)) Plot_Visu_Barplot(input,resDiff)
# })
  
  output$heatmap <- renderD3heatmap({
    resDiff = ResDiffAnal()
    filesize = file.info(namesfile)[,"size"]
    if(is.na(filesize)){filesize=0}
    resplot = NULL
    if(!is.null(resDiff$dds))
    { 
      if(input$HeatMapType=="Counts") resplot = withProgress(message="Loading...",Plot_Visu_Heatmap(input,resDiff))
      if(input$HeatMapType=="Log2FC" && filesize!=0)
      { 
        BaseContrast = read.table(namesfile,header=TRUE)
        resplot = withProgress(message="Loading...",Plot_Visu_Heatmap_FC(input,BaseContrast,resDiff))
      }
    }
    return(resplot)
  },env=new.env())
  
  
  output$ScatterplotD3 <- renderScatterD3({
    resDiff = ResDiffAnal()
    if(!is.null(resDiff$dds)) withProgress(message="Loading...",Plot_Visu_Scatterplot(input,resDiff))
  })

  output$Scatterplotgg <- renderPlot({
    resDiff = ResDiffAnal()
    if(!is.null(resDiff$dds)) withProgress(message="Loading...",Plot_Visu_Scatterplot(input,resDiff,lmEst=FALSE))
  })
  
  ## Regression coefficients Table
  output$lmRegScatter <- renderDataTable(
    Plot_Visu_Scatterplot(input,ResDiffAnal(),lmEst=TRUE)$regCoef, 
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE
    ))
  
  ## Correlation coefficients Table
  output$CorTable <- renderDataTable(
    Plot_Visu_Scatterplot(input,ResDiffAnal(),CorEst=TRUE)$cor.est, 
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE
    ))
  
  output$lmEquation <- renderPrint({ 
    res = Plot_Visu_Scatterplot(input,ResDiffAnal(),lmEst=TRUE)
    coef = res$regCoef
    Rsq = res$Rsq
    
    div(HTML(paste(h4(strong("Linear equation: ")),
               "y =", round(coef[2,1],2),'x ',ifelse(coef[1,1]>=0,"+",""), round(coef[1,1],2),'<br/>','<br/>',
               h4(strong("Adjusted R squared:")),round(Rsq,5)*100," %")))
  })
  
  output$Boxplot <- renderPlot({
    resDiff = ResDiffAnal()
    if(!is.null(resDiff$dds)) withProgress(message="Loading...",Plot_Visu_Boxplot(input,resDiff))
  },height=reactive(input$heightVisu))


  output$DiversityPlot <- renderPlot({
    resDiff = ResDiffAnal()
    if(!is.null(resDiff$dds)) withProgress(message="Loading...",Plot_Visu_Diversity(input,resDiff,type="point"))
  })
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  output$RarefactionPlot <- renderPlot({
    resDiff = ResDiffAnal()
    taxo = input$TaxoSelect
    if(!is.null(resDiff)) withProgress(message="Loading...",Plot_Visu_Rarefaction(input,resDiff,ranges$x,ranges$y,ylab=taxo))
  }, height = reactive(input$heightVisu))
  
  observeEvent(input$RarefactionPlot_dblclick, {
    brush <- input$RarefactionPlot_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })


# 
#   output$SelectVarBoxDiv <- renderUI({
#     
#     selectVar = input$VisuVarInt
#     
#     if(!is.null(selectVar)) 
#     {
#       selectInput("VarBoxDiv", h6(strong("By")),selectVar)
#     }
#     
#   })
#   
  output$plotVisu <- renderUI({
    
    res=NULL
    if(input$PlotVisuSelect=="Barplot") res =  showOutput("PlotVisuBar")
    if(input$PlotVisuSelect=="Heatmap") res =  d3heatmapOutput("heatmap", height = input$heightVisu+10)
    if(input$PlotVisuSelect=="Boxplot") res = plotOutput("Boxplot", height = input$heightVisu+10)
    if(input$PlotVisuSelect=="Scatterplot" && !input$AddRegScatter) res = scatterD3Output("ScatterplotD3", height = input$heightVisu+10)
    if(input$PlotVisuSelect=="Scatterplot" && input$AddRegScatter) res = plotOutput("Scatterplotgg", height = input$heightVisu+10)
    if(input$PlotVisuSelect=="Diversity") res =  plotOutput("DiversityPlot", height = input$heightVisu+10)
    if(input$PlotVisuSelect=="Rarefaction") res = plotOutput("RarefactionPlot",dblclick = "RarefactionPlot_dblclick",brush = brushOpts(id = "RarefactionPlot_brush",resetOnNew = TRUE), height = input$heightVisu+10)
    
    return(res)
  })


  output$ColBoxplot <- renderUI({
    
    VarInt = input$VisuVarInt
    res = NULL
    if(length(VarInt)>1) res = selectizeInput("BoxColorBy",h6(strong(paste("Color by"))),VarInt, selected = VarInt,multiple = TRUE)
    return(res)
  })
  
  
#   output$DiversityBoxPlot <- renderPlot({
#     resDiff = ResDiffAnal()
#     if(!is.null(resDiff$dds)) Plot_Visu_Diversity(input,resDiff,type="box")
#   })




  output$TaxoToPlotVisu <- renderUI({
    
    data = dataInput()$data 
    taxo = input$TaxoSelect
    resDiff = ResDiffAnal()
    res = NULL
    
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0 && !is.null(taxo) && taxo!="...") 
    {
      counts = dataMergeCounts()$counts
      sumTot = rowSums(counts)
      ord = order(sumTot,decreasing=TRUE)
      Available_taxo = rownames(counts)[ord]
      selTaxo = Available_taxo[1:min(12,length(Available_taxo))]
      
      if(input$SelectSpecifTaxo=='Most')  res = selectizeInput("selectTaxoPlot",h6(strong(paste("Select the",input$TaxoSelect, "to plot"))),Available_taxo, selected = selTaxo,multiple = TRUE)
      if(input$SelectSpecifTaxo=="Diff")
      {
        filesize = file.info(namesfile)[,"size"]
        if(is.na(filesize)){filesize=0}
        if(filesize!=0)
        { 
          BaseContrast = read.table(namesfile,header=TRUE)
          SelContrast = input$ContrastList_table_Visu
          padj = Get_log2FC_padj(input,BaseContrast,resDiff, info = NULL)$padj
          Feature_names = rownames(padj)
          if(ncol(as.matrix(padj))>1)
          { 
            cont = which(colnames(padj)%in%SelContrast)
            padj = padj[,cont] 
          }
          ind = which(padj<=as.numeric(input$AlphaVal))
          if(length(ind)>0) selTaxo = Feature_names[ind]
          else selTaxo = NULL
          res = selectizeInput("selectTaxoPlot",h6(strong(paste("Select the",input$TaxoSelect, "to plot"))),Available_taxo, selected = selTaxo,multiple = TRUE,options = list(minItems = 2))
        }    
      }
      if(input$SelectSpecifTaxo=="NoDiff")
      {
        filesize = file.info(namesfile)[,"size"]
        if(is.na(filesize)){filesize=0}
        if(filesize!=0)
        { 
          BaseContrast = read.table(namesfile,header=TRUE)
          SelContrast = input$ContrastList_table_Visu
          padj = Get_log2FC_padj(input,BaseContrast,resDiff, info = NULL)$padj
          Feature_names = rownames(padj)
          if(ncol(as.matrix(padj))>1)
          { 
            cont = which(colnames(padj)%in%SelContrast)
            padj = padj[,cont] 
          }
          ind = which(padj>as.numeric(input$AlphaVal))
          if(length(ind)>0) selTaxo = Feature_names[ind]
          else selTaxo = NULL
          res = selectizeInput("selectTaxoPlot",h6(strong(paste("Select the",input$TaxoSelect, "to plot"))),Available_taxo, selected = selTaxo,multiple = TRUE,options = list(minItems = 2))
        }    
      }
      if(input$SelectSpecifTaxo=="All") res = selectizeInput("selectTaxoPlot",h6(strong(paste("Select the",input$TaxoSelect, "to plot"))),Available_taxo, selected = Available_taxo,multiple = TRUE)
    }
    return(res)
})


  
  output$VarIntVisu <- renderUI({
    
#     int = input$InterestVar
#     if(length(int)>=2) intSel = int[c(1,2)]
#     else intSel = int[1]
        
      target=dataInputTarget()$target
      
      if(!is.null(target)) 
      {
        namesTarget = colnames(target)[2:ncol(target)]
        selectizeInput("VisuVarInt",h6(strong("Select the variables of interest")),namesTarget, selected = namesTarget[1],multiple = TRUE)
      }
    
  })
  
  
  output$VarIntVisuScatter <- renderUI({
    
    target=dataInputTarget()
    data = dataInput()$data 
    taxo = input$TaxoSelect
    resDiff = ResDiffAnal()
    res = list()
    namesTarget = colnames(target)[2:ncol(target)]
    
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0 && !is.null(taxo) && taxo!="..." && !is.null(target)) 
    {
      counts = dataMergeCounts()$counts
      
      ## Get numeric variables from target
      typesTarget = sapply(target,class)
      numInd = (typesTarget=="numeric")[2:ncol(target)]
      Available_x = sort(rownames(counts))
      if(any(numInd)) Available_x = c(Available_x,namesTarget[numInd])
      Available_y = Available_x
      
      res[[1]] = selectizeInput("Xscatter",h6(strong("X variable")),Available_x, selected = Available_x[1],multiple = FALSE)
      res[[2]] = selectizeInput("Yscatter",h6(strong("Y variable")),Available_y, selected = Available_x[2],multiple = FALSE)
      res[[3]] = selectizeInput("ColorBy",h6(strong("Color variable")),c("None"="None",namesTarget[!numInd]),multiple = FALSE)
      res[[4]] = selectizeInput("PchBy",h6(strong("Symbol variable")),c("None"="None",namesTarget[!numInd]),multiple = FALSE)
      res[[5]] = selectizeInput("PointSize",h6(strong("Point size according to")),c("None"="None",Available_x), selected = NULL,multiple = FALSE)
    }
    
    return(res)
    
  })
  
  
  #####################################################
  ##
  ##                KRONA
  ##
  #####################################################
  output$krona <- renderTable({
    data = dataInput()$data 
    taxo = input$TaxoSelect
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0 && !is.null(taxo) && taxo!="...") 
    {
      #print(data$counts)
    krona_table=tempfile(pattern = "krona", tmpdir = tempdir(), fileext = "")
    url=paste(krona_table, ".html", sep="")
    #system(paste("export PERL5LIB=/home/aghozlan/workspace/SHAMAN_App/KronaTools-2.6/lib:$PERL5LIB; /home/aghozlan/workspace/META10S_App/krona_bin/bin/ktImportText", krona_table))
    system(paste("ktImportText", krona_table))
    refs <- paste0("<a href='",  url, "' target='_blank'>krona</a>")
    
    data.frame(refs)
    }
  }, sanitize.text.function = function(x) x)

#   output$VarIntVisuBP <- renderUI({
#     
#     int = input$InterestVar
#     if(length(int)>=2) intSel = int[c(1,2)]
#     else intSel = int[1]
#     
#     selectizeInput("VisuVarIntBP",h6(strong("Select the variables of interest")),int, selected = intSel,multiple = TRUE,options = list(minItems = 1))
#     
#   })
# 
#   output$VarIntVisuHM <- renderUI({
#     
#     int = input$InterestVar
#     if(length(int)>=2) intSel = int[c(1,2)]
#     else intSel = int[1]
#     
#     selectizeInput("VisuVarIntHM",h6(strong("Select the variables of interest")),int, selected = intSel,multiple = TRUE,options = list(minItems = 1))
#     
#   })
# 
#   output$VarIntVisuBoxP <- renderUI({
#     
#     int = input$InterestVar
#     intSel = int[1]
#     
#     selectizeInput("VisuVarIntBoxP",h6(strong("X variable")),int, selected = intSel,multiple = TRUE)
#     
#   })
# 
#   output$VarIntVisuDiv <- renderUI({
#     
#     int = input$InterestVar
#     intSel = int[1]
#     
#     selectizeInput("VisuVarIntDiv",h6(strong("X variable")),int, selected = intSel,multiple = TRUE)
#     
#   })
  
#   output$DiversityGroupBy <- renderUI({
#     
#     int = input$InterestVar
#     intSel = int[1]
#     
#     selectizeInput("GroupBy",h6(strong("Group by")),int, selected = intSel)
#     
#   })


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