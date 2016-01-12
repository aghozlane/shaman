library(shinydashboard)
if (!require(rNVD3)) {
  install.packages('rNVD3')
  library(rNVD3)
}
library(plotly)
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
  #paste("/srv/shiny-server/sample-apps/meta16s/BaseContrast_",rand,".txt",sep="")
  file.create(namesfile,showWarnings=FALSE)

  #namesfile = "www/All_Contrast.txt"

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
      if(!TRUE%in%duplicated(data[,1])) rownames(data)=data[,1];data=data[,-1]
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
    
    counts = NULL
    CheckTarget = FALSE
    normFactors = NULL
    CT_noNorm = NULL
    data = dataInput()$data
    target = dataInputTarget()
    
    taxo = input$TaxoSelect
    
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0 && !is.null(taxo) && taxo!="..." && !is.null(target)) 
    {
      design = GetDesign(input)
      tmp = GetCountsMerge(input,data,taxo,target,design)
      counts = tmp$counts
      CheckTarget = tmp$CheckTarget
      normFactors = tmp$normFactors
      CT_noNorm = tmp$CT_noNorm
    }
    return(list(counts=counts,CheckTarget=CheckTarget,normFactors=normFactors,CT_noNorm=CT_noNorm))
  })


  # Infobox Error counts
  output$InfoErrorCounts <- renderInfoBox({
    
    tmp = dataInput()
    data = tmp$data
    check = tmp$check
    cond = (!is.null(data$counts) && nrow(data$counts)>0)
    res =infoBox(h6(strong("Counts table")), subtitle = h6("Load the counts table") ,color = "light-blue",width=NULL,fill=TRUE, icon = icon("upload"))
    
    if(cond)
    {
      if(!is.null(check$CheckCounts$Warning)) res = infoBox(h6(strong("Counts table")), subtitle = h6(check$CheckCounts$Warning), icon = icon("warning"),color = "orange",width=NULL,fill=TRUE)
      if(!is.null(check$CheckCounts$Error)) res = infoBox(h6(strong("Counts table")), subtitle = h6(check$CheckCounts$Error), icon = icon("thumbs-o-down"),color = "red",width=NULL,fill=TRUE)
      if(is.null(check$CheckCounts$Error) && is.null(check$CheckCounts$Warning)) res = infoBox(h6(strong("Counts table")), subtitle = h6(paste("Format of the counts table seems to be OK")), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
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
        menuItem("Krona plot", icon = icon("pie-chart"), tabName = "Krona")

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
    counts = dataInput()$data$counts
    if (is.null(inFile)) return(NULL)
    
    
    data = read.csv(inFile$datapath,sep="\t",header=TRUE)
    rownames(data) <- as.character(data[, 1])
    #ord = order(rownames(data))
    #data = data[ord,]
    ### A SUPPRIMER 
    #rownames(data) <- colnames(counts)
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
    else infoBox(h6(strong("Target file")), subtitle = h6("Label of the target file must correspond to counts table column names") ,color = "light-blue",width=NULL,fill=TRUE, icon = icon("warning"))
  })




  ## target table
  output$DataTarget <- renderDataTable(
  dataInputTarget(),
  options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                 pageLength = 10,scrollX=TRUE
  ))

  ## Counts table for the selected taxonomy level
  output$CountsMerge <- renderDataTable(
    dataMergeCounts()$counts,
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE
    ))


  ## Box for target visualisation
  output$BoxCountsMerge <- renderUI({
    
    counts = dataMergeCounts()$counts
    taxo = input$TaxoSelect
    
    if(!is.null(counts) && taxo != "...")
    {
      box(title=paste("Counts table (",taxo,")",sep=""),width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = TRUE,
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


#################################################
##        FOR PIERRE
#################################################
#   
#   ## Merge counts data
#   dataMergeCounts_pierre <-reactive({ 
#     resDiff = ResDiffAnal()
#     dds = resDiff$dds
#     counts = round(counts(dds, normalized = TRUE))
#     
#     VarInt = input$VisuVarIntBoxP
#     ind_taxo = rownames(counts)
#     
#     tmp_merge = GetDataToPlot(resDiff,VarInt,ind_taxo,aggregate=TRUE)
#     counts_tmp_combined = tmp_merge$counts
#     
#     return(counts_tmp_combined)
#   })
# 
# 
# 
# 
# output$CountsMerge_pierre <- renderDataTable(
#   dataMergeCounts_pierre(),
#   options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
#                  pageLength = 10,scrollX=TRUE
#   ))
# 
# 
# ## Box for target visualisation
# output$BoxCountsMerge_pierre <- renderUI({
#     box(title="Counts table",width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = TRUE,
#         dataTableOutput("CountsMerge_pierre"),
#         downloadButton('ExportPloted', 'Export')
#     )
# })
# 
# 
# ## Export in .csv
# output$ExportPloted <- downloadHandler(
#   filename = function() { 'CountsMerge.csv' },
#   content = function(file){write.csv(dataMergeCounts_pierre(), file, sep='\t')}
# )

#################################################
##       END FOR PIERRE
#################################################


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
    tmp = read.table(namesfile,header=TRUE)
    Contrast = colnames(as.matrix(tmp))
    updateSelectInput(session, "ContrastList","Contrasts",Contrast)
    updateSelectInput(session, "ContrastList_table","Contrasts",Contrast)
    updateSelectInput(session, "ContrastList_table_Visu","For which contrasts",Contrast)
    updateSelectInput(session, "ContrastList_table_FC","Contrasts (Min = 2)",Contrast)
  })

  ## Add contrast 
  observeEvent(input$AddContrast,{  
    
    AddCont()
    
  })


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
  })


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
    
  })


  # Infobox Contrast
  output$InfoContrast <- renderInfoBox({
    
    test = FALSE
    res = infoBox("Contrasts", subtitle = h6("At least one contrast (non null) must be defined"), icon = icon("warning"),color = "light-blue",width=NULL,fill=TRUE)
    input$AddContrast
    input$RemoveContrast
    input$fileContrast
    
    filesize = isolate(file.info(namesfile)[,"size"])
    if(is.na(filesize)){filesize=0}
    if(filesize!=0) 
    {
      tmp = read.table(namesfile,header=TRUE)
      if(any(as.vector(tmp)!=0)) test = TRUE
    }
    
    if(test) 
    {
      res = infoBox("Contrasts", subtitle = h6("Contrasts OK"), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
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
      box(title="Contrasts (advanced user)",width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = TRUE,
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
    target = dataInputTarget()
    design = GetDesign(input)
   
    Get_dds_object(input,counts,target,design,normFactors,CT_noNorm)

    
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
    
    int = input$InterestVar
    if(length(int)>=2) intSel = int[c(1,2)]
    else intSel = int[1]
    
    selectizeInput("VarInt",h6(strong("Select the variables of interest (max 2)")),int, selected = intSel,multiple = TRUE,options = list(maxItems = 2))
    
  })
  
  
  output$PlotDiag <- renderPlot({
    input$RunDESeq
    
    resDiff = isolate(ResDiffAnal())
    Plot_diag(input,resDiff)
  })

  output$PlotpcoaEigen <- renderPlot({
    
    resDiff = ResDiffAnal()
    Plot_diag_pcoaEigen(input,resDiff)
  })

  output$PlotEigen <- renderPlot({
    
    resDiff = ResDiffAnal()
    Plot_diag_Eigen(input,resDiff)
  })

  SizeFactor_table <-reactive({ 
    res = ResDiffAnal()
    return(t(data.frame(Factor=res$normFactors)))
    
  })

  output$SizeFactTable <- renderDataTable(
    SizeFactor_table(),
    options = list(scrollX=TRUE,searching = FALSE
  ))


  ## Select Modality PCOA

  output$ModMat <- renderUI({
    
    VarInt = input$VarInt
    target = dataInputTarget()
    
    Mod = list()
    
    for(i in 1:length(VarInt)){
      value = as.character(unique(as.factor(target[,VarInt[i]])))
      Mod[[i]] = selectizeInput(paste("Mod",VarInt[i],sep=""),VarInt[i],value,selected=value, multiple = TRUE)
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
    filename <- function() { paste(input$DiagPlot,paste('meta16S',input$Exp_format,sep="."),sep="_") },
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
  filename <- function() { paste("test",'meta16S.ps',sep="_")},
  content <- function(file) {
    resDiff = ResDiffAnal()
    BaseContrast = read.table(namesfile,header=TRUE)
    #ggsave(filename = filename, Plot_Visu_Heatmap_FC(input,BaseContrast,resDiff),width = input$widthHeat, height = input$heightHeat)
    postscript(file, width = input$widthHeat, height = input$heightHeat)
    if(input$HeatMapType=="Counts")  Plot_Visu_Heatmap(input,resDiff)
    if(input$HeatMapType=="Log2FC")     Plot_Visu_Heatmap_FC(input,BaseContrast,resDiff)
    dev.off()
  }
)
  
  
  ## PNG
  output$exportPNGVisu <- downloadHandler(
    filename <- function() { paste("test",'meta16S.png',sep="_") },
    content <- function(file) {
      png(file, width = 600, height = 400)
      Plot_Visu(input,ResDiffAnal())
      dev.off()
    }
  )

#### Export Visu
output$exportVisu <- downloadHandler(
  filename <- function() { paste(input$PlotVisuSelect,paste('meta16S',input$Exp_format_Visu,sep="."),sep="_") },
  content <- function(file) {
    BaseContrast = read.table(namesfile,header=TRUE)
    taxo = input$TaxoSelect
    
    if(input$Exp_format_Visu=="png") png(file, width = input$widthVisuExport, height = input$heightVisuExport)
    if(input$Exp_format_Visu=="pdf") pdf(file, width = input$widthVisuExport/96, height = input$heightVisuExport/96)
    if(input$Exp_format_Visu=="eps") postscript(file, width = input$widthVisuExport/96, height = input$heightVisuExport/96)
    if(input$Exp_format_Visu=="svg") svg(file, width = input$widthVisuExport/96, height = input$heightVisuExport/96)

    if(input$PlotVisuSelect=="Barplot") print(Plot_Visu_Barplot(input,ResDiffAnal())$gg)
    if(input$PlotVisuSelect=="Heatmap"){
      if(input$HeatMapType=="Counts") print(Plot_Visu_Heatmap(input,ResDiffAnal(),export=TRUE))
      if(input$HeatMapType=="Log2FC") print(Plot_Visu_Heatmap_FC(input,BaseContrast,ResDiffAnal(),export=TRUE))
    } 
    
    if(input$PlotVisuSelect=="Boxplot") print(Plot_Visu_Boxplot(input,ResDiffAnal()))
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
  filename = function() {paste(input$WhichExportTable,input$ContrastList_table,'meta16S.csv',sep="_")},
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
  target = dataInputTarget()
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
    if(!is.null(resDiff$dds)) Plot_Visu_Barplot(input,resDiff)$plotd3
  },env=new.env())

# output$PlotVisu <- renderPlotly({
#   resDiff = ResDiffAnal()
#   if(!is.null(resDiff$dds)) Plot_Visu_Barplot(input,resDiff)
# })
  
  output$heatmap <- renderD3heatmap({
    resDiff = ResDiffAnal()
    BaseContrast = read.table(namesfile,header=TRUE)
    resplot = NULL
    if(!is.null(resDiff$dds))
    { 
      if(input$HeatMapType=="Counts") resplot = Plot_Visu_Heatmap(input,resDiff)
      if(input$HeatMapType=="Log2FC") resplot = Plot_Visu_Heatmap_FC(input,BaseContrast,resDiff)
    }
    return(resplot)
  },env=new.env())
  


  output$Boxplot <- renderPlot({
    resDiff = ResDiffAnal()
    if(!is.null(resDiff$dds)) Plot_Visu_Boxplot(input,resDiff)
  },height=reactive(input$heightVisu))


  output$DiversityPlot <- renderPlot({
    resDiff = ResDiffAnal()
    if(!is.null(resDiff$dds)) Plot_Visu_Diversity(input,resDiff,type="point")
  })
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  output$RarefactionPlot <- renderPlot({
    resDiff = ResDiffAnal()
    taxo = input$TaxoSelect
    if(!is.null(resDiff)) Plot_Visu_Rarefaction(input,resDiff,ranges$x,ranges$y,ylab=taxo)
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



  output$SelectVarBoxDiv <- renderUI({
    
    selectVar = input$VisuVarInt
    
    if(!is.null(selectVar)) 
    {
      selectInput("VarBoxDiv", h6(strong("By")),selectVar)
    }
    
  })
  
  output$plotVisu <- renderUI({
    
    res=NULL
    if(input$PlotVisuSelect=="Barplot") res =  showOutput("PlotVisuBar")
    if(input$PlotVisuSelect=="Heatmap") res =  d3heatmapOutput("heatmap", height = input$heightVisu+10)
    if(input$PlotVisuSelect=="Boxplot") res =  plotOutput("Boxplot", height = input$heightVisu+10)
    if(input$PlotVisuSelect=="Diversity") res =  plotOutput("DiversityPlot", height = input$heightVisu+10)
    if(input$PlotVisuSelect=="Rarefaction") res = plotOutput("RarefactionPlot",dblclick = "RarefactionPlot_dblclick",brush = brushOpts(id = "RarefactionPlot_brush",resetOnNew = TRUE), height = input$heightVisu+10)
    
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
    BaseContrast = read.table(namesfile,header=TRUE)
    
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

        SelContrast = input$ContrastList_table_Visu
        padj = Get_log2FC_padj(input,BaseContrast,resDiff, info = NULL)$padj
        cont = which(colnames(padj)%in%SelContrast)
        padj = padj[,cont] 
        selTaxo = names(padj[which(padj<=input$AlphaVal)])
        res = selectizeInput("selectTaxoPlot",h6(strong(paste("Select the",input$TaxoSelect, "to plot"))),Available_taxo, selected = selTaxo,multiple = TRUE,options = list(minItems = 2))
      }      
      if(input$SelectSpecifTaxo=="All") res = selectizeInput("selectTaxoPlot",h6(strong(paste("Select the",input$TaxoSelect, "to plot"))),Available_taxo, selected = Available_taxo,multiple = TRUE)
    }
    return(res)
})


  
  output$VarIntVisu <- renderUI({
    
    int = input$InterestVar
    if(length(int)>=2) intSel = int[c(1,2)]
    else intSel = int[1]
    
    selectizeInput("VisuVarInt",h6(strong("Select the variables of interest")),int, selected = intSel,multiple = TRUE)
    
  })


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