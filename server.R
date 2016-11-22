source('LoadPackages.R')

shinyServer(function(input, output,session) {

  hide(id = "loading-content", anim = TRUE, animType = "fade",time=1.5)
  hide(id = "loading-content-bar", anim = TRUE, animType = "fade",time=1.5)
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
    
    data = NULL
    inFile <- input$fileCounts
    
    if (is.null(inFile)) return(NULL)
    
    
    try(read.csv(inFile$datapath,sep=input$sepcount,header=TRUE,check.names=FALSE)->data,silent=T)

    if(!is.null(data)){
      colnames(data) = gsub("-",".",colnames(data))
      ## Rownames
      if(!TRUE%in%duplicated(data[,1])) rownames(data)=data[,1];data=data[,-1]
    }
    return(as.data.frame(data))
  })
  
  
  
  ## Taxo File
  dataInputTaxo <-reactive({ 
    
    inFile <- input$fileTaxo
    
    if (is.null(inFile)) return(NULL)
    
    if(input$TypeTaxo=="Table") 
    {
      try(read.csv(inFile$datapath,sep=input$septaxo,header=TRUE)->data,silent=T)
    
      ## Rownames
      if(!is.null(data))
      {
        if(!TRUE%in%duplicated(data[,1])){ 
          DataNames=data[,1]
          colNames=colnames(data)[-1]
          data=as.matrix(data[,-1])
          rownames(data)=DataNames
          colnames(data) = colNames
        }
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
    
    data = NULL
    inFile <- input$fileBiom
    
    if (is.null(inFile)) return(NULL)
    try(read_biom(inFile$datapath)->data,silent=T)
    
    return(data)
  })
  
  
  
  ## Input data
  dataInput <-reactive({ 
    
    data = NULL
    check = NULL
    percent = NULL
    Taxo = NULL
    Counts = NULL
    if(input$FileFormat=="fileCounts")
    {
      Counts = dataInputCounts()
      if(!input$NoTaxoFile) Taxo = dataInputTaxo()
      if(!is.null(Counts) && input$NoTaxoFile) {Taxo = data.frame(rownames(Counts),row.names = rownames(Counts));names(Taxo)=NA}
      
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
  
  
  
  
  
  ## Size factor file (optional)
  dataSizeFactors <-reactive({ 
    
    inFile <- input$fileSizeFactors
    
    if (is.null(inFile)) return(NULL)
    
    data = read.csv(inFile$datapath,sep=input$sepsize,header=TRUE)
    return(as.data.frame(data))
  })
  
  
  ## Size factor file (optional)
  SizeFactors_fromFile <-reactive({ 
    
    Error = NULL
    Check = TRUE
    
    data = dataSizeFactors()
    normFactors = dataMergeCounts()$normFactors
    
    if(!is.null(data)){
      ## Check the format
      
      tmp = as.numeric(data)
      names(tmp) = colnames(data)
      
      if(length(tmp)!=length(normFactors)){Error = "The number of samples is not the same than in the target file, size factors will be estimated"; Check = FALSE}
      if(!identical(names(tmp),names(normFactors))){Error = "The names are not the same or in the same order than in the target file, size factors will be estimated"; Check = FALSE}
      
      if(Check) normFactors = tmp
    }
    
    return(list(Check = Check,Error = Error,normFactors=normFactors))
  })
  
  
  
  
  
  
  ## Merge counts data
  dataMergeCounts <-reactive({
    input$RunDESeq
    
    counts = NULL
    CheckTarget = FALSE
    normFactors = NULL
    CT_noNorm = NULL
    CT_Norm = NULL
    ChTM = NULL
    data = isolate(dataInput()$data)
    target = isolate(dataInputTarget()$target)
    labeled= isolate(dataInputTarget()$labeled)
    taxo = isolate(input$TaxoSelect)
    withProgress(
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0 && !is.null(taxo) && taxo!="..." && !is.null(target)) 
    {
      design = GetDesign(isolate(input))
      ChTM = CheckTargetModel(input,target,labeled,data$counts)$Error
      if(!is.null(design) && is.null(ChTM))
      {
        tmp = isolate(GetCountsMerge(input,data,taxo,target,design))
        counts = tmp$counts
        ## Filtering the counts
        if(isolate(input$AddFilter) && !is.null(isolate(input$SliderThSamp)) && !is.null(isolate(input$SliderThAb)))
        {
          ind.filter =Filtered_feature(counts,isolate(input$SliderThSamp),isolate(input$SliderThAb))$ind
          counts = counts[-ind.filter,]
        }
        CheckTarget = tmp$CheckTarget
        #target = tmp$target
        #labeled = tmp$labeled
        normFactors = tmp$normFactors
        ## OTU table, norm and no norm
        CT_noNorm = tmp$CT_noNorm
        CT_Norm = tmp$CT_Norm
      }
    }
    ,message="Merging the counts ...")
    return(list(counts=counts,CheckTarget=CheckTarget,normFactors=normFactors,CT_noNorm=CT_noNorm, CT_Norm=CT_Norm))
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
    res = infoBox(h6(strong("Taxonomy table")), subtitle = h6("Load the taxonomy table") ,color = "light-blue",width=NULL,fill=TRUE, icon = icon("upload"))
    
    if(cond)
    {
      if(!is.null(check$CheckTaxo$Warning)) res = infoBox(h6(strong("Taxonomy table")), subtitle = h6(check$CheckTaxo$Warning), icon = icon("warning"),color = "orange",width=NULL,fill=TRUE)
      if(!is.null(check$CheckTaxo$Error)) res = infoBox(h6(strong("Taxonomy table")), subtitle = h6(check$CheckTaxo$Error), icon = icon("thumbs-o-down"),color = "red",width=NULL,fill=TRUE)
      if(is.null(check$CheckTaxo$Error) && is.null(check$CheckTaxo$Warning)) res = infoBox(h6(strong("Taxonomy table")), subtitle = h6(paste("Format of the taxonomy table seems to be OK")), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
    }
    
    if(input$NoTaxoFile && input$FileFormat=="fileCounts") res = infoBox(h6(strong("Taxonomy table")), subtitle = h6("No taxonomy table has been uploaded, the analysis can only be done at the OTU/gene level"), icon = icon("warning"),color = "orange",width=NULL,fill=TRUE)
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
  
  
  ####### Filtering the counts (sliders)
  
  output$ThAb <- renderUI({
    input$AddFilter
    
    res = NULL
    counts = isolate(dataMergeCounts()$counts)
    tot = rowSums(counts)
    #save(counts,tot,file="testFilter.RData")
    withProgress({tmp = SelectThreshAb(counts,lambda=max(round(sum(counts)/nrow(counts)*0.05),min(tot)+1),graph=FALSE)},message="Loading...")
    
    res = sliderInput("SliderThAb","Threshold on the total abundance (in log)",min=0,max=round(max(log(tot+1)),1),value = log(tmp+1))
    return(res)
  })
  
  
  output$ThSamp <- renderUI({
    input$AddFilter
    
    res = NULL
    counts = isolate(dataMergeCounts()$counts)
    counts.bin = as.matrix(counts)
    counts.bin[which(counts>0)] = 1
    nbSampByFeat = rowSums(counts.bin)
    
    ## Default value
    val = round(max(nbSampByFeat)*0.2)
    
    res = sliderInput("SliderThSamp","Threshold on the minimal number of samples",min=0,max=max(nbSampByFeat),value = val)
    return(res)
  })
  
  
  ## Plot for the filtering step$
  
  # plot_filter(counts,th.samp,th.abund,type="Scatter")
  
  output$Plot_ThAb <- renderPlot({
    counts = dataMergeCounts()$counts
    ## output of plot_filter is ggplot class
    plot_filter(counts,input$SliderThSamp,input$SliderThAb,type="Abundance")
    
  })
    
  output$Plot_ThSamp <- renderPlot({
    counts = dataMergeCounts()$counts
    ## output of plot_filter is ggplot class
    plot_filter(counts,input$SliderThSamp,input$SliderThAb,type="Samples")
  })
  
  output$Plot_Scatter_Filter <- renderScatterD3({
    counts = dataMergeCounts()$counts
    ## output of plot_filter is ggplot class
    plot_filter(counts,input$SliderThSamp,input$SliderThAb,type="Scatter")
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

      sidebarMenu(id = "side",
        menuItem("Statistical analysis",
                 menuSubItem("Run differential analysis",tabName="RunDiff"),
                 menuSubItem("Diagnostic plots",tabName="DiagPlotTab"),
                 menuSubItem("Tables",tabName="TableDiff"),
                 icon = icon("bar-chart-o"), tabName = "AnaStat"
        ),
        menuItem("Visualization",icon = icon("area-chart"),
                 menuSubItem("Global views",tabName="GlobVisu"),
                 menuSubItem("Comparison plots",tabName="CompPlot"),
                 tabName = "Visu"),
        menuItem("Perspective plots", icon = icon("pie-chart"), tabName = "Krona")
      )
    } else{ sidebarMenu(id = "side",NULL)}

  })
  
  
  
  #####################################################
  ##
  ##                DATA TABLE
  ##
  #####################################################
  
  ## Counts Table
  output$DataCounts <- DT::renderDataTable(
    dataInput()$data$counts, 
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE, processing=FALSE
    ))
  
  ## Counts Table
  output$DataVenn<- DT::renderDataTable(#{
    #SelContrast = input$ContrastList_table_FC
    #resDiff = ResDiffAnal()
    #BaseContrast = read.table(namesfile,header=TRUE)
    GetData_venn(input,input$ContrastList_table_FC,read.table(namesfile,header=TRUE),ResDiffAnal())$df.tot,
  #}
   options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                    pageLength = 10,scrollX=TRUE, processing=FALSE
  ))
  
  
  ## Taxonomy table
  output$DataTaxo <- DT::renderDataTable(
    dataInput()$data$taxo, 
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE, processing=FALSE
    ))
  
  
  ## Tab box for data visualisation
  output$TabBoxData <- renderUI({
    
    data=dataInput()$data
    
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0)
    {
      tabBox(width = NULL, selected = "Count table",
             tabPanel("Count table",DT::dataTableOutput("DataCounts")),
             tabPanel("Taxonomy",DT::dataTableOutput("DataTaxo")),
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
      tmpPercent = round(apply(is.na(taxo),2,table)["FALSE",]/(nrow(taxo)-1)*100,2)
      
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
    labeled = 0
    
    if (is.null(inFile)) return(NULL)
    
    ## Read the data
    try(read.csv(inFile$datapath,sep=input$septarget,header=TRUE)->data,silent=TRUE)
    
    if(!is.null(data))
    {
      data = as.data.frame(data)
      names = colnames(data)
      
      ## Change the rownames
      rownames(data) <- as.character(data[, 1])
      
      ## Keep only the row which are in the count table
      ind = which(rownames(data)%in%colnames(counts))
      data = as.data.frame(data[ind,])
      colnames(data) = names
      ## Replace "-" by "."
      if(ncol(data)>1 && nrow(data)>1){
        ind_num = which(sapply(as.data.frame(data[,-1]),is.numeric)) + 1
        if(length(ind_num)>0){
          data_tmp =cbind( as.data.frame(apply(as.data.frame(data[,-ind_num]),2,gsub,pattern = "-",replacement = ".")),data[,ind_num])
          colnames(data_tmp) = c(colnames(data)[-ind_num],colnames(data)[ind_num])
          data = data_tmp
        }
        if(length(ind_num)==0){data = as.data.frame(apply(data,2,gsub,pattern = "-",replacement = "."))}
      }
      target = as.data.frame(data)
      
      # target = as.data.frame(apply(target,2,gsub,pattern = "-",replacement = "."))
      
      #ord = order(rownames(data))
      #data = data[ord,]
      ### A SUPPRIMER 
      #rownames(data) <- colnames(counts)
      
      # Percent annotated
      #     print(ind)
      #     print(colnames(counts))
      #     print(rownames(data))
      labeled = length(ind)/length(colnames(counts))*100.0
    }
    
    return(list(target = target, labeled=labeled))
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
    VarInt = input$InterestVar
    res = NULL
    if(!is.null(target) && length(input$InterestVar)>1) 
    {
      Interac = GetInteraction2(target,VarInt)
      res = selectInput("Interaction2",h6(strong("Add interactions")),Interac,selected=NULL,multiple=TRUE)
    }
    if(length(input$InterestVar)==1) res = NULL
    
    return(res)
  })


  ## Var for normalization
  output$SelectVarNorm <- renderUI({
    
    target=dataInputTarget()$target
    res = selectInput("VarNorm",h6(strong("Normalization by:")),NULL,multiple=TRUE)
    if(!is.null(target)) 
    {
      namesTarget = colnames(target)[2:ncol(target)]
      ind = which(apply(as.data.frame(target[,namesTarget]),2,is.numeric))
      if(length(ind)>=1) namesTarget = namesTarget[-ind]
      res = selectInput("VarNorm",h6(strong("Normalization by:")),c(NULL,namesTarget),multiple=TRUE)
    }
    return(res)
    
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
  output$DataTarget <- DT::renderDataTable(
    dataInputTarget()$target,
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE, processing=FALSE
    ))
  
  ## Counts table for the selected taxonomy level
  output$CountsMerge <- DT::renderDataTable(
    round(counts(ResDiffAnal()$dds,normalized=TRUE)),
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE, processing=FALSE
    ))
  
  
  ## Box for merged counts
  output$BoxCountsMerge <- renderUI({
    input$RunDESeq
    counts = isolate(dataMergeCounts()$counts)
    taxo = input$TaxoSelect
    
    if(!is.null(counts) && taxo != "...")
    {
      box(title=paste("Count table (",taxo,")",sep=""),width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = TRUE,
          DT::dataTableOutput("CountsMerge"),
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
  
  ## Export size factors
  output$ExportSizeFactor <- downloadHandler(
    filename = function() { if (input$sepsizef == "\t") 'SHAMAN_sizefactors.tsv' else 'SHAMAN_sizefactors.csv' },
    content = function(file){write.table(SizeFactor_table(), file,quote=FALSE,row.names = FALSE,sep=input$sepsizef)}
  )
  
  
  ## Box for target visualisation
  output$BoxTarget <- renderUI({
    
    target = dataInputTarget()$target
    
    if(!is.null(target) &&  nrow(target)>0)
    {
      box(title="Target file overview",width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = TRUE,
          DT::dataTableOutput("DataTarget")
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
    updateSelectInput(session, "ContrastList_table_VisuComp","For which contrasts",Contrast)
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
    target = resDiff$target
    names = resultsNames(dds)
    
    BaseContrastEasy(input,names,namesfile,target)
    filesize = file.info(namesfile)[,"size"]
    if(is.na(filesize)){filesize=0}
    if(filesize!=0) tmp = read.table(namesfile,header=TRUE)
    Contrast = colnames(as.matrix(tmp))
    
    updateSelectInput(session, "ContrastList","Contrasts",Contrast)
    updateSelectInput(session, "ContrastList_table","Contrasts",Contrast)
    updateSelectInput(session, "ContrastList_table_Visu","For which contrasts",Contrast)
    updateSelectInput(session, "ContrastList_table_VisuComp","For which contrasts",Contrast)
    updateSelectInput(session, "ContrastList_table_FC","Contrasts (Min = 2)",Contrast)
  })
  
  ## Add contrast 
  observeEvent(input$AddContrastEasy,{  
    
    AddContEasy()
    
  },priority=1)
  
  
  AddContFromFile <-eventReactive(input$fileContrast,{ 
    
    res = ReadContrastFile()
    resDiff = ResDiffAnal()
    dds = resDiff$dds
    CheckCont = CheckContrast(res,dds)
    
    if(is.null(CheckCont$Error))
    {
      res = CheckCont$contrastFile
      createdCont = NULL
      if(!is.null(res))
      { 
        filesize = file.info(namesfile)[,"size"]
        if(is.na(filesize)){filesize=0}
        if(filesize!=0){ createdCont = read.table(namesfile,header=TRUE) }
      
        if(!is.null(createdCont)) res = cbind(res,createdCont)
        updateSelectInput(session, "ContrastList","Contrasts",colnames(res))
        updateSelectInput(session, "ContrastList_table","Contrasts",colnames(res))
        updateSelectInput(session, "ContrastList_table_Visu","For which contrasts",colnames(res))
        updateSelectInput(session, "ContrastList_table_VisuComp","For which contrasts",colnames(res))
        updateSelectInput(session, "ContrastList_table_FC","Contrasts (Min = 2)",colnames(res))
        write.table(res,namesfile,row.names=FALSE)
      }
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
      updateSelectInput(session, "ContrastList_table_VisuComp","For which contrasts",ContrastKept)
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
    updateSelectInput(session, "ContrastList_table_VisuComp","For which contrasts",NULL)
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
    input$AddContrastEasy
    input$RemoveContrast
    input$fileContrast
    resDiff = ResDiffAnal()
    dds = resDiff$dds
    res=NULL
    tmpFile = NULL
    
    if(!is.null(resDiff)){
      
      res = infoBox("Contrasts", subtitle = h6("At least one contrast (non null) must be defined"), icon = icon("warning"),color = "light-blue",width=NULL,fill=TRUE)
      
      filesize = isolate(file.info(namesfile)[,"size"])

      if(is.na(filesize)){filesize=0}
      if(filesize!=0) tmpFile = read.table(namesfile,header=TRUE)
      if(!is.null(tmpFile))
      {
        CheckCont = CheckContrast(tmpFile,dds)
        if(!is.null(CheckCont$Warning)) res = infoBox(h6(strong("Contrasts")), subtitle = h6(CheckCont$Warning), icon = icon("warning"),color = "orange",width=NULL,fill=TRUE)
        if(!is.null(CheckCont$Error)) res = infoBox(h6(strong("Contrasts")), subtitle = h6(CheckCont$Error), icon = icon("thumbs-o-down"),color = "red",width=NULL,fill=TRUE)
        if(is.null(CheckCont$Error) && is.null(CheckCont$Warning))  res = infoBox("Contrasts", subtitle = h6("Contrasts OK"), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
      
      }
      ## if user load a bad contrast file after having define one or more good contrasts
#       if(!is.null(input$fileContrast)){
#         tmpRead = ReadContrastFile()
#         CheckCont_new = CheckContrast(tmpRead,dds)
#         if(!is.null(CheckCont_new$Warning)) info("test1")
#         if(!is.null(CheckCont_new$Error)) info("test2")
#       }
    }
    return(res)
  })
  
  
  # Infobox Contrast
  output$InfoContrast_box <- renderUI({
    resDiff = ResDiffAnal()
    dds = resDiff$dds
    
    if(!is.null(resDiff)){

      if(!is.null(input$fileContrast)){
        tmpRead = ReadContrastFile()
        CheckCont_new = CheckContrast(tmpRead,dds)
        if(!is.null(CheckCont_new$Warning)){      
          box(title = "Warning", status = "warning",width = NULL,
              h6(strong(CheckCont_new$Warning)))
        }
        if(!is.null(CheckCont_new$Error)){      
          box(title = "Warning", status = "warning",width = NULL,
              h6(strong(CheckCont_new$Error)))
        }
      }
    }
  })
  
  
  output$contrastBox <- renderUI({
    
    resDiff = ResDiffAnal()
    int = input$Interaction2
    
    if(!is.null(resDiff))
    { 
      ## Check the R version
       if(as.numeric(R.Version()$major)<=3 && as.numeric(R.Version()$minor) <=1.2){
        box(title="Contrasts (New)",width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = FALSE,
            fluidRow(
              column(width=3,selectInput("Select1_contrast","Compare","")),
              column(width=3,selectInput("Select2_contrast","To","")),
              if(length(int)>=1) column(width=3,selectInput("Select3_contrast",label=h6(strong("For")),"")),
              column(width=3,br(),actionButton("AddContrastEasy","Add",icon = icon("plus")))
            )
        )
      }
    }
    
    
    
  })
  
  ModifMod_ContEasy <-eventReactive(input$Select1_contrast,{
    input$RunDESeq
    resDiff = ResDiffAnal()
    int = input$Interaction2
    target = as.data.frame(resDiff$target)
    
    InterVar = input$InterestVar
    
    
    
    ## Get the selected variable from the selected modality
    Sel_Var = InterVar[which(unlist(lapply(as.data.frame(target[,InterVar]),FUN = function(x){input$Select1_contrast%in%x})))]
    
    ModInterestCond = levels(sapply(target[,Sel_Var],as.factor))
    ModInterestCond = ModInterestCond[-which(ModInterestCond==input$Select1_contrast)]
    
    updateSelectInput(session,"Select2_contrast","To",ModInterestCond)
  })
  
  
  observeEvent(input$Select1_contrast,{ 
    
    ModifMod_ContEasy()
  })
  
  
  ModifMod_ContEasyFrom <-eventReactive(input$RunDESeq,{
    
    resDiff = ResDiffAnal()
    int = input$Interaction2
    target = as.data.frame(resDiff$target)
    
    InterVar = input$InterestVar
    
    ## Remove numeric variable
    ind = unlist(lapply(as.data.frame(target[,InterVar]),is.numeric))
    InterVar = InterVar[!ind]
    target_int = lapply(as.data.frame(target[,InterVar]),as.factor)
    ModInterestAll = unique(unlist(lapply(target_int,levels)))
    
    updateSelectInput(session, "Select1_contrast",label="Compare",ModInterestAll)
  })
  
  
  observeEvent(input$RunDESeq,{ 
    
    ModifMod_ContEasyFrom()
  })
  
  
  
  ModifMod_ContEasyFor <-eventReactive(input$Select1_contrast,
    {
    
    resDiff = ResDiffAnal()
    int = input$Interaction2
    ModInterestFor = "All"
    target = as.data.frame(resDiff$target)
    
    InterVar = input$InterestVar
    
    ## Get the selected variable from the selected modality
    Sel_Var = InterVar[which(unlist(lapply(as.data.frame(target[,InterVar]),FUN = function(x){input$Select1_contrast%in%x})))]
   
    
    ## Keep only the variables in interactoin with Sel_Var
    if(!is.null(Sel_Var) && length(int)>0 && length(Sel_Var)>0){
      indInter = grep(Sel_Var,int)
      if(length(indInter)>0) int = int[indInter]
      var_Inter = unique(unlist(strsplit(int,":"))) 
      var_Inter = var_Inter[-which(var_Inter%in%Sel_Var)]
      
      ## remove if numeric
      if(length(var_Inter)>1){ind = unlist(lapply(as.data.frame(target[,var_Inter]),is.numeric));var_Inter = var_Inter[!ind]}
      if(length(var_Inter)==1){ind = is.numeric(target[,var_Inter]);var_Inter = var_Inter[!ind]}
      

      if(length(var_Inter)>=1)  ModInterestFor = c("All",unique(unlist(lapply(as.data.frame(target[,var_Inter]),levels))))

    }
    
    updateSelectInput(session,"Select3_contrast","For",ModInterestFor)
  })
  
  
  observeEvent(input$Select1_contrast,{ 
    
    ModifMod_ContEasyFor()
  })
  
  
  output$contrastBoxAdvanced <- renderUI({
    
    resDiff = ResDiffAnal()
    
    if(!is.null(resDiff))
    { 
      box(title="Contrasts (advanced user)",width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = FALSE,
          fluidRow(
            column(width=9,
                   fileInput('fileContrast', h6(strong('Select a file of contrasts')),width="80%")
            ),
            column(width=3,
                   column(width=12,selectInput("sepContFile", h6(strong("Separator:")), c("Tab" = "\t", "Comma" = ",", "Semicolon" = ";","Space"= " "),selected = " "))
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
    
    try(read.csv(inFile$datapath,header=TRUE,sep=input$sepContFile)->res,silent=T)
    
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
  ResDiffAnal <-eventReactive(input$RunDESeq,withProgress({
    
    target = dataInputTarget()$target
    design = GetDesign(input)
    dMC = dataMergeCounts()
    counts = dMC$counts
    CT_noNorm = dMC$CT_noNorm
    CT_Norm = dMC$CT_Norm
    
    ## If no file, size factors are estimated
    normFactors = SizeFactors_fromFile()$normFactors
    
    Get_dds_object(input,counts,target,design,normFactors,CT_noNorm,CT_Norm)
    
    
  },message = "Analysis in progress"))
  
  
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
      ind = which(is.na(tmp))
      if(length(ind)>0) tmp = tmp[-which(is.na(tmp))]
      selectInput("TaxoSelect",h6(strong("Select the taxonomy level")),c("...",tmp,"OTU/Gene"))
    }
    else selectInput("TaxoSelect",h6(strong("Select the taxonomy level")),c("..."))
    
  })
  
  
  
  
  # Infobox taxo
  output$InfoTaxo <- renderInfoBox({
    input$RunDESeq
    data = isolate(dataInput()$data)
    taxo = isolate(input$TaxoSelect)
    
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0 && !is.null(isolate(taxo)) && taxo!="...") 
    {
      counts = isolate(dataMergeCounts()$counts)
      nfeature = nrow(counts)
      infoBox(h6(strong("Taxonomy")), subtitle = h6(paste(taxo, ", nb features: ",nfeature,sep="")), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
    }
    else infoBox(h6(strong("Taxonomy")), subtitle = h6("Select the taxonomy for the analysis") ,color = "light-blue",width=NULL,fill=TRUE, icon = icon("warning"))
  })
  
  
  # Infobox model/target
#   output$InfoModel<- renderInfoBox({
#     res = infoBox(h6(strong("Target file and variables")), 
#                   subtitle = h6(strong("Your target file must contain at least 2 columns and 2 rows. NA's values are not allowed and the variables must not be collinear.")), 
#                   icon = icon("book"),color = "green",width=NULL,fill=TRUE)
#     
#     target = dataInputTarget()$target
#     taxo = input$TaxoSelect
#     ChTM = NULL
#     
#     ## Return NULL if there is no error
#     if(!is.null(target)) ChTM = CheckTargetModel(input,target)
#   
#     if(!is.null(ChTM)) res = infoBox(h6(strong("Error")), subtitle = h6(ChTM), icon = icon("thumbs-o-down"),color = "red",width=NULL,fill=TRUE)
#     
#     return(res)
#     
#     })
  
  output$InfoModel<- renderUI({
    
    CT = dataInput()$data$counts
    target = dataInputTarget()$target
    labeled = dataInputTarget()$labeled
    taxo = input$TaxoSelect
    ChTM = NULL
    
    ## Return NULL if there is no error
    if(!is.null(target)) ChTM = CheckTargetModel(input,target,labeled,CT)
    if(!is.null(ChTM$Error)) {   
      box(title = "Error", status = "danger",width = 6,
          h6(strong(ChTM$Error)),
          footer = em("Reminder: Your target file must contain at least 2 columns and 2 rows. NA's values are not allowed and the variables must not be collinear.")
      )
    } else return(NULL)
    
  })
  
  output$InfoModelHowTo<- renderUI({
    
    CT = dataInput()$data$counts
    target = dataInputTarget()$target
    labeled = dataInputTarget()$labeled
    taxo = input$TaxoSelect
    ChTM = NULL
    
    ## Return NULL if there is no error
    if(!is.null(target)) ChTM = CheckTargetModel(input,target,labeled,CT)
    
    if(!is.null(ChTM$HowTo)) {   
      box(title = "How To", status = "success",width = 6,
          h6(strong(ChTM$HowTo))
      )
    }
    
  })
  
  
  output$InfoBIOM<- renderUI({
    
    if(input$FileFormat=="fileBiom")
    {
      inFile <- input$fileBiom
      tmpBIOM = dataInputBiom()
      
      if(!is.null(inFile) && is.null(tmpBIOM)) {   
        box(title = "Error", status = "danger",width = 3,
            h5(strong("This file can not be loaded.")),br(),
            em("The loaded file is not in the biom format or its format is not currently supported by SHAMAN software")
        )
      }
    }
  })
  
  
  output$InfoCountsFile<- renderUI({
    
    if(input$FileFormat=="fileCounts")
    {
      inFile <- input$fileCounts
      Counts = dataInputCounts()
      
      if(!is.null(inFile) && is.null(Counts)) {   
        box(title = "Error", status = "danger",width = 3,
            h5(strong("This file can not be loaded.")),br(),
            em("The count table file is not in the correct format for SHAMAN software")
        )
      }
    }
  })
  
  
  output$InfoTaxoFile<- renderUI({
    
    if(input$FileFormat=="fileCounts")
    {
      inFile <- input$fileTaxo
      Taxo = dataInputTaxo()
      
      if(!is.null(inFile) && !input$NoTaxoFile && is.null(Taxo)) {   
        box(title = "Error", status = "danger",width = 3,
            h5(strong("This file can not be loaded.")),br(),
            em("The taxonomy table file is not in the correct format for SHAMAN software")
        )
      }
    }
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
  
  
  
  ## Select PCA/PCOA axis
  output$PC1_sel <-renderUI ({
    res = NULL
    resDiff = ResDiffAnal()
    
    if(!is.null(resDiff)){ 
      pca_tab = Plot_diag(input,resDiff,getTable=TRUE)
      if(!is.null(pca_tab)) 
      {
        pc_axes = paste("PC",seq(1,ncol(pca_tab)),sep="")
        res = selectizeInput("PCaxe1","X-axis",pc_axes)
      }
    }
    return(res)
  })
  
  output$PC2_sel <-renderUI ({
    res = NULL
    resDiff = ResDiffAnal()
    
    if(!is.null(resDiff)){ 
      pca_tab = Plot_diag(input,resDiff,getTable=TRUE)
      if(!is.null(pca_tab)) 
      {
        pc_axes = paste("PC",seq(1,ncol(pca_tab)),sep="")
        res = selectizeInput("PCaxe2","Y-axis",pc_axes,selected=pc_axes[min(2,ncol(pca_tab))])
      }
    }
    return(res)
  })
  
  
  
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
  
  output$SizeFactTable <- DT::renderDataTable(
    SizeFactor_table(),
    options = list(scrollX=TRUE,searching = FALSE, processing=FALSE
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
      else if(input$Exp_format=="pdf") pdf(file, width = input$widthDiagExport/96, height = input$heightDiagExport/96)
      else if(input$Exp_format=="eps") postscript(file, width = input$widthDiagExport/96, height = input$heightDiagExport/96)
      else if(input$Exp_format=="svg") svg(file, width = input$widthDiagExport/96, height = input$heightDiagExport/96)
      
      print(Plot_diag(input,ResDiffAnal()))
      dev.off()
    }
  )
  
  #####################################################
  ##
  ##                EXPORT VISU GRAPH
  ##
  #####################################################
  
  
  #### Export Visu
  output$exportVisu <- downloadHandler(
    filename <- function() { paste(input$PlotVisuSelect,paste('SHAMAN',input$Exp_format_Visu,sep="."),sep="_") },
    content <- function(file) {
      
      taxo = input$TaxoSelect
      
      if(input$Exp_format_Visu=="png") png(file, width = input$widthVisuExport, height = input$heightVisuExport)
      else if(input$Exp_format_Visu=="pdf") pdf(file, width = input$widthVisuExport/96, height = input$heightVisuExport/96)
      else if(input$Exp_format_Visu=="eps") postscript(file, width = input$widthVisuExport/96, height = input$heightVisuExport/96,paper="special")
      else if(input$Exp_format_Visu=="svg") svg(file, width = input$widthVisuExport/96, height = input$heightVisuExport/96)
      
      if(input$PlotVisuSelect=="Barplot") print(Plot_Visu_Barplot(input,ResDiffAnal())$gg)
      else if(input$PlotVisuSelect=="Heatmap") Plot_Visu_Heatmap(input,ResDiffAnal(),export=TRUE)
      else if(input$PlotVisuSelect=="Boxplot") print(Plot_Visu_Boxplot(input,ResDiffAnal(),alpha=ifelse(input$Exp_format_Visu=="eps",1,0.7)))
      else if(input$PlotVisuSelect=="Scatterplot") print(Plot_Visu_Scatterplot(input,ResDiffAnal(),export=TRUE,lmEst = FALSE))
      else if(input$PlotVisuSelect=="Diversity") print(Plot_Visu_Diversity(input,ResDiffAnal())$plot)
      else if(input$PlotVisuSelect=="Rarefaction") print( Plot_Visu_Rarefaction(input,ResDiffAnal(),ranges$x,ranges$y,ylab=taxo))
      dev.off()
      
    }
  )
  
  
  
  #### Export Visu
  output$exportVisuComp <- downloadHandler(
    filename <- function() { paste(input$PlotVisuSelectComp,paste('SHAMAN',input$Exp_format_Visu,sep="."),sep="_") },
    content <- function(file) {
      
      taxo = input$TaxoSelect
      
      if(input$Exp_format_VisuComp=="png") png(file, width = input$widthVisuExportComp, height = input$heightVisuExportComp)
      else if(input$Exp_format_VisuComp=="pdf") pdf(file, width = input$widthVisuExportComp/96, height = input$heightVisuExportComp/96)
      else if(input$Exp_format_VisuComp=="eps") postscript(file, width = input$widthVisuExportComp/96, height = input$heightVisuExportComp/96,paper="special")
      else if(input$Exp_format_VisuComp=="svg") svg(file, width = input$widthVisuExportComp/96, height = input$heightVisuExportComp/96)
      
      BaseContrast = read.table(namesfile,header=TRUE)
      filesize = file.info(namesfile)[,"size"]
      if(is.na(filesize)){filesize=0}
      
      if(input$PlotVisuSelectComp=="Venn"){ 
        if(filesize!=0) print(Plot_Visu_Venn(input,BaseContrast,ResDiffAnal(),export=TRUE))
      }
      if(input$PlotVisuSelectComp=="Heatmap_comp"){
        if(filesize!=0) Plot_Visu_Heatmap_FC(input,BaseContrast,ResDiffAnal(),export=TRUE)
      } 
      
      
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
    input$AddContrastEasy
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
  
  ## Significant diff table
  output$DataDiffsignificant <- DT::renderDataTable(
    datatable(dataDiff()$significant,rownames = FALSE),
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE, processing=FALSE
    ))
  ## Complete diff table
  output$DataDiffcomplete <- DT::renderDataTable(
    datatable(dataDiff()$complete,rownames = FALSE),
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE, processing=FALSE
    ))
  
  ## Up diff table
  output$DataDiffup <- DT::renderDataTable(
    datatable(dataDiff()$up,rownames = FALSE),
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE, processing=FALSE
    ))
  
  ## Down diff table
  output$DataDiffdown <- DT::renderDataTable(
    datatable(dataDiff()$down,rownames = FALSE),
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE, processing=FALSE
    ))
  
  
  ## TabBox for diff table
  output$TabBoxDataDiff <- renderUI({
    
    data = dataDiff()
    
    if(!is.null(data))
    {
      
      tabBox(width = NULL, selected = "Significant",
             tabPanel("Significant",DT::dataTableOutput("DataDiffsignificant")),
             tabPanel("Complete",DT::dataTableOutput("DataDiffcomplete")),
             tabPanel("Up",DT::dataTableOutput("DataDiffup")),
             tabPanel("Down",DT::dataTableOutput("DataDiffdown"))
             
      )
    }
    
  })
  
  
  #####################
  ###
  ###
  ###################
  
  
  #### Export diff table
  output$exportDiffTable <- downloadHandler(
    filename = function() {
      extension='SHAMAN.csv'
      if(input$sepexpdiff == "\t") extension='SHAMAN.tsv' 
      paste(input$WhichExportTable,input$ContrastList_table,extension,sep="_")
      },
    
    content = function(file){
      switch(input$WhichExportTable,
             "Significant" = write.table(dataDiff()$significant, file,row.names = FALSE, sep=input$sepexpdiff),
             "Complete" = write.table(dataDiff()$complete, file,row.names = FALSE, sep=input$sepexpdiff),
             "Up" =  write.table(dataDiff()$up, file,row.names = FALSE, sep=input$sepexpdiff),
             "Down" =  write.table(dataDiff()$down, file,row.names = FALSE, sep=input$sepexpdiff)
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
    ChTM = "Error"
    target = dataInputTarget()$target
    labeled = dataInputTarget()$labeled
    CT = dataInput()$data$counts
    taxo = input$TaxoSelect
    VarInt = input$InterestVar
    
    ## Return NULL if there is no error
    if(!is.null(target) && length(VarInt)>=1) ChTM = CheckTargetModel(input,target,labeled,CT)$Error

    if(!is.null(target) && taxo!="..." && is.null(ChTM) && length(VarInt)>=1) res = actionButton("RunDESeq",strong("Run analysis"),icon = icon("caret-right"))
    
    return(res)
  })
  
  
  
  #####################################################
  ##
  ##                VISUALISATION
  ##
  #####################################################
  
  output$PlotVisuTree <- renderTreeWeightD3({
    resDiff = ResDiffAnal()
    taxo_table = dataInput()$data$taxo
    CT_Norm_OTU = dataMergeCounts()$CT_Norm
    res = NULL
    if(!is.null(resDiff$dds) && length(input$VisuVarInt)>=1) res = Plot_Visu_Tree(input,resDiff,CT_Norm_OTU,taxo_table)
    return(res)

  })

  
  output$PlotVisuBar <- renderChart({
    resDiff = ResDiffAnal()
    res = NULL
    tmp = NULL
    if(!is.null(resDiff$dds) && length(input$VisuVarInt)>=1) tmp = Plot_Visu_Barplot(input,resDiff)
    if(!is.null(tmp)) res = tmp$plotd3
    return(res)
  })
  
  
  output$heatmap <- renderD3heatmap({
    resDiff = ResDiffAnal()
    resplot = NULL
    if(!is.null(resDiff$dds)) resplot = withProgress(message="Loading...",Plot_Visu_Heatmap(input,resDiff))
    
    return(resplot)
  })
  
  
  output$heatmap_comp <- renderD3heatmap({
    resDiff = ResDiffAnal()
    ## Just for reactivity
    SelContrast = input$ContrastList_table_FC
    
    resplot = NULL
    filesize = file.info(namesfile)[,"size"]
    if(is.na(filesize)){filesize=0}
    if(filesize!=0)
    {
      BaseContrast = read.table(namesfile,header=TRUE)
      if(!is.null(resDiff$dds)) resplot = withProgress(message="Loading...",Plot_Visu_Heatmap_FC(input,BaseContrast,resDiff))
    }
    return(resplot)
  })
  
  
  
  
  output$ScatterplotD3 <- renderScatterD3({
    resDiff = ResDiffAnal()
    if(!is.null(resDiff$dds)) withProgress(message="Loading...",Plot_Visu_Scatterplot(input,resDiff))
  })
  
  output$Scatterplotgg <- renderPlot({
    resDiff = ResDiffAnal()
    if(!is.null(resDiff$dds)) withProgress(message="Loading...",Plot_Visu_Scatterplot(input,resDiff,lmEst=FALSE))
  })
  
  
  output$VennD3 <- renderD3vennR({
    resDiff = ResDiffAnal()
    ## Just for reactivity
    SelContrast = input$ContrastList_table_FC
    filesize = file.info(namesfile)[,"size"]
    if(is.na(filesize)){filesize=0}
    if(filesize!=0){
      BaseContrast = read.table(namesfile,header=TRUE)
      if(!is.null(resDiff$dds)) withProgress(message="Loading...",Plot_Visu_Venn(input,BaseContrast,resDiff))
    }
  })
  
  
  ## Regression coefficients Table
  output$lmRegScatter <- DT::renderDataTable(
    Plot_Visu_Scatterplot(input,ResDiffAnal(),lmEst=TRUE)$regCoef, 
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE, processing=FALSE
    ))
  
  ## Correlation coefficients Table
  output$CorTable <- DT::renderDataTable(
    Plot_Visu_Scatterplot(input,ResDiffAnal(),CorEst=TRUE)$cor.est, 
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE, processing=FALSE
    ))
  
  output$lmEquation <- renderPrint({ 
    res = Plot_Visu_Scatterplot(input,ResDiffAnal(),lmEst=TRUE)
    coef = res$regCoef
    Rsq = res$Rsq
    
    div(HTML(paste(h4(strong("Linear equation: ")),
                   "y =", round(coef[2,1],2),'x ',ifelse(coef[1,1]>=0,"+",""), round(coef[1,1],2),'<br/>','<br/>',
                   h4(strong("Adjusted R squared:")),round(Rsq,5)*100," %")))
  })
  
  
  
  # Infobox Contrast
  output$InfoSizeFactor <- renderPrint({
    input$RunDESeq
    res = div(HTML(""))
    
    tmpFact = SizeFactors_fromFile()
    
    if(!tmpFact$Check) res = div(HTML(paste('<p><font color="red">',tmpFact$Error,'</font></p>')))
    
    return(res)
  })
  
  
  output$Boxplot <- renderPlot({
    resDiff = ResDiffAnal()
    if(!is.null(resDiff$dds)) withProgress(message="Loading...",Plot_Visu_Boxplot(input,resDiff))
  },height=reactive(input$heightVisu))
  
  
  output$DiversityPlot <- renderPlot({
    resDiff = ResDiffAnal()
    if(!is.null(resDiff$dds)) withProgress(message="Loading...",Plot_Visu_Diversity(input,resDiff)$plot)
  })
  
  
  output$Diversitytable <- DT::renderDataTable(
    datatable({
    resDiff = ResDiffAnal()
    tmp = Plot_Visu_Diversity(input,resDiff)$dataDiv
    tmp$VarX=NULL; tmp$VarCol=NULL
    tmp[,c(4,5,1,2,3)]},rownames= FALSE),
  options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                 pageLength = 10,scrollX=TRUE, processing=FALSE
  ))
  
  ## Export Diversitytable in .csv
  output$ExportDiversitytable <- downloadHandler(
    filename = function() { 
      if(input$sepdiversity) 'SHAMAN_Diversity.tsv'
      else 'SHAMAN_Diversity.csv'
    },
    content = function(file){
      resDiff = ResDiffAnal()
      tmp = Plot_Visu_Diversity(input,resDiff)$dataDiv
      tmp$VarX=NULL; tmp$VarCol=NULL
      datatable(tmp[,c(4,5,1,2,3)],rownames= FALSE)
      write.table(tmp, file,row.names = FALSE, sep=input$sepdiversity)
    }
  )
  
  
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
  
  #### Select color and split for diversity
  
  output$SelectVarBoxDiv <- renderUI({
    
    selectVar = input$VisuVarInt
    selectInput("VarBoxDiv", h6(strong("Split by")),selectVar,selectVar[1],multiple = TRUE)
    
  })
  
  
  #   output$SelectVarDivCol <- renderUI({
  #     
  #     selectVar = input$VisuVarInt
  #     VarB = input$VarBoxDiv
  #     
  #     if(length(selectVar)>1) 
  #     {
  #       selectInput("VarDivCol", h6(strong("Color by")),c(NULL,selectVar[-which(selectVar%in%VarB)]),multiple = FALSE)
  #     }
  #     
  #   })
  #   
  
  
  
  output$plotVisu <- renderUI({
    
    res=NULL
    if(input$PlotVisuSelect=="Barplot") res =  showOutput("PlotVisuBar")
    if(input$PlotVisuSelect=="Heatmap") res =  d3heatmapOutput("heatmap", height = input$heightVisu+10)
    if(input$PlotVisuSelect=="Boxplot") res = plotOutput("Boxplot", height = input$heightVisu+10)
    if(input$PlotVisuSelect=="Tree") res = treeWeightD3Output('PlotVisuTree', height = input$heightVisu+10,width="100%")
    if(input$PlotVisuSelect=="Scatterplot" && !input$AddRegScatter) res = scatterD3Output("ScatterplotD3", height = input$heightVisu+10)
    if(input$PlotVisuSelect=="Scatterplot" && input$AddRegScatter) res = plotOutput("Scatterplotgg", height = input$heightVisu+10)
    
    if(input$PlotVisuSelect=="Diversity") res =  plotOutput("DiversityPlot", height = input$heightVisu+10)
    if(input$PlotVisuSelect=="Rarefaction") res = plotOutput("RarefactionPlot",dblclick = "RarefactionPlot_dblclick",brush = brushOpts(id = "RarefactionPlot_brush",resetOnNew = TRUE), height = input$heightVisu+10)
    return(res)
  })
  
  ## Comparison plots
  output$plotVisuComp <- renderUI({
    
    res=NULL
    if(input$PlotVisuSelectComp=="Heatmap_comp") res =  d3heatmapOutput("heatmap_comp", height = input$heightVisuComp+10)
    if(input$PlotVisuSelectComp=="Venn") res =  d3vennROutput("VennD3", height = input$heightVisuComp+10)
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
      if(input$SelectSpecifTaxo=="Diff" && length(input$ContrastList_table_Visu)>=1)
      {
        filesize = file.info(namesfile)[,"size"]
        if(is.na(filesize)){filesize=0}
        if(filesize!=0)
        { 
          BaseContrast = read.table(namesfile,header=TRUE)
          SelContrast = input$ContrastList_table_Visu
          #padj = Get_log2FC_padj(input,BaseContrast,resDiff, info = NULL)$padj
          selcontrast_matrix = as.matrix(BaseContrast[,SelContrast])
          colnames(selcontrast_matrix) = SelContrast
          padj = Get_log2FC_padj(input,selcontrast_matrix,resDiff, info = NULL)$padj
          Feature_names = rownames(padj)
          if(ncol(as.matrix(padj))>1)
          { 
            cont = which(colnames(padj)%in%SelContrast)
            padj = padj[,cont] 
          }
          if(ncol(padj)<2) ind = which(padj<=as.numeric(input$AlphaVal))
          if(ncol(padj) >= 2 && input$UnionInterContrasts=="Union") ind = which(apply(apply(padj,2,FUN = function(x){x<=as.numeric(input$AlphaVal)}),1,any))
          if(ncol(padj) >= 2 && input$UnionInterContrasts=="Inter") ind = which(!apply(!apply(padj,2,FUN = function(x){x<=as.numeric(input$AlphaVal)}),1,any))
          
          if(length(ind)>0 && !is.null(ind)) selTaxo = Feature_names[ind]
          else selTaxo = NULL
          res = selectizeInput("selectTaxoPlot",h6(strong(paste("Select the",input$TaxoSelect, "to plot"))),Available_taxo, selected = selTaxo,multiple = TRUE,options = list(minItems = 2))
        }    
      }
      if(input$SelectSpecifTaxo=="NoDiff" && length(input$ContrastList_table_Visu)>=1)
      {
        filesize = file.info(namesfile)[,"size"]
        if(is.na(filesize)){filesize=0}
        if(filesize!=0)
        { 
          BaseContrast = read.table(namesfile,header=TRUE)
          SelContrast = input$ContrastList_table_Visu
          #padj = Get_log2FC_padj(input,BaseContrast,resDiff, info = NULL)$padj
          selcontrast_matrix = as.matrix(BaseContrast[,SelContrast])
          colnames(selcontrast_matrix) = SelContrast
          padj = Get_log2FC_padj(input,selcontrast_matrix,resDiff, info = NULL)$padj
          Feature_names = rownames(padj)
          
          if(ncol(as.matrix(padj))>1)
          { 
            cont = which(colnames(padj)%in%SelContrast)
            padj = as.matrix(padj[,cont])
          }
          if(ncol(padj)<2) ind = which(padj>as.numeric(input$AlphaVal))
          if(ncol(padj) >= 2 && input$UnionInterContrasts=="Union") ind = which(apply(apply(padj,2,FUN = function(x){x>as.numeric(input$AlphaVal)}),1,any))
          if(ncol(padj) >= 2 && input$UnionInterContrasts=="Inter") ind = which(!apply(!apply(padj,2,FUN = function(x){x>as.numeric(input$AlphaVal)}),1,any))
          
          if(length(ind)>0 && !is.null(ind)) selTaxo = Feature_names[ind]
          else selTaxo = NULL
          res = selectizeInput("selectTaxoPlot",h6(strong(paste("Select the",input$TaxoSelect, "to plot"))),Available_taxo, selected = selTaxo,multiple = TRUE,options = list(minItems = 2))
        }    
      }
      if(input$SelectSpecifTaxo=="All") res = selectizeInput("selectTaxoPlot",h6(strong(paste("Select the",input$TaxoSelect, "to plot"))),Available_taxo, selected = Available_taxo,multiple = TRUE)
    }
    return(res)
  })
  
  
  ## For comp plot
  output$TaxoToPlotVisuComp <- renderUI({
    
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
      
      if(input$SelectSpecifTaxoComp=='Most')  res = selectizeInput("selectTaxoPlotComp",h6(strong(paste("Select the",input$TaxoSelect, "to plot"))),Available_taxo, selected = selTaxo,multiple = TRUE)
      if(input$SelectSpecifTaxoComp=="Diff" && length(input$ContrastList_table_VisuComp)>=1)
      {
        
        filesize = file.info(namesfile)[,"size"]
        if(is.na(filesize)){filesize=0}
        if(filesize!=0)
        { 
          BaseContrast = read.table(namesfile,header=TRUE)
          SelContrast = input$ContrastList_table_VisuComp
          #padj = Get_log2FC_padj(input,BaseContrast,resDiff, info = NULL)$padj
          selcontrast_matrix = as.matrix(BaseContrast[,SelContrast])
          colnames(selcontrast_matrix) = SelContrast
          padj = Get_log2FC_padj(input,selcontrast_matrix,resDiff, info = NULL)$padj
          Feature_names = rownames(padj)
          if(ncol(as.matrix(padj))>1)
          { 
            cont = which(colnames(padj)%in%SelContrast)
            padj = padj[,cont] 
          }
          if(ncol(padj)<2) ind = which(padj<=as.numeric(input$AlphaVal))
          if(ncol(padj) >= 2 && input$UnionInterContrastsComp=="Union") ind = which(apply(apply(padj,2,FUN = function(x){x<=as.numeric(input$AlphaVal)}),1,any))
          if(ncol(padj) >= 2 && input$UnionInterContrastsComp=="Inter") ind = which(!apply(!apply(padj,2,FUN = function(x){x<=as.numeric(input$AlphaVal)}),1,any))
          
          if(length(ind)>0 && !is.null(ind)) selTaxo = Feature_names[ind]
          else selTaxo = NULL
          res = selectizeInput("selectTaxoPlotComp",h6(strong(paste("Select the",input$TaxoSelect, "to plot"))),Available_taxo, selected = selTaxo,multiple = TRUE,options = list(minItems = 2))
        }    
      }
      if(input$SelectSpecifTaxoComp=="NoDiff"  && length(input$ContrastList_table_VisuComp)>=1)
      {
        filesize = file.info(namesfile)[,"size"]
        if(is.na(filesize)){filesize=0}
        if(filesize!=0)
        { 
          BaseContrast = read.table(namesfile,header=TRUE)
          SelContrast = input$ContrastList_table_VisuComp
          #padj = Get_log2FC_padj(input,BaseContrast,resDiff, info = NULL)$padj
          selcontrast_matrix = as.matrix(BaseContrast[,SelContrast])
          colnames(selcontrast_matrix) = SelContrast
          padj = Get_log2FC_padj(input,selcontrast_matrix,resDiff, info = NULL)$padj
          Feature_names = rownames(padj)
          
          if(ncol(as.matrix(padj))>1)
          { 
            cont = which(colnames(padj)%in%SelContrast)
            padj = padj[,cont] 
          }
          if(ncol(padj)<2) ind = which(padj>as.numeric(input$AlphaVal))
          if(ncol(padj) >= 2 && input$UnionInterContrastsComp=="Union") ind = which(apply(apply(padj,2,FUN = function(x){x>as.numeric(input$AlphaVal)}),1,any))
          if(ncol(padj) >= 2 && input$UnionInterContrastsComp=="Inter") ind = which(!apply(!apply(padj,2,FUN = function(x){x>as.numeric(input$AlphaVal)}),1,any))
          
          if(length(ind)>0 && !is.null(ind)) selTaxo = Feature_names[ind]
          else selTaxo = NULL
          res = selectizeInput("selectTaxoPlotComp",h6(strong(paste("Select the",input$TaxoSelect, "to plot"))),Available_taxo, selected = selTaxo,multiple = TRUE,options = list(minItems = 2))
        }    
      }
      if(input$SelectSpecifTaxoComp=="All") res = selectizeInput("selectTaxoPlotComp",h6(strong(paste("Select the",input$TaxoSelect, "to plot"))),Available_taxo, selected = Available_taxo,multiple = TRUE)
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
    
    target=dataInputTarget()$target
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
      ## Using list slows down the application if the number of rows too high
      if(nrow(counts)<300) 
      {
        Available_x = list(x1 = c(sort(rownames(counts))),"Diversity" = c("Alpha div","Shannon div","Inv.Simpson div","Simpson div"))
        names(Available_x)[1] = taxo
        if(any(numInd)) Available_x$Variables = namesTarget[numInd]
      } else{
        Available_x = c(sort(rownames(counts)),"Alpha div","Shannon div","Inv.Simpson div","Simpson div")
        if(any(numInd)) Available_x = c(Available_x,namesTarget[numInd])
      }
      Available_y = Available_x
    
      res[[1]] = selectizeInput("Xscatter",h6(strong("X variable")),Available_x, selected = Available_x[1],multiple = FALSE)
      res[[2]] = selectizeInput("Yscatter",h6(strong("Y variable")),Available_y, selected = Available_x[2],multiple = FALSE)
      res[[3]] = selectizeInput("ColorBy",h6(strong("Color variable")),c("None"="None",namesTarget[!numInd]),multiple = FALSE)
      res[[4]] = selectizeInput("PchBy",h6(strong("Symbol variable")),c("None"="None",namesTarget[!numInd]),multiple = FALSE)
      res[[5]] = selectizeInput("PointSize",h6(strong("Point size according to")),c("None"="None",Available_x), selected = NULL,multiple = FALSE)
    }
    
    return(res)
    
  })
  
  
  output$VarIntVisuTree <- renderUI({

    target=dataInputTarget()$target
    data = dataInput()$data
    taxo = input$TaxoSelect
    resDiff = ResDiffAnal()
    res = NULL

    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0 && !is.null(taxo) && taxo!="..." && !is.null(target))
    {
      counts = dataMergeCounts()$counts

      Available_x = sort(rownames(counts))

      res = selectizeInput("TaxoTree",h6(strong(paste("Select a specific",taxo,sep=" "))),c("...",Available_x),multiple = FALSE)

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
  
  
  #####################################################
  ##
  ##      Disable/Enable actions
  ##
  #####################################################
  
  
  ## Disable the actionbutton if the number of feature is lower than 2
  observe({
    
    input$TaxoSelect
    counts = dataMergeCounts()$counts
    ChTM = ""
    CT = dataInput()$data$counts
    target = dataInputTarget()$target
    labeled= dataInputTarget()$labeled
    ChTM = CheckTargetModel(input,target,labeled,CT)$Error
    
    if(input$AddFilter && !is.null(input$SliderThSamp) && !is.null(input$SliderThAb) && is.null(ChTM))
    {
      ind.filter =Filtered_feature(counts,input$SliderThSamp,input$SliderThAb)$ind
      counts = counts[-ind.filter,]
    }

    
    if (!is.null(counts)){
      if (nrow(counts)>=2){
        shinyjs::enable("RunDESeq")
      }
      if (nrow(counts)<2) {
        shinyjs::disable("RunDESeq")
      }
    } 
    if (is.null(counts) ) {
      shinyjs::disable("RunDESeq")
    }

  })
  
  
  observe({
    taxo = input$TaxoSelect
    target=dataInputTarget()$target
    
    if(is.null(target) || taxo =="...") 
    {
      shinyjs::disable("AddFilter")
    } else {
      shinyjs::enable("AddFilter")
    }
  })
  
  

  
})


