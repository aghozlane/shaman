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
  target = NULL
  
  
  ## JSON name for masque
  curdir  = getwd()
  json_name = tempfile(pattern = "file", tmpdir = paste(curdir,"www","masque","todo",sep= .Platform$file.sep),  fileext = ".json")
  
  ## Pass for MASQUE
  pass = gsub("file","",basename(file_path_sans_ext(json_name)))
  
  ## Popup messages
  observe(if(input$AddRegScatter) info("By adding the regression line, you will lose interactivity."))
  
  ## Reactive target
  values <- reactiveValues(TargetWorking = target,labeled=NULL,fastq_names_only=NULL,R1fastQ=NULL,R2fastQ=NULL,
                           json_name=json_name,num=0,pass=pass,login_email = NULL,is.valid =NULL,
                           biom_masque = NULL,tree_masque=NULL, masque_key = NULL, count_table_masque = NULL, 
                           rdp_annot_masque = NULL, rdp_thres_masque = NULL,
                           paths_fastq_tmp=NULL,curdir=curdir, error_progress=FALSE)
  
  ## Counts file
  dataInputCounts <-reactive({ 
    
    data = NULL
    inFile <- input$fileCounts
    if (is.null(inFile) && is.null(values$count_table_masque)) return(NULL)
    #if (is.null(inFile)) return(NULL)
    
    if (!is.null(values$count_table_masque) && file.exists(values$count_table_masque)){
      tryCatch(read.csv(values$count_table_masque,sep="\t",header=TRUE,check.names=FALSE)->data,
               error=function(e) sendSweetAlert(messageId="ErrorCounts",
                                                title = "Oops",
                                                text=paste("The count file can not be read in SHAMAN.\n \n",e),type ="error"))
    }
    else{
      tryCatch(read.csv(inFile$datapath,sep=input$sepcount,header=TRUE,check.names=FALSE)->data,
               error=function(e) sendSweetAlert(messageId="ErrorCounts",
                                                title = "Oops",
                                                text=paste("Your file can not be read in SHAMAN.\n \n",e),type ="error"))
    }
    #print(data)
    if(!is.null(data)){
      colnames(data) = gsub("-",".",colnames(data))
      ## Rownames
      if(!TRUE%in%duplicated(data[,1])) rownames(data)=data[,1];data=data[,-1]
      try(round(data, 0)->data, silent=T)
    }
    
    return(as.data.frame(data))
  })
  
  
  
  ## Taxo File
  dataInputTaxo <-reactive({ 
    
    inFile <- input$fileTaxo
    
    if (is.null(inFile) && is.null(values$rdp_annot_masque)) return(NULL)
    #if (is.null(inFile)) return(NULL)
    
    if(input$TypeTaxo=="Table" && !is.null(inFile)) 
    {
      tryCatch(read.csv(inFile$datapath,sep=input$septaxo,header=TRUE)->data,
               error=function(e) sendSweetAlert(messageId="ErrorTaxo",
                                                title = "Oops",
                                                text=paste("Your file can not be read in SHAMAN.\n \n",e),type ="error"))
      
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
    
    if(input$TypeTaxo=="RDP" && !is.null(inFile) || !is.null(values$rdp_annot_masque)) 
    {
      if (!is.null(values$rdp_annot_masque) && file.exists(values$rdp_annot_masque)){
        tryCatch(read_rdp(values$rdp_annot_masque,values$rdp_thres_masque)->data,
                 error=function(e) sendSweetAlert(messageId="ErrorRDP",
                                                  title = "Oops",
                                                  text=paste("The annotation file can not be read in SHAMAN.\n \n",e),type ="error"))
      }
      else{
        tryCatch(read_rdp(inFile$datapath,input$RDP_th)->data,
                 error=function(e) sendSweetAlert(messageId="ErrorRDP",
                                                  title = "Oops",
                                                  text=paste("Your file can not be read in SHAMAN.\n \n",e),type ="error"))
      }
      
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
    
    if (!is.null(inFile) && is.null(values$biom_masque)){
      tryCatch(read_biom(inFile$datapath)->data,
               error=function(e) sendSweetAlert(messageId="ErrorBiom1",
                                                title = "Oops",
                                                text=paste("Your file can not be read in SHAMAN.\n \n",e),type ="error"))
      
    }
    if (!is.null(values$biom_masque) && file.exists(values$biom_masque)){ 
      tryCatch(read_biom(values$biom_masque)->data,
               error=function(e) sendSweetAlert(messageId="ErrorBiom2",
                                                title = "Oops",
                                                text=paste("Your file can not be read in SHAMAN.\n \n",e),type ="error"))
    }
    return(data)
  })
  
  observeEvent(input$fileCounts,{
    values$count_table_masque=NULL;
  })
  observeEvent(input$fileTaxo,{
    values$rdp_annot_masque=NULL;
    values$rdp_thres_masque=NULL;
  })
  observeEvent(input$fileBiom,{
    values$biom_masque=NULL;
  })
  
  
  ## Unifrac File (tree)
  dataInputTree <-reactive({ 
    
    data = NULL
    inFile <- input$fileTree
    
    if (!is.null(inFile) && is.null(values$tree_masque)) {
      try(read.tree(inFile$datapath)->data, silent=T)
      CheckTree = CheckTreeFile(data)
      data = CheckTree$tree
      try(readLines(inFile$datapath)->treeseq, silent=T)
      return(list(data=data, Error=CheckTree$Error, Warning=CheckTree$Warning, treeseq=treeseq))
    }
    
    if (!is.null(values$tree_masque) && file.exists(values$tree_masque)) {
      try(read.tree(values$tree_masque)->data, silent=T)
      CheckTree = CheckTreeFile(data)
      data = CheckTree$tree
      try(readLines(values$tree_masque)->treeseq, silent=T)
      return(list(data=data, Error=CheckTree$Error, Warning=CheckTree$Warning, treeseq=treeseq))
    }
    
  })
  
  
  observeEvent(input$fileTree,{
    values$tree_masque=NULL;
  })
  
  
  # Infobox Tree (Unifrac)
  output$InfoTreePhylo_box <- renderInfoBox({
    input$fileTree
    tree_tmp = isolate(dataInputTree())
    tree = tree_tmp$data
    
    res = infoBox(h6(strong("Phylogenetic tree")), subtitle = h6(strong("Load the phylogenetic tree (optional)")) ,color = "light-blue",width=NULL,fill=TRUE, icon = icon("upload"))
    if(!is.null(tree)){
      if(!is.null(isolate(input$fileTree))){
        res = infoBox(h6(strong("Phylogenetic tree")), subtitle = h6("The phylogenetic has been loaded") ,color = "green",width=NULL,fill=TRUE, icon = icon("thumbs-o-up"))
        if(!is.null(tree_tmp$Warning)){      
          res = infoBox(h6(strong("Phylogenetic tree")), subtitle = h6(tree_tmp$Warning) ,color = "orange",width=NULL,fill=TRUE, icon = icon("warning"))
        }
        if(!is.null(tree_tmp$Error)){      
          res = infoBox(h6(strong("Phylogenetic tree")), subtitle = h6(tree_tmp$Error),color = "red",width=NULL,fill=TRUE, icon = icon("thumbs-o-down"))
        }
      }
    } 
    return(res)
  })
  
  observe({
    val <- input$annotationKingdomthreshold
    # Control the value, min, max, and step.
    # Step size is 2 when input value is even; 1 when value is odd.
    updateSliderInput(session, "annotationPhylumthreshold", value = input$annotationPhylumthreshold,
                      min = val, max = 1, step = 0.005)
  })
  observe({
    val <- input$annotationPhylumthreshold[2]
    # Control the value, min, max, and step.
    # Step size is 2 when input value is even; 1 when value is odd.
    updateSliderInput(session, "annotationClassthreshold", value = input$annotationClassthreshold,
                      min = val, max = 1, step = 0.005)
  })
  observe({
    val <- input$annotationClassthreshold[2]
    # Control the value, min, max, and step.
    # Step size is 2 when input value is even; 1 when value is odd.
    updateSliderInput(session, "annotationOrderthreshold", value = input$annotationOrderthreshold,
                      min = val, max = 1, step = 0.005)
  })
  observe({
    val <- input$annotationOrderthreshold[2]
    # Control the value, min, max, and step.
    # Step size is 2 when input value is even; 1 when value is odd.
    updateSliderInput(session, "annotationFamilythreshold", value = input$annotationFamilythreshold,
                      min = val, max = 1, step = 0.005)
  })
  observe({
    val <- input$annotationFamilythreshold[2]
    # Control the value, min, max, and step.
    # Step size is 2 when input value is even; 1 when value is odd.
    updateSliderInput(session, "annotationGenusthreshold", value = input$annotationGenusthreshold,
                      min = val, max = 1, step = 0.005)
  })
  observe({
    val <- input$annotationGenusthreshold[2]
    # Control the value, min, max, and step.
    # Step size is 2 when input value is even; 1 when value is odd.
    updateSliderInput(session, "annotationSpeciethreshold", value = input$annotationSpeciethreshold,
                      min = val, max = 1, step = 0.005)
  })
  ## Input data
  dataInput <-reactive({ 
    
    data = NULL
    check = NULL
    percent = NULL
    Taxo = NULL
    Counts = NULL
    inputData = NULL
    
    if(input$FileFormat=="fileCounts")
    {
      Counts = dataInputCounts()
      if(!input$NoTaxoFile) Taxo = dataInputTaxo()
      if(!is.null(Counts) && input$NoTaxoFile) {Taxo = data.frame(rownames(Counts),row.names = rownames(Counts));names(Taxo)=NA}
      
      if(!is.null(Counts) && !is.null(Taxo))
      { 
        tmp = GetDataFromCT(Counts,Taxo, ifelse(input$TypeTable=="MGS" && input$FileFormat!="fileBiom", TRUE, FALSE))
        data = list(counts=tmp$counts,taxo=tmp$taxo)
        ## Remove row with only O
        # data[["counts"]] = data[["counts"]][rowSums(data[["counts"]])>1,]
        
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
        ## Remove row with only O
        # data[["counts"]] = data[["counts"]][rowSums(data[["counts"]])>1,]
        
        check = list(CheckCounts=tmp$CheckCounts,CheckTaxo=tmp$CheckTaxo,CheckPercent=tmp$CheckPercent)
        percent = tmp$Percent
      }
    }
    
    
    
    #     if(input$FileFormat=="fileRData")
    #     {
    #       inFile <- input$fileRData
    #       load(inFile)
    #       if(!is.null(inputData)){
    #         data = inputData$data
    #         check = inputData$check
    #         percent = inputData$percent
    #       }
    #     }
    
    return(list(data=data,check=check,percent=percent))
  })
  
  
  ## Size factor file (optional)
  dataSizeFactors <-reactive({ 
    
    inFile <- input$fileSizeFactors
    
    if (is.null(inFile)) return(NULL)
    
    tryCatch(read.csv(inFile$datapath,sep=input$sepsize,header=TRUE)->data,
             error=function(e) sendSweetAlert(messageId="ErrorSizeFactor",
                                              title = "Oops",
                                              text=paste("Your file can not be read in SHAMAN.\n \n",e),type ="error"))
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
    ChMC = NULL
    data = isolate(dataInput()$data)
    target = isolate(values$TargetWorking)
    labeled= isolate(values$labeled)
    taxo = isolate(input$TaxoSelect)
    withProgress(
      if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0 && !is.null(taxo) && taxo!="..." && !is.null(target)) 
      {
        design = GetDesign(isolate(input),target)
        print(design)
        ChTM = CheckTargetModel(input,target,labeled,data$counts)$Error
        if(!is.null(design) && is.null(ChTM))
        {
          tmp = isolate(GetCountsMerge(input,data,taxo,target,design))
          #ChMC = tmp$Error
          #if (!is.null(ChMC))
          #{
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
          #}
        }
      }
      ,message="Merging the counts ...")
    return(list(counts=counts,CheckTarget=CheckTarget,normFactors=normFactors,CT_noNorm=CT_noNorm, CT_Norm=CT_Norm, Error = ChMC))
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
    res = shinydashboard::valueBox(paste0(0, "%"),h6(strong("Annotated features")), color = "light-blue",width=NULL,icon = icon("list"))
    
    if(cond)
    {
      percent = round(100*tmp$percent,2)
      if(percent==0) res = shinydashboard::valueBox(paste0(percent, "%"),h6(strong("Annotated features")), color = "red",width=NULL,icon = icon("list"))  
      if(percent!=0) res = shinydashboard::valueBox(paste0(percent, "%"),h6(strong("Annotated features")), color = "green",width=NULL,icon = icon("list"))  
      
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
                           tabName = "Visu")
                  #menuItem("Perspective plots", icon = icon("pie-chart"), tabName = "Krona")
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
    tree = dataInputTree()$data
    data=dataInput()$data
    res=NULL
    if(!is.null(tree))
    {
      res = tabBox(width = NULL, selected = "Count table",
                   tabPanel("Count table",DT::dataTableOutput("DataCounts")),
                   tabPanel("Taxonomy",DT::dataTableOutput("DataTaxo")),
                   tabPanel("Summary",h5(strong("Percentage of annotation")),htmlOutput("SummaryView"),
                            br(),h5(strong("Number of features by level:")),plotOutput("SummaryViewBarplot",width = 1200,height=500)),
                   tabPanel("Phylogeny", PhyloTreeMetaROutput('PhyloTreeMetaR'))
      )
      
    }
    else if(is.null(tree))
    {
      res = tabBox(width = NULL,selected = "Count table",
                   tabPanel("Count table",DT::dataTableOutput("DataCounts")),
                   tabPanel("Taxonomy",DT::dataTableOutput("DataTaxo")),
                   tabPanel("Summary",h5(strong("Percentage of annotation")),htmlOutput("SummaryView"),
                            br(),h5(strong("Number of features by level:")),plotOutput("SummaryViewBarplot",width = 1200,height=500))
      )
      
    }
    return(res)
  })
  
  observe({
    data=dataInput()$data
    if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0)
    {
      showElement("tabboxdata_col",anim=TRUE)
    } else hideElement("tabboxdata_col",anim=TRUE)
    
  })
  
  output$PhyloTreeMetaR <- renderPhyloTreeMetaR({
    PhyloTreeMetaR(dataInputTree()$treeseq,NULL)
  })
  
  output$SummaryView <- renderGvis({
    tmp = dataInput()
    data = tmp$data
    taxo = data$taxo
    counts = data$counts
    check = tmp$check
    cond = (!is.null(data$counts) && nrow(data$counts)>0 && !is.null(data$taxo) && nrow(data$taxo)>0 && is.null(check$CheckTaxo$Error) && is.null(check$CheckCounts$Error))
    
    res = NULL
    if(cond)
    {
      taxo = rbind(taxo,rep(NA,ncol(taxo)))
      #tmpPercent = round(apply(is.na(taxo),2,table)["FALSE",]/(nrow(taxo)-1)*100,2)
      
      tmp = apply(is.na(taxo),2,table)
      
      if (class(tmp) == "list") {
        tmp2 = sapply(tmp, function (x) {if (! "FALSE" %in% names(x)) {x["FALSE"] = 0} ; return(x["FALSE"])})
      }
      else
      {
        tmp2 = tmp["FALSE",]
      }
      
      tmpPercent = round(tmp2/(nrow(taxo)-1)*100,2)
      
      
      df <- data.frame(Label = colnames(taxo),Value = tmpPercent)
      
      # res = gvisGauge(df,options=list(min=0, max=100, greenFrom=80,
      #                                 greenTo=100, yellowFrom=60, yellowTo=80,
      #                                 redFrom=0, redTo=60, width=1200, height=300))
      res = gvisGauge(df,options=list(min=0, max=100, greenFrom=80,
                                      greenTo=100, yellowFrom=60, yellowTo=80,
                                      redFrom=0, redTo=60, width=800, height=200))
    }
    return(res)
  })
  
  
  output$SummaryViewBarplot <- renderPlot({
    tmp = dataInput()
    data = tmp$data
    taxo = data$taxo
    counts = data$counts
    check = tmp$check
    cond = (!is.null(data$counts) && nrow(data$counts)>0 && !is.null(data$taxo) && nrow(data$taxo)>0 && is.null(check$CheckTaxo$Error) && is.null(check$CheckCounts$Error))
    
    res = NULL
    if(cond)
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
  observe({ 
    
    inFile <- input$fileTarget
    counts = dataInput()$data$counts
    labeled = 0
    data = NULL
    
    if (is.null(inFile)) return(NULL)
    
    ## Read the data
    try(read.csv(inFile$datapath,sep=input$septarget,header=TRUE)->data,silent=TRUE)

    if(!is.null(data))
    {
      data = as.data.frame(data)
      names = colnames(data)
      
      ## Change the rownames
      if(!TRUE%in%duplicated(data[,1])) rownames(data)=gsub(pattern = "-",replacement = ".",as.character(data[,1]))
      
      ## Keep only the row which are in the count table
      ind = which(rownames(data)%in%colnames(counts))
      data = as.data.frame(data[ind,])
      colnames(data) = names
      
      
      ## Replace "-" by "."
      if(ncol(data)>1 && nrow(data)>1){
        ind_num = which(sapply(as.data.frame(data[,-1]),is.numeric)) + 1
        if(length(ind_num)>0){
          data_tmp =cbind( as.data.frame(apply(as.data.frame(data[,-ind_num]),2,gsub,pattern = "-",replacement = ".")),data[,ind_num])
          #data_tmp =cbind( as.data.frame(as.data.frame(data[,-ind_num])),data[,ind_num])
          colnames(data_tmp) = c(colnames(data)[-ind_num],colnames(data)[ind_num])
          data = data_tmp
        }
        if(length(ind_num)==0){data = as.data.frame(apply(data,2,gsub,pattern = "-",replacement = "."))}
      }
      
      values$TargetWorking = as.data.frame(data)
      #       ind_sel = Target_selection()
      #       if(length(ind))
      # target = as.data.frame(apply(target,2,gsub,pattern = "-",replacement = "."))
      
      #ord = order(rownames(data))
      #data = data[ord,]
      ### A SUPPRIMER 
      #rownames(data) <- colnames(counts)
      
      # Percent annotated
      #     print(ind)
      #     print(colnames(counts))
      #     print(rownames(data))
      values$labeled = length(ind)/length(colnames(counts))*100.0
    }
    
    # return(list(target = target, labeled=labeled))
  })
  
  
  
  
  #############################################################
  ##
  ##                        MASQUE
  ##
  #############################################################
  
  
  
  observeEvent(input$dir,{
    
    inFiles <- input$dir
    
    if (!is.null(inFiles)){
      # values$fastq_names_only = unique(paste(values$fastq_names_only,inFiles$name))
      values$paths_fastq_tmp = rbind(isolate(values$paths_fastq_tmp),inFiles)
      values$fastq_names_only = isolate(unique(values$paths_fastq_tmp[,"name"]))
    }
  })
  
  
  
  ## Create a fasta file containing the contaminant
  CreateFasta <- reactive({
    seq = NULL
    tmp = tempdir()
    fastaName = paste(tmp,paste(basename(file_path_sans_ext(json_name)),"_contaminant.fasta",sep=""),sep = .Platform$file.sep)
    
    if(!file.exists(fastaName)) file.create(fastaName,showWarnings=FALSE)
    if(input$PairedOrNot=="y"){seq =paste("#Seq1\n",input$R1primer,"\n \n","#Seq2\n",input$R2primer,sep="")}
    if(input$PairedOrNot=="n"){seq =input$primerSingle}
    if(!is.null(seq))  write(seq, file=fastaName)
    
  })
  
  
  ## Action with submit button
  MasqueSubmit <- eventReactive(input$submit,{
    #activate check_mail
    CMP = CheckMasque(input, values,check_mail = TRUE)
    Error = CMP$Error
    
    isJSONalreadyExist = file.exists(paste(values$curdir,"www","masque","doing",basename(json_name),sep= .Platform$file.sep))
    
    if(is.null(Error) && !isJSONalreadyExist)
    {
      CreateFasta()
      values$num = 1
      tmp = tempdir()
      # home <- normalizePath("~")
      home <- ""
      # path_glob = file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
      
      
      ## Paired-end
      if(input$PairedOrNot=="y"){
        cmp = 0
        nfiles = length(values$R1fastQ)+length(values$R2fastQ)
        
        withProgress(message = 'Uploading files...', value = 0, {
          pathToR1 = paste(tmp,"Masque_files_R1",sep=.Platform$file.sep)
          pathToR2 = paste(tmp,"Masque_files_R2",sep=.Platform$file.sep)
          
          if(dir.exists(pathToR1)){file.remove(list.files(pathToR1,full.names =TRUE))} else dir.create(pathToR1)
          if(dir.exists(pathToR2)){file.remove(list.files(pathToR2,full.names =TRUE))} else dir.create(pathToR2)
          for(i in values$R1fastQ){
            ind=which(i==values$paths_fastq_tmp[,"name"])[1]
            file.rename(from=values$paths_fastq_tmp[,"datapath"][ind], to=paste(tmp,"Masque_files_R1",i,sep= .Platform$file.sep))
            cmp = cmp +1
            incProgress(cmp/nfiles, detail = "Forward fastq files...")
          }
          for(i in values$R2fastQ){
            ind=which(i==values$paths_fastq_tmp[,"name"])[1]
            file.rename(from=values$paths_fastq_tmp[,"datapath"][ind], to=paste(tmp,"Masque_files_R2",i,sep= .Platform$file.sep))
            cmp = cmp +1
            incProgress(cmp/nfiles, detail = "Reverse fastq files...")
          }
        })
        
      } else{
        
        cmp = 0
        nfiles = length(values$fastq_names_only)
        
        withProgress(message = 'Uploading files...', value = 0, {
          
          pathTo = paste(tmp,"Masque_files",sep=.Platform$file.sep)
          
          if(dir.exists(pathTo)){file.remove(list.files(pathTo,full.names =TRUE))} else dir.create(pathTo)
          
          for(i in values$fastq_names_only){
            ind=which(i==values$paths_fastq_tmp[,"name"])[1]
            file.rename(from=values$paths_fastq_tmp[,"datapath"][ind], to=paste(tmp,"Masque_files",i,sep= .Platform$file.sep));cmp = cmp +1;incProgress(cmp/nfiles)}
        })
      }
      
      ## Create JSON file
      withProgress(message = 'Creating JSON file...',{CreateJSON(input,values)})
      if(file.exists(values$json_name)) values$num = 1
      sendSweetAlert(messageId="SuccessMasque",
                     title = "Success",
                     text = paste("Your data have been submitted. You will receive an e-mail once the computation over. <br /> This can take few hours.
                                  <br /> 
                                  <br /> 
                                  <br /> 
                                  <em> Remind: You can close shaman and use your key to check the progression and get your results: </em>",values$pass),
                     type = "success",
                     html=TRUE
                     )
    }
    
  })
  
  
  observeEvent(input$submit,{  
    
    tryCatch(MasqueSubmit(),
             error=function(e) sendSweetAlert(messageId="ErrorMasque",
                                              title = "Oops",
                                              text=paste("Something wrong when submitting.\n \n",e),type ="error"))
    
    
  },priority = 1)
  
  
  ## FastQ list
  output$FastQList_out <- renderUI({
    res = NULL
    if(!is.null(input$dir)){
      NullBox = h3(strong("0 FastQ file detected"),style="color:red;  text-align: center")
      res = NullBox
      
      if(length(values$fastq_names_only)>0)
      {
        res =list(selectInput("FastQList",label = "List of the fastq files in the selected directory",isolate(values$fastq_names_only),multiple =TRUE,selectize=FALSE,size = 6),
                  actionButton("RemoveFastQbut",'Remove file(s)',icon=icon("remove")))
      } else res = NullBox
    }
    return(res)
  })
  
  
  
  
  ## Remove FastQ function
  RemoveFastQ <-eventReactive(input$RemoveFastQbut,{
    
    if(length(input$FastQList)>0)
    {
      ind = which(values$fastq_names_only%in% input$FastQList)
      values$fastq_names_only = values$fastq_names_only[-ind]
      values$paths_fastq_tmp = values$paths_fastq_tmp[-ind,]
      updateSelectInput(session, "FastQList","List of the fastq files in the selected directory",values$fastq_names_only)
    }
  })
  
  
  
  
  
  ## Remove FastQ
  observeEvent(input$RemoveFastQbut,{  
    
    RemoveFastQ()
    if(input$MatchFiles_button>=1) MatchFiles()
    
  },priority=1)
  
  
  
  ## Update R1 and R2 lists
  MatchFiles <-reactive({
    
    if(length(values$fastq_names_only)>0 && input$PairedOrNot=='y')
    {
      indR1 = grep(input$R1files,values$fastq_names_only)
      indR2 = grep(input$R2files,values$fastq_names_only)
      
      ## If some are R1 and R2, removed from both list
      b12 = intersect(indR1,indR2)
      if(length(b12)>0) {indR1 = indR1[-which(indR1%in%b12)]; indR2 = indR2[-which(indR2%in%b12)]}
      
      
      if(length(indR1)>0 && length(indR2)>0){
        
        values$R1fastQ = values$fastq_names_only[indR1]
        values$R2fastQ = values$fastq_names_only[indR2]
        
        ## If some are only R1 or R2
        tmpR1 = gsub(input$R1files,x=values$R1fastQ,""); tmpR2 = gsub(input$R2files,x=values$R2fastQ,"")
        values$R1fastQ = values$R1fastQ[tmpR1%in%tmpR2];values$R2fastQ = values$R2fastQ[tmpR2%in%tmpR1]
        
        ## Update the files lists
        updateSelectInput(session, "R1filesList","",values$R1fastQ); updateSelectInput(session, "R2filesList","",values$R2fastQ)
      } else{updateSelectInput(session, "R1filesList","","");updateSelectInput(session, "R2filesList","","")}
    } else{updateSelectInput(session, "R1filesList","","");updateSelectInput(session, "R2filesList","","")}
  })
  
  
  observeEvent(input$MatchFiles_button,{  
    
    MatchFiles()
    
  },priority=1)
  
  
  
  ## Remove FastQ function directly from R1, R2
  RemoveFastQ_R1R2 <-eventReactive(input$RemoveFastQbut_R1R2,{
    
    if(length(input$R1filesList)>0)
    {
      ind = which(values$R1fastQ%in% input$R1filesList)
      values$R1fastQ = values$R1fastQ[-ind]
      updateSelectInput(session, "R1filesList","",values$R1fastQ)
    }
    
    if(length(input$R2filesList)>0)
    {
      ind = which(values$R2fastQ%in% input$R2filesList)
      values$R2fastQ = values$R2fastQ[-ind]
      updateSelectInput(session, "R2filesList","",values$R2fastQ)
    }
    
    
  })
  
  
  
  RemoveFastQ_R1R2_all <-eventReactive(input$LoadFiles,{
    
    values$R1fastQ = NULL
    updateSelectInput(session, "R1filesList","","")
    values$R2fastQ = NULL
    updateSelectInput(session, "R2filesList","","")
    
  })
  
  ## Remove FastQ from R1, R2
  observeEvent(input$RemoveFastQbut_R1R2,{  
    
    RemoveFastQ_R1R2()
    
  },priority=1)
  
  
  ## Remove FastQ from R1, R2 (load button)
  observeEvent(input$LoadFiles,{  
    
    RemoveFastQ_R1R2_all()
    
  },priority=1)
  
  
  
  # observe({
  #   CMP = CheckMasque(input, values)
  #   toggleState("submit",condition = is.null(CMP$Error))
  # })
  
  observe({
    toggleState("box-match",condition = (input$PairedOrNot=="y"))
  })
  
  
  output$InfoMasque<- renderUI({
    input$submit
    
    CMP = isolate(CheckMasque(input, values, check_mail = FALSE))
    
    if(!is.null(CMP$Error) && input$submit>0) {
      toastr_error(title="Error",message=HTML(CMP$Error),closeButton = TRUE,position ="bottom-right",preventDuplicates = TRUE,newestOnTop = TRUE,
                   progressBar = FALSE,showDuration = 300,showMethod="show",timeOut = 10000)
    }
    
  })
  
  
  output$InfoMasqueHowTo<- renderUI({
    input$submit
    
    CMP = isolate(CheckMasque(input, values, check_mail = FALSE))
    
    if(!is.null(CMP$HowTo) && input$submit>0) {
      toastr_success(title="How to",message=HTML(CMP$HowTo),closeButton = TRUE,position ="bottom-right",preventDuplicates = TRUE,newestOnTop = TRUE,
                     progressBar = FALSE,showDuration = 300,showMethod="show",timeOut = 10000)
    }
    
  })
  
  
  
  # observeEvent(input$submit,{
  #   
  #   CMP = isolate(CheckMasque(input, values))
  #   if(!is.null(CMP$HowTo)) { 
  #     toastr_success(title="How to",message=HTML(CMP$HowTo),closeButton = TRUE,position ="bottom-right",preventDuplicates = TRUE,newestOnTop = TRUE,
  #                  progressBar = FALSE,showDuration = 300,showMethod="show",timeOut = 10000)
  #   }
  #   
  # })
  # 
  # 
  # observeEvent(input$submit,{
  #   
  #   CMP = isolate(CheckMasque(input, values))
  #   if(!is.null(CMP$Error)) { 
  #     toastr_error(title="Error",message=HTML(CMP$Error),closeButton = TRUE,position ="bottom-right",preventDuplicates = TRUE,newestOnTop = TRUE,
  #                  progressBar = FALSE,showDuration = 300,showMethod="show",timeOut = 10000)
  #   }
  #   
  # })
  
  #########    ICONS   ################ 
  
  
  output$spinner_anim <- renderUI(
    htmltools::HTML('<i class="fa fa-spinner fa-pulse fa-fw" style="color:white" ></i><span class="sr-only">Loading...</span>')
  )
  
  
  output$spinner_icon <- renderUI(
    htmltools::HTML('<i class="fa fa-spinner" aria-hidden="true" style="color:white" ></i><span class="sr-only">Loading...</span>')
  )
  
  output$pause_icon <- renderUI(
    htmltools::HTML('<i class="fa fa-pause" aria-hidden="true" style="color:white" ></i><span class="sr-only">Loading...</span>')
  )
  
  output$check_icon <- renderUI(
    htmltools::HTML('<i class="fa fa-check" style="color:white" ></i><span class="sr-only">Loading...</span>')
  )
  
  
  output$key_icon <- renderUI(
    htmltools::HTML('<i class="fa fa-key" style="color:white" ></i><span class="sr-only">Loading...</span>')
  )
  
  output$test_icon <- renderUI(
    htmltools::HTML('<img src="icon.png" alt="dna" style="width:80px;height:80px;">')
  )
  
  output$amplicon_icon <- renderUI(
    htmltools::HTML('<img src="icons/amplicon.png" alt="dna" style="width:80px;height:80px;">')
  )
  output$dereplication_icon <- renderUI(
    htmltools::HTML('<img src="icons/dereplication.png" alt="dna" style="width:80px;height:80px;">')
  )
  
  output$singleton_icon <- renderUI(
    htmltools::HTML('<img src="icons/singleton.png" alt="dna" style="width:80px;height:80px;">')
  )
  output$chimera_icon <- renderUI(
    htmltools::HTML('<img src="icons/chimera.png" alt="dna" style="width:80px;height:80px;">')
  )
  output$otu_icon <- renderUI(
    htmltools::HTML('<img src="icons/otu.png" alt="dna" style="width:80px;height:80px;">')
  )
  output$silva_icon <- renderUI(
    htmltools::HTML('<img src="icons/silva.png" alt="dna" style="width:199px;height:80px;">')
  )
  output$greengenes_icon <- renderUI(
    htmltools::HTML('<img src="icons/greengenes.png" alt="dna" style="width:143px;height:80px;">')
  )
  output$rdp_icon <- renderUI(
    htmltools::HTML('<img src="icons/rdp.png" alt="dna" style="width:107px;height:80px;">')
  )
  output$findley_icon <- renderUI(
    htmltools::HTML('<img src="icons/findley.png" alt="dna" style="width:131px;height:80px;">')
  )
  output$unite_icon <- renderUI(
    htmltools::HTML('<img src="icons/unite.png" alt="dna" style="width:163px;height:80px;">')
  )
  output$underhill_icon <- renderUI(
    htmltools::HTML('<img src="icons/underhill.png" alt="dna" style="width:80px;height:80px;">')
  )
  
  #####################################
  
  
  
  output$progressBoxMasque <- renderValueBox({
    res = NULL
    error_progress = values$error_progress
    num = round(as.numeric(values$num),1)
    res = shinydashboard::valueBox("0 %",h6(strong("Waiting for the data...")), color = "light-blue",width=NULL,icon = uiOutput("spinner_icon"))  
    CMP = isolate(CheckMasque(input, values, check_mail=FALSE))
    Error = CMP$Error
    if(num>=100) res = shinydashboard::valueBox(paste("100 %"),h6(strong("Analysis completed ! Check your mail.")), color = "green",width=NULL,icon =  uiOutput("check_icon"))
    else if(is.null(Error) || num>=1) res = shinydashboard::valueBox(paste(values$num,"%"),h6(strong("Analysis in progress...")), color = "green",width=NULL,icon = uiOutput("spinner_anim"))  
    else if(!is.null(Error) && num>=1 || error_progress) res = shinydashboard::valueBox(paste(values$num,"%"),h6(strong("Workflow failed during progression...")), color = "red",width=NULL,icon = uiOutput("spinner_anim"))  
    return(res)
  })
  
  
  # output$infoBoxPass <- renderInfoBox({
  #   
  #   res = NULL
  #   # pass = toupper(gsub(" ","",input$password))
  #   # passOK = identical(pass,toupper(values$pass))
  #   # 
  #   if(input$password =="") res = infoBox("Get a key","Require a valid email address", color = "light-blue",width=NULL,icon = uiOutput("key_icon"),fill = TRUE)
  #   
  #   if(passOK)  res = infoBox("Key created !",paste("Your key is ",values$pass), color = "green",width=NULL,icon = uiOutput("key_icon"),fill = TRUE)
  #   if(!passOK && input$password !="")  res = infoBox("Invalid key","Use the key that you have received by mail", color = "red",width=NULL,icon = uiOutput("key_icon"),fill = TRUE)
  #   return(res)
  # })
  
  
  
  output$infoBoxPass <- renderInfoBox({
    
    res = NULL
    # pass = toupper(gsub(" ","",input$password))
    # passOK = identical(pass,toupper(values$pass))
    # 
    res = infoBox("Get a key","Require a valid email address", color = "light-blue",width=NULL,icon = uiOutput("key_icon"),fill = TRUE)
    
    if(input$checkMail>=1 && isValidEmail(input$to))  res = infoBox("Key created !",paste("Your key is ",values$pass), color = "green",width=NULL,icon = uiOutput("key_icon"),fill = TRUE)
    if(input$checkMail>=1 && !isValidEmail(input$to))  res = infoBox("Invalid email","Enter a valid email address to get your key", color = "red",width=NULL,icon = uiOutput("key_icon"),fill = TRUE)
    return(res)
  })
  
  
  
  output$infoBoxFastQ <- renderInfoBox({
    # FastqLoad()
    res = NULL
    res = infoBox("Fastq files","Load the fastq files ", color = "light-blue",width=NULL,icon = icon("play"),fill = TRUE)
    
    if(!is.null(input$dir)){
      if(length(unique(values$fastq_names_only))==0) res = infoBox("Fastq files","Select at least one fastq file", color = "red",width=NULL,icon = icon("play"),fill = TRUE)
      if(length(unique(values$fastq_names_only))>0) res = infoBox("Fastq files",paste(length(unique(values$fastq_names_only)), "files are loaded"), color = "green",width=NULL,icon = icon("play"),fill = TRUE)
    }
    return(res)
  })
  
  
  output$infoBoxFastQ_match <- renderInfoBox({
    
    res = NULL
    SM = SamplesMasque(input,values)
    if(input$PairedOrNot=="n"){  res = infoBox("Match the pairs","Only for paired-end sequencing", color = "black",width=NULL,icon = icon("exchange"),fill = TRUE)}
    
    if(input$PairedOrNot=='y'){
      if(input$MatchFiles_button==0) res = infoBox("Match the pairs","Identify forward and reverse files and then click the match button", color = "light-blue",width=NULL,icon = icon("exchange"),fill = TRUE)
      if(input$MatchFiles_button>0 && length(SM$samples)>=1){res = infoBox("Pairs are matched",paste(length(SM$samples), "samples are detected"), color = "green",width=NULL,icon = icon("exchange"),fill = TRUE)}
      if(input$MatchFiles_button>0 && length(SM$samples)<1){res = infoBox("Match the pairs","Failed. 0 samples detected", color = "red",width=NULL,icon = icon("exchange"),fill = TRUE)}
    }
    
    return(res)
  })
  
  
  
  
  ## plot gauge
  # output$gaugeMasque <-renderGauge({
  #   input$submit
  #   
  #   res = NULL;
  #   num = round(as.numeric(values$num),1)
  #   
  #   CMP = isolate(CheckMasque(input, values))
  #   Error = CMP$Error
  #   if(is.null(Error) || num>1) res = gauge(min(num,100), 0,100,symbol = '%',label= "Progress...")
  # 
  #   return(res)
  # })
  # 
  
  ## Timer for the gauge
  Timer <- reactiveTimer(10000)
  
  ## Check masque progress
  observe({
    
    Timer()
    CMP = isolate(CheckMasque(input, values, check_mail=FALSE))
    Error = CMP$Error
    
    if(is.null(Error) && isolate(values$num)<100){
      
      progress_file = paste(values$curdir,"www","masque","doing",paste(basename(file_path_sans_ext(json_name)),"_progress",".txt",sep=""),sep= .Platform$file.sep)
      if(file.exists(progress_file))
      {
        pf = read_lines(progress_file)
        if(!is.null(pf)){
          pf = as.numeric(pf)
          if(!is.na(pf)){
            pf = min(pf,100); pf = max(pf,0)
            if(isolate(values$num)<pf) {values$num = round(pf,1)}
          }
        }
      }
    }
    #CHANGEMENT DE COULEUR BORDEL
    else if(!is.null(Error) && isolate(values$num)>=1) values$error_progress = TRUE
  })
  
  
  observe({
    toggleState("checkMail",condition = isValidEmail(input$to))
  })
  
  
  
  
  # 
  # Project_status <-reactive({
  #   input$Check_project_over
  #   input$Check_project
  #   
  #   passOK = FALSE;status = NULL;file = NULL
  #  print("OK")
  #   json_files = list.files(paste(values$curdir,"www","masque",sep= .Platform$file.sep),pattern = "json",recursive = TRUE)
  #   allpass = gsub(gsub(json_files,pattern = ".*file",replacement = ""),pattern = ".json",replacement = "")
  #   
  #   print(allpass)
  #   
  #   if(length(allpass)>0){
  #     passOK = any(isolate(input$password)==allpass)
  #     if(passOK){
  #       ind = which(isolate(input$password)==allpass)
  #       file = paste(values$curdir,"www","masque",json_files[ind],sep= .Platform$file.sep)
  #       status = gsub(json_files[ind],pattern = "/.*",replacement = "")
  #     }
  #   }
  #   
  #   return(list(status=status,file=file,passOK=passOK))
  # })
  
  
  # 
  # Project_current <-reactive({ 
  #   input$Check_project_over
  #   
  #   passOK = FALSE;status = NULL;file = NULL
  #   
  #   json_files = list.files(paste(values$curdir,"www","masque",sep= .Platform$file.sep),pattern = "json",recursive = TRUE)
  #   allpass = gsub(gsub(json_files,pattern = ".*file",replacement = ""),pattern = ".json",replacement = "")
  #   
  #   print(allpass)
  #   
  #   if(length(allpass)>0){
  #     passOK = any(isolate(values$pass)==allpass)
  #     if(passOK){
  #       ind = which(values$pass==allpass)
  #       file = paste(values$curdir,"www","masque",json_files[ind],sep= .Platform$file.sep)
  #       status = gsub(json_files[ind],pattern = "/.*",replacement = "")
  #     }
  #   }
  #   
  #   return(list(status=status,file=file,passOK=passOK))
  # })
  # 
  
  
  observeEvent(input$Check_project,{
    values$masque_key = input$password
    PS = Project_status(values$masque_key,values$curdir)
    #resbox = Project_box_result(values$masque_key,values$curdir)
    #PS = resbox$PS
    #passOK = PS$passOK
    #if(!is.null(input$password) && input$password!="" && !passOK){
    if(!is.null(input$password) && input$password!="" && !PS$passOK){ 
      removeCssClass(class = 'pwdGREEN', selector = '#password')
      addCssClass(class = 'pwdRED', selector = '#password')
    }
    
    #if(!is.null(input$password) && input$password!="" && passOK){  
    if(!is.null(input$password) && input$password!="" && PS$passOK){
      removeCssClass(class = 'pwdRED', selector = '#password')
      addCssClass(class = 'pwdGREEN', selector = '#password')
      hideElement("masque-form",anim=TRUE)
      hideElement("masque-infobox",anim=TRUE)
      hideElement("boxsum",anim=TRUE)
      showElement("reload-project",anim=TRUE)
      hideElement("project_over",anim=TRUE)
      hideElement("pass",anim=TRUE)
      #showElement("MasqueToShaman",anim=TRUE)
    }
    if(is.null(input$password) || input$password==""){    
      removeCssClass(class = 'pwdRED', selector = '#password')
      removeCssClass(class = 'pwdGREEN', selector = '#password')
    }
  })
  
  ## the same from home
  observeEvent(input$Check_project_home,{
    values$masque_key = input$password_home
    
    PS = Project_status(values$masque_key,values$curdir)
    
    if(!is.null(input$password_home) && input$password_home!="" && !PS$passOK){ 
      removeCssClass(class = 'pwdGREEN', selector = '#password_home')
      addCssClass(class = 'pwdRED', selector = '#password_home')
    }
    
    if(!is.null(input$password_home) && input$password_home!="" && PS$passOK){
      removeCssClass(class = 'pwdRED', selector = '#password_home')
      addCssClass(class = 'pwdGREEN', selector = '#password_home')
      hideElement("masque-form",anim=TRUE)
      hideElement("masque-infobox",anim=TRUE)
      hideElement("boxsum",anim=TRUE)
      showElement("reload-project",anim=TRUE)
      hideElement("project_over",anim=TRUE)
      hideElement("pass",anim=TRUE)
      #showElement("MasqueToShaman",anim=TRUE)
    }
    if(is.null(input$password_home) || input$password_home==""){    
      removeCssClass(class = 'pwdRED', selector = '#password_home')
      removeCssClass(class = 'pwdGREEN', selector = '#password_home')
    }
  })
  
  
  #eventReactive(input$refresh, {
  #  shinyjs::js$refresh()
  #hideElement("gaugeMasque_progress")
  #reset("gaugeMasque_progress")
  #showElement("gaugeMasque_progress")
  #shinyjs::reset("gaugeMasque_progress")
  #shinyjs::reset("reload-project")
  #showElement("reload-project")
  #})
  
  observeEvent(input$Check_project_over,{
    values$masque_key = values$pass
    PS = Project_status(values$masque_key,values$curdir)
    #resbox = Project_box_result(values$masque_key,values$curdir)
    #PS = resbox$PS
    #passOK = PS$passOK
    
    if(PS$passOK){  
      hideElement("masque-form",anim=TRUE)
      hideElement("masque-infobox",anim=TRUE)
      hideElement("boxsum",anim=TRUE)
      showElement("reload-project",anim=TRUE)
      hideElement("project_over",anim=TRUE)
      #showElement("MasqueToShaman",anim=TRUE)
      ##TODO test
      hideElement("pass",anim=TRUE)
    }
  })
  
  
  
  
  
  ### Check button once computation are over
  observe({
    if(values$num>=100){
      hideElement("masque-form",anim=TRUE)
      hideElement("boxsum",anim=TRUE)
      hideElement("reload-project",anim=TRUE)
      showElement("project_over",anim=TRUE)
      showElement("current-project",anim=TRUE)
    }
    
  })
  
  
  # output$masque_results<- renderUI({
  #   
  #   res=NULL
  #   PS = Project_current()
  #   
  #   if(PS$status=="done")
  #   {
  #     hideElement("project_over",anim=TRUE)
  #     showElement("MasqueToShaman",anim=TRUE)
  #     json_file = PS$file
  #     folder_name = basename(file_path_sans_ext(json_file))
  #     print(folder_name)
  #     ### Paste file name as folder
  #     annot_process = paste(values$curdir,"www","masque","done",folder_name,"shaman_annotation_process.tsv",sep= .Platform$file.sep)
  #     
  #     if(file.exists(annot_process))
  #     {
  #       ap = read.csv(annot_process,sep="\t")
  #       res = fluidRow(
  #         HTML('<center><h1><strong>Your project is done !</strong></h1> <br/> <em><h4> Hereafter is a summary of the building and annotation processes</h4> </em> </center>'),
  #         br(),
  #         column(width=4,
  #                valueBox(ap$Count[1],tags$strong(tags$h5("Number of amplicons", style = "width: 70%;")), color = "light-blue",width=NULL,icon = uiOutput("amplicon_icon")),
  #                valueBox(ap$Count[2],tags$strong(tags$h5("Remaining amplicons after dereplication", style = "width: 70%;")), color = "light-blue",width=NULL,icon = uiOutput("dereplication_icon")),
  #                valueBox(ap$Count[3],tags$strong(tags$h5("Remaining amplicons after removing singletons", style = "width: 70%;")), color = "light-blue",width=NULL,icon = uiOutput("singleton_icon")),
  #                valueBox(ap$Count[4],tags$strong(tags$h5("Remaining amplicons after removing chimera", style = "width: 70%;")), color = "light-blue",width=NULL,icon = uiOutput("chimera_icon"))
  #         )
  #       )
  #     } else{res =HTML('<center><h1><strong>Your project is done !</strong></h1> <br/> <em><h4> Hereafter is a summary of the building and annotation processes</h4> </em> </center>')}
  #   }
  #   
  #   if(PS$status=="doing"){
  #     hideElement("project_over",anim=TRUE)
  #     res = fluidRow(
  #       HTML('<center><h1><strong>Your project is currently running !</strong></h1> <br/> <br/> </center>'),
  #       inlineCSS(gaugeCSS),
  #       gaugeOutput("gaugeMasque_progress", width = "100%", height = "100%")
  #     )
  #   }
  #   
  #   if(PS$status=="error"){
  #     hideElement("project_over",anim=TRUE)
  #     
  #     json_file = PS$file
  #     error_file = paste(values$curdir,"www","masque","error",paste(basename(file_path_sans_ext(json_file)),"_error",".txt",sep=""),sep= .Platform$file.sep)
  #     print(error_file)
  #     if(file.exists(error_file)){error_message = read_lines(error_file)}
  #     
  #     res = fluidRow(
  #       HTML('<center><h1><strong>Sorry, the workflow failed during progression</strong></h1> <br/> <em><h4> Hereafter is the message error.</h4> </em> <br/> </center>'),
  #       
  #       column(width = 12,
  #              div(style = "background-color: white; margin: 0 auto;width: 50%; text-align:center;border:1px solid red",
  #                  h4(strong("Error message")),
  #                  hr(style = "width: 70%;"),
  #                  HTML(paste(error_message,collapse = " <br/> ")),
  #                  br()
  #              )
  #       )
  #     )
  #   }
  #   
  #   return(res)
  # })
  # 
  # 
  
  ## Action of the comeback button
  observeEvent(input$comeback,{
    showElement("masque-infobox",anim=TRUE)
    showElement("masque-form",anim=TRUE)
    showElement("boxsum",anim=TRUE)
    showElement("pass",anim=TRUE)
    hideElement("reload-project",anim=TRUE)
    hideElement("project-over-wait",anim=TRUE)
    hideElement("project_over",anim=TRUE)
  })
  
  
  # observeEvent(input$Check_project_over,{
  #   values$masque_key = values$pass
  #   print(values$masque_key)
  #   # hideElement("project_over",anim=TRUE)
  # })
  # 
  # observeEvent(input$Check_project,{
  #   values$masque_key = isolate(input$password)
  #   print(values$masque_key)
  #   # hideElement("project_over",anim=TRUE)
  # })
  # 
  output$masque_status_key <-renderUI({
    input$refresh
    res = NULL
    resbox = Project_box_result(values$masque_key,values$curdir)
    
    # if(resbox$PS$status=='done') showElement("MasqueToShaman",anim=TRUE)
    # if(resbox$PS$status!='done') hideElement("MasqueToShaman",anim=TRUE)
    return(resbox$box)
    
    # res=NULL
    # PS = Project_status()
    # print("OK2")
    # if(PS$status=="done")
    # {
    #   showElement("MasqueToShaman",anim=TRUE)
    #   json_file = PS$file
    #   folder_name = basename(file_path_sans_ext(json_file))
    #   print(folder_name)
    #   ### Paste file name as folder
    #   annot_process = paste(values$curdir,"www","masque","done",folder_name,"shaman_annotation_process.tsv",sep= .Platform$file.sep)
    #   
    #   if(file.exists(annot_process))
    #   {
    #     ap = read.csv(annot_process,sep="\t")
    #     res = fluidRow(
    #         HTML('<center><h1><strong>Your project is done !</strong></h1> <br/> <em><h4> Hereafter is a summary of the building and annotation processes</h4> </em> </center>'),
    #         br(),
    #         column(width=4,
    #                 valueBox(ap$Count[1],tags$strong(tags$h5("Number of amplicons", style = "width: 70%;")), color = "light-blue",width=NULL,icon = uiOutput("amplicon_icon")),
    #                 valueBox(ap$Count[2],tags$strong(tags$h5("Remaining amplicons after dereplication", style = "width: 70%;")), color = "light-blue",width=NULL,icon = uiOutput("dereplication_icon")),
    #                 valueBox(ap$Count[3],tags$strong(tags$h5("Remaining amplicons after removing singletons", style = "width: 70%;")), color = "light-blue",width=NULL,icon = uiOutput("singleton_icon")),
    #                 valueBox(ap$Count[4],tags$strong(tags$h5("Remaining amplicons after removing chimera", style = "width: 70%;")), color = "light-blue",width=NULL,icon = uiOutput("chimera_icon"))
    #         )
    #       )
    #   } else{res =HTML('<center><h1><strong>Your project is done !</strong></h1> <br/> <em><h4> Hereafter is a summary of the building and annotation processes</h4> </em> </center>')}
    # }
    # 
    # if(PS$status=="doing"){
    #       res = fluidRow(
    #             HTML('<center><h1><strong>Your project is currently running !</strong></h1> <br/> <br/> </center>'),
    #             inlineCSS(gaugeCSS),
    #             gaugeOutput("gaugeMasque_progress", width = "100%", height = "100%")
    #       )
    # }
    # 
    # if(PS$status=="error"){
    #   
    #   json_file = PS$file
    #   error_file = paste(values$curdir,"www","masque","error",paste(basename(file_path_sans_ext(json_file)),"_error",".txt",sep=""),sep= .Platform$file.sep)
    #   print(error_file)
    #   if(file.exists(error_file)){error_message = read_lines(error_file)}
    #     
    #   res = fluidRow(
    #     HTML('<center><h1><strong>Sorry, the workflow failed during progression</strong></h1> <br/> <em><h4> Hereafter is the message error.</h4> </em> <br/> </center>'),
    #     
    #     column(width = 12,
    #            div(style = "background-color: white; margin: 0 auto;width: 50%; text-align:center;border:1px solid red",
    #                h4(strong("Error message")),
    #                hr(style = "width: 70%;"),
    #               HTML(paste(error_message,collapse = " <br/> ")),
    #               br()
    #            )
    #     )
    #   )
    # }
    # 
    # return(res)
    
  })
  
  
  output$gaugeMasque_progress <- renderGauge({
    res = NULL
    PS = Project_status(values$masque_key,values$curdir)
    if(PS$passOK){
      
      if(PS$status=="doing"){
        json_file = PS$file
        progress_file = paste(values$curdir,"www","masque","doing",paste(basename(file_path_sans_ext(json_file)),"_progress",".txt",sep=""),sep= .Platform$file.sep)
        if(file.exists(progress_file))
        {
          pf = read_lines(progress_file)
          pf = round(as.numeric(pf),1)
          if(!is.na(pf)){
            pf = min(pf,100); pf = max(pf,0)
          }
          res = gauge(min(pf,100), 0,100,symbol = '%',label= "Progress...")
        }
      }
    }
    return(res)
  })
  
  
  
  
  output$build_process_table <- DT::renderDataTable({
    folder_name = paste('file',values$masque_key,sep="")
    build_process = paste(values$curdir,"www","masque","done",folder_name,"shaman_process_build.tsv",sep= .Platform$file.sep)
    if(file.exists(build_process))
    {
      bp = read.csv(build_process,sep="\t",header=TRUE)
    }
    return(bp)},
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE, processing=FALSE,server=TRUE))
  # ## plot gauge
  # output$gaugeMasque_progress <-renderGauge({
  # 
  #   res = NULL;
  #   num = round(as.numeric(values$num),1)
  # 
  #   CMP = isolate(CheckMasque(input, values))
  #   Error = CMP$Error
  #   if(is.null(Error) || num>1) res = gauge(min(num,100), 0,100,symbol = '%',label= "Progress...")
  # 
  #   return(res)
  # })
  # 
  # 
  textMASQUE <- reactive({
    
    samp = SamplesMasque(input,values)
    
    text = paste("<b>Type of data:</b>",input$DataTypeMasque,"<br /> <br /> ",
                 "<b>Paired-end sequencing:</b>",input$PairedOrNot)
    text = paste(text,"<br /> <br /> ","<b>Number of samples:</b>",length(samp$samples))
    text = paste(text,"<br /> <br /> ","<b>Removed samples:</b>",length(samp$samples_removed))
    
    if(isValidEmail(input$to)) text = paste(text,"<br /> <br /> ","<b>Email:</b>",input$to)
    
    if(input$HostName!="") text = paste(text,"<br /> <br /> ","<b>Host:</b>",input$HostName)
    if(input$primer && input$PairedOrNot=="n") {text = paste(text,"<br /> <br /> ","<b>Primer:</b>",input$primerSingle)}
    if(input$primer && input$PairedOrNot=="y") {text = paste(text,"<br /> <br /> ","<b>Primer forward:</b>",input$R1primer,
                                                             "<br /> <br /> ","<b>Primer reverse:</b>",input$R2primer)}
    
    return(text)
  })
  
  
  output$summary_box_masque <- renderUI({
    
    text = textMASQUE()
    
    res = div(id="boxsum",style = "word-wrap: break-word;",box(id="boxsum",
                                                               title = strong("Summary of your analysis"), width = NULL, background = "light-blue",
                                                               # HTML("<h4><b>Enter the key:</b></h4>",'<hr color="white"> <br /> <br /> '),
                                                               
                                                               HTML(text),
                                                               div(style = "text-align:right;",
                                                                   downloadButton("printMasque_summary", "Save"),
                                                                   tags$style(type='text/css', "#printMasque_summary {margin-top: 15px;}")
                                                               )
    )
    )
    return(res)
  })
  
  ## Export MASQUE summary in .txt
  output$printMasque_summary <- downloadHandler(
    filename = function() { 'Summary.txt' },
    content = function(file){
      txt = textMASQUE()
      txt = gsub("<br />",replacement = "\n",txt)
      txt = gsub("</b>",replacement = "",txt)
      txt = gsub("<b>",replacement = "",txt)
      write(paste("Summary of your analysis  \n \n ",txt), file)
    }
  )
  
  
  ## Send mail with the password
  observeEvent(input$checkMail,{
    
    # observe( info(paste("You will received a password by email at :",isolate(input$to), '\nThis can take few seconds.')))
    to <- isolate(input$to)
    subject <- "SHAMAN Analysis"
    body <- paste("Hello, \n You are using SHAMAN to run a quantitative metagenomic analysis. Hereafter is the key you need in SHAMAN :
                  \n",values$pass," \n \n Best regards, \n SHAMAN team")
    mailControl=list(smtpServer="smtp.pasteur.fr")
    from="<shaman@pasteur.fr>"
    ## Send mail
    if(Sys.info()["nodename"] == "ShinyPro"){
      sendmail(from=from,to=to,subject=subject,msg=body,control=mailControl)
    }
    
    ## Update the key value.
    updateTextInput(session,"password","",value = values$pass)
    
    ## Store the email 
    values$login_email = to
  })
  
  
  ## Create button once MASQUE computation is over
  # output$MasqueToShaman_button <- renderUI({
  #   res = NULL
  # 
  #   CMP = CheckMasque(input, values)
  #   Error = CMP$Error
  #   if(is.null(Error) && values$num>=100){
  # 
  # 
  #    }
  # 
  #   return(res)    
  # })
  #   
  
  # observe({
  #   if(values$num<100) disable("RunResMasque")
  #   if(values$num<100) disable("masque_database")
  #   
  #   if(values$num<100){ addPopover(session,"load-masque-res", 
  #                                 title= "Waiting for the results",
  #                                 content = paste("Once the computation is over, you will received a password by email at:",values$login_email)
  #                                 )
  #     } else removePopover(session, "load-masque-res")
  # })
  
  
  #observeEvent(input$RunResMasque,{
  #  sendSweetAlert(messageId="WTF2", title = "Success", text = "EUHHH", type = "success", html=TRUE)
  #updateSelectInput(session, "FileFormat","",selected = "fileBiom")
  #reset("fileBiom"); reset("fileTree"); print("WTG");
  #reset("fileBiom")
  #reset("fileTree")
  #values$biom_masque = paste(values$curdir,"www","masque","done",paste("file",values$masque_key,sep=""),paste("shaman_",input$masque_database,".biom",sep=""),sep= .Platform$file.sep)
  #values$tree_masque = paste(values$curdir,"www","masque","done",paste("file",values$masque_key,sep=""),paste("shaman_",input$masque_database,"_tree.nhx",sep=""),sep= .Platform$file.sep)
  #})
  
  observeEvent(input$LoadResMasque, {
    
    if (input$masque_db == "rdp"){
      updateSelectInput(session, "FileFormat","",selected = "fileCounts")
      updateSelectInput(session, "TypeTaxo","",selected = "RDP")
      reset("fileTaxo")
      reset("fileCounts")
      values$rdp_thres_masque = as.numeric(input$rdp_thres)
      values$rdp_annot_masque = paste(values$curdir,"www","masque","done",paste("file",values$masque_key,sep=""),"shaman_rdp_annotation.tsv",sep= .Platform$file.sep)
      values$count_table_masque = paste(values$curdir,"www","masque","done",paste("file",values$masque_key,sep=""),"shaman_otu_table.tsv",sep= .Platform$file.sep)
    }
    else{
      updateSelectInput(session, "FileFormat","",selected = "fileBiom")
      reset("fileBiom")
      reset("fileTree")
      values$biom_masque = paste(values$curdir,"www","masque","done",paste("file",values$masque_key,sep=""),paste("shaman_",input$masque_db,".biom",sep=""),sep= .Platform$file.sep)
      values$tree_masque = paste(values$curdir,"www","masque","done",paste("file",values$masque_key,sep=""),paste("shaman_",input$masque_db,"_tree.nhx",sep=""),sep= .Platform$file.sep)
    }
    sendSweetAlert(messageId="LoadResMasque", title = "Success", text = "Processed data were successfully loaded. You can go to statistical analysis part to analyze your data with SHAMAN.", type = "success", html=TRUE)
  })
  
  
  
  #from home
  observeEvent(input$LoadResMasque_home, {
    
    if (input$masque_db_home == "rdp"){
      updateSelectInput(session, "FileFormat","",selected = "fileCounts")
      updateSelectInput(session, "TypeTaxo","",selected = "RDP")
      reset("fileTaxo")
      reset("fileCounts")
      values$rdp_thres_masque = as.numeric(input$rdp_thres_home)
      values$rdp_annot_masque = paste(values$curdir,"www","masque","done",paste("file",values$masque_key,sep=""),"shaman_rdp_annotation.tsv",sep= .Platform$file.sep)
      values$count_table_masque = paste(values$curdir,"www","masque","done",paste("file",values$masque_key,sep=""),"shaman_otu_table.tsv",sep= .Platform$file.sep)
    }
    else{
      updateSelectInput(session, "FileFormat","",selected = "fileBiom")
      reset("fileBiom")
      reset("fileTree")
      values$biom_masque = paste(values$curdir,"www","masque","done",paste("file",values$masque_key,sep=""),paste("shaman_",input$masque_db_home,".biom",sep=""),sep= .Platform$file.sep)
      values$tree_masque = paste(values$curdir,"www","masque","done",paste("file",values$masque_key,sep=""),paste("shaman_",input$masque_db_home,"_tree.nhx",sep=""),sep= .Platform$file.sep)
    }
    sendSweetAlert(messageId="LoadResMasque_home", title = "Success", text = "Processed data were successfully loaded. You can go to statistical analysis part to analyze your data with SHAMAN.", type = "success", html=TRUE)
  })
  
  ## Export results in .zip 
  output$Download_masque_zip <- downloadHandler(
    filename = function() { paste("SHAMAN_",values$masque_key,'.zip',sep="")},
    content = function(file){
      zip_file = paste(values$curdir,"www","masque","done",paste("shaman_",values$masque_key,".zip",sep=""),sep= .Platform$file.sep)
      if(file.exists(zip_file)){file.copy(zip_file, file)}
    }
  )
  
  ## Export results in .zip from home
  output$Download_masque_zip_home <- downloadHandler(
    filename = function() { paste("SHAMAN_",values$masque_key,'.zip',sep="")},
    content = function(file){
      zip_file = paste(values$curdir,"www","masque","done",paste("shaman_",values$masque_key,".zip",sep=""),sep= .Platform$file.sep)
      if(file.exists(zip_file)){file.copy(zip_file, file)}
    }
  )
  
  
  
  output$Project_box_home <- renderUI({
    
    res = NULL
    PS = Project_status(values$masque_key,values$curdir)
    
    if(PS$passOK){
      
      if(PS$status=="done")
      {
        res = list()
        json_file = PS$file
        folder_name = paste('file',values$masque_key,sep="")
        json_file = paste(curdir,"www","masque","done",paste(folder_name, ".json", sep=""),sep= .Platform$file.sep)
        json_data = rjson::fromJSON(file=json_file)
        
        ### Paste file name as folder
        annot_process = paste(curdir,"www","masque","done",folder_name,"shaman_process_annotation.tsv",sep= .Platform$file.sep)
        
        ## Waiting for file creation (max 3min)
        start = Sys.time(); diff = 0
        while(!file.exists(annot_process) && diff<180){
          annot_process = paste(curdir,"www","masque","done",folder_name,"shaman_process_annotation.tsv",sep= .Platform$file.sep)
          tmp = Sys.time()
          diff = tmp-start
        }
        
        if(file.exists(annot_process))
        {
          
          ap = read.csv(annot_process,sep="\t")
          if(json_data$type == "16S") db_choices = c("Silva" = "silva","Greengenes" = "greengenes", "RDP"= "rdp")
          else if(json_data$type == "23S_28S" || json_data$type == "18S") db_choices = c("Silva" = "silva","RDP"= "rdp")
          else db_choices = c("Findley" = "findley", "Underhill"= "underhill", "Unite"= "unite", "RDP"= "rdp")
          
          res[[1]] = fluidRow(
            box(title="Load the results",width = 3, background = "light-blue",
                selectInput("masque_db_home","Select the database",choices=db_choices),
                conditionalPanel(condition="input.masque_db_home=='rdp'",numericInput("rdp_thres_home",h6(strong("Threshold:")),0.5,step=0.01,min=0.01,max=1)),
                actionButton("LoadResMasque_home", "Upload the results",icon=icon('upload')),
                tags$style(type='text/css', "#LoadResMasque_home { width:100%; }"),
                receiveSweetAlert(messageId = "LoadResMasque_home")
            ),
            box(title="Download .zip file",width = 3, status = "success",
                downloadButton('Download_masque_zip_home', 'Download the results'),
                tags$style(type='text/css', "#Download_masque_zip_home { width:100%;}")
            )
          )
        }
      }
    }
    
    return(res)
  })
  
  ######################## END MASQUE #################################
  
  
  
  
  observeEvent(input$deleteRows,{
    
    if (!is.null(input$DataTarget_rows_selected)) {
      
      if(nrow(values$TargetWorking)!=0) values$labeled <- (nrow(values$TargetWorking)-length(input$DataTarget_rows_selected))*values$labeled/nrow(values$TargetWorking)
      else values$labeled <- 0
      values$TargetWorking <- values$TargetWorking[-as.numeric(input$DataTarget_rows_selected),]
      
    }
  })
  
  
  
  
  
  ## Interest Variables
  output$SelectInterestVar <- renderUI({
    
    target=values$TargetWorking
    if(!is.null(target)) 
    {
      namesTarget = colnames(target)[2:ncol(target)]
      selectInput("InterestVar",h6(strong("Select the variables")),namesTarget,selected=namesTarget,multiple=TRUE)
    }
    
  })
  
  ## Interactions
  output$SelectInteraction2 <- renderUI({
    
    target = values$TargetWorking
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
    
    target=values$TargetWorking
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
    
    target = values$TargetWorking
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
    labeled = values$labeled
    if(!is.null(labeled)) 
    {
      #### Ajout fontion check target
      #infoBox(h6(strong("Target file")), subtitle = h6("Your target file is OK"), icon = icon("thumbs-o-up"),color = "green",width=NULL,fill=TRUE)
      labeled = round(labeled,2)
      if(labeled>0) res = shinydashboard::valueBox(paste0(labeled, "%"),h6(strong("Labeled features")), color = "green",width=NULL,icon = icon("list"))
      else res = shinydashboard::valueBox(paste0(labeled, "%"),h6(strong("Labeled features")), color = "red",width=NULL,icon = icon("list"))  
    }
    #else infoBox(h6(strong("Target file")), subtitle = h6("Label of the target file must correspond to count table column names") ,color = "light-blue",width=NULL,fill=TRUE, icon = icon("warning"))
    else res = shinydashboard::valueBox(paste0(0, "%"),h6(strong("Labeled features")), color = "light-blue",width=NULL,icon = icon("list"))
    return(res)
  }
  )
  
  
  
  
  ## target table
  output$DataTarget <- DT::renderDataTable(
    values$TargetWorking,
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE, processing=FALSE,server=TRUE
    ))
  
  Target_selection <- reactive ({
    input$DataTarget_rows_selected
  })
  
  
  
  ## Counts table for the selected taxonomy level
  output$CountsMerge <- DT::renderDataTable(
    round(counts(ResDiffAnal()$dds,normalized=TRUE)),
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                   pageLength = 10,scrollX=TRUE, processing=FALSE
    ))
  
  
  ## Box for merged counts
  output$BoxCountsMerge <- renderUI({
    input$RunDESeq
    #counts = isolate(dataMergeCounts()$counts)
    counts = isolate(dataMergeCounts()$CT_Norm)
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
    #content = function(file){write.csv(dataMergeCounts()$counts, file)}
    content = function(file){write.csv(round(counts(ResDiffAnal()$dds,normalized=TRUE)), file)}
  )
  
  ## Export in .csv
  output$ExportRelative <- downloadHandler(
    filename = function() { 'RelativeAb.csv' },
    content = function(file){write.csv(sweep(round(counts(ResDiffAnal()$dds,normalized=TRUE)),2,colSums(round(counts(ResDiffAnal()$dds,normalized=TRUE))),`/`), file)}
  )
  
  ## Export size factors
  output$ExportSizeFactor <- downloadHandler(
    filename = function() { if (input$sepsizef == "\t") 'SHAMAN_sizefactors.tsv' else 'SHAMAN_sizefactors.csv' },
    content = function(file){write.table(SizeFactor_table(), file,quote=FALSE,row.names = FALSE,sep=input$sepsizef)}
  )
  
  ## Export in .csv
  output$ExportTarget <- downloadHandler(
    filename = function() { 'Target.csv' },
    content = function(file){write.table(values$TargetWorking, file,sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE)}
  )
  
  ## Box for target visualisation
  output$BoxTarget <- renderUI({
    
    target = values$TargetWorking
    
    if(!is.null(target) &&  nrow(target)>0)
    {
      box(title="Target file overview",width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = TRUE,
          DT::dataTableOutput("DataTarget"),
          actionButton("deleteRows", "Delete samples"),
          downloadButton('ExportTarget', 'Export target file')
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
  
  ## Add Phylogenetic tree in the result
  observeEvent(input$fileTree,{
    updateSelectizeInput(session, "PlotVisuSelect", "", 
                         c("Barplot"="Barplot","Heatmap"="Heatmap","Boxplot"="Boxplot",
                           "Tree"="Tree","Scatterplot"="Scatterplot","Diversity"="Diversity",
                           "Rarefaction"="Rarefaction","Krona"="Krona","Phylogeny"="Phylogeny"))
  }, priority=1)
  
  
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
    
    target = values$TargetWorking
    design = GetDesign(input,target)
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
    create_forked_task(ResDiffAnal())
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
  #     target = values$TargetWorking
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
    target = values$TargetWorking
    labeled = values$labeled
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
    target = values$TargetWorking
    labeled = values$labeled
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
        toastr_error(title="Error",message=paste("<h5>","This file can not be loaded.","</h5>","<br/> <br/>","<em>The loaded file is not in the biom format or its format is not currently supported by SHAMAN software.</em>"),
                     closeButton = TRUE,position ="bottom-right",preventDuplicates = TRUE,newestOnTop = TRUE,
                     progressBar = TRUE,showDuration = 300,showMethod="show",timeOut = 20000,extendedTimeOut = 2000)
      }
    }
  })
  
  
  output$InfoCountsFile<- renderUI({
    
    if(input$FileFormat=="fileCounts")
    {
      inFile <- input$fileCounts
      Counts = dataInputCounts()
      
      if(!is.null(inFile) && is.null(Counts)) {   
        toastr_error(title="Error",message=paste("<h5>","This file can not be loaded.","</h5>","<br/> <br/>","<em>The count table file is not in the correct format for SHAMAN software.</em>"),
                     closeButton = TRUE,position ="bottom-right",preventDuplicates = TRUE,newestOnTop = TRUE,
                     progressBar = TRUE,showDuration = 300,showMethod="show",timeOut = 20000,extendedTimeOut = 2000)
      }
    }
  })
  
  
  output$InfoTaxoFile<- renderUI({
    
    if(input$FileFormat=="fileCounts")
    {
      inFile <- input$fileTaxo
      Taxo = dataInputTaxo()
      
      if(!is.null(inFile) && !input$NoTaxoFile && is.null(Taxo)) {   
        toastr_error(title="Error",message=paste("<h5>","This file can not be loaded.","</h5>","<br/> <br/>","<em>The taxonomy table file is not in the correct format for SHAMAN software.</em>"),
                     closeButton = TRUE,position ="bottom-right",preventDuplicates = TRUE,newestOnTop = TRUE,
                     progressBar = TRUE,showDuration = 300,showMethod="show",timeOut = 20000,extendedTimeOut = 2000)
      }
    }
  })
  
  
  
  
  
  
  
  #####################################################
  ##
  ##                Diagnostic plots
  ##
  #####################################################
  
  
  
  output$VarIntDiag <- renderUI({
    
    target=values$TargetWorking
    
    if(!is.null(target)) 
    {
      namesTarget = colnames(target)[2:ncol(target)]
      selectizeInput("VarInt",h6(strong("Select the variables of interest (max 2)")),namesTarget, selected = namesTarget[1],multiple = TRUE,options = list(maxItems = 2))
    }
    
  })
  
  
  output$PlotDiag <- renderPlot({
    input$RunDESeq
    
    resDiff = isolate(ResDiffAnal())
    ## Phylogenetic tree
    tree = dataInputTree()$data
    
    Plot_diag(input,resDiff,tree)
  },height = reactive(input$heightDiag), width = reactive(ifelse(input$modifwidthDiag,input$widthDiag,"auto")))
  
  
  
  ## Select PCA/PCOA axis
  output$PC1_sel <-renderUI ({
    res = NULL
    resDiff = ResDiffAnal()
    ## Phylogenetic tree
    tree = dataInputTree()$data
    
    if(!is.null(resDiff)){ 
      pca_tab = Plot_diag(input,resDiff,tree,getTable=TRUE)
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
    ## Phylogenetic tree
    tree = dataInputTree()$data
    
    if(!is.null(resDiff)){ 
      pca_tab = Plot_diag(input,resDiff,tree,getTable=TRUE)
      if(!is.null(pca_tab)) 
      {
        pc_axes = paste("PC",seq(1,ncol(pca_tab)),sep="")
        res = selectizeInput("PCaxe2","Y-axis",pc_axes,selected=pc_axes[min(2,ncol(pca_tab))])
      }
    }
    return(res)
  })
  
  output$PlotnmdsStress <- renderPlot({
    
    resDiff = ResDiffAnal()
    ## Phylogenetic tree
    tree = dataInputTree()$data
    
    Plot_diag_nmdsStress(input,resDiff,tree)
  },height = 400)
  
  
  output$PlotpcoaEigen <- renderPlot({
    
    resDiff = ResDiffAnal()
    ## Phylogenetic tree
    tree = dataInputTree()$data
    
    Plot_diag_pcoaEigen(input,resDiff,tree)
  },height = 400)
  
  
  output$PlotEigen <- renderPlot({
    
    resDiff = ResDiffAnal()
    Plot_diag_Eigen(input,resDiff)
  },height =400)
  
  
  SizeFactor_table <-reactive({ 
    res = ResDiffAnal()
    return(t(data.frame(Factor=res$normFactors)))
    
  })
  
  
  output$ResPermaTestBox <- renderUI({
    
    resDiff = ResDiffAnal()
    ## Phylogenetic tree
    tree = dataInputTree()$data
    
    resTest = Perma_test_Diag(input,resDiff,tree)
    resBox = NULL
    if(!is.null(resDiff) && !is.null(resTest))
    {    
      res = list()
      ## Title
      res$title = paste("<center><b><font size='+3'>","Permanova test","</font></b></center><br/>")
      ## Subtitle
      res$subtitle = paste("<center><em>","Analysis of variance using distance matrices","</em></center><br/>")
      ## Pvalue   
      res$ccl = paste("<center><b><font size='+1'>p-value :",round(resTest,5),"</font></b></center><br/>")
      
      resBox = box(title="Permanova ",width = 6, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = FALSE,
                   HTML(unlist(res))
      )
    } 
    return(resBox)
  })
  
  
  
  output$DistList <-renderUI ({
    tree = dataInputTree()$data
    ErrorTree = dataInputTree()$Error
    TaxoSelect = input$TaxoSelect
    
    dist_phyl = getDistMethods()[!getDistMethods() %in% c("additive_symm", "jensen-shannon", "jensen_difference", "minkowski", "topsoe")]
    
    
    res = selectInput("DistClust","Distance", unique(sort(c("altGower", "binomial", "bray", "canberra", "cao", "chao", "euclidean","gower", "horn",
                                                            "jaccard", "kulczynski",  "mahalanobis", "morisita", "mountford","raup",
                                                            "SERE"="sere", dist_phyl))),selected="bray")
    
    ## Add the unifrac distance
    if(!is.null(tree) && !is.null(input$fileTree) && is.null(ErrorTree)  && TaxoSelect %in% c("OTU/Gene", "MGS"))
    {
      res = selectInput("DistClust","Distance",unique(sort(c("altGower", "binomial", "bray", "canberra", "cao", "chao", "euclidean","gower", "horn",
                                                             "jaccard","kulczynski",  "mahalanobis", "morisita", "mountford", "raup",
                                                             "SERE"="sere","Unifrac", dist_phyl))),selected="bray")
    }
    
    return(res)
  })
  
  
  output$SizeFactTable <- DT::renderDataTable(
    SizeFactor_table(),
    options = list(scrollX=TRUE,searching = FALSE, processing=FALSE
    ))
  
  
  ## Select Modality DiagPlot
  
  output$ModMat <- renderUI({
    
    VarInt = input$VarInt
    target = values$TargetWorking
    
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
    target = values$TargetWorking
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
      resDiff = ResDiffAnal()
      tree = isolate(dataInputTree()$data)
      
      
      print(Plot_diag(input,resDiff,tree))
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
    target = values$TargetWorking
    labeled = values$labeled
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
  
  
  output$PhyloTreeMetaR2 <- renderPhyloTreeMetaR({
    resDiff = ResDiffAnal()
    taxo_table = dataInput()$data$taxo
    treeseq = dataInputTree()$treeseq
    
    
    if(!is.null(resDiff$dds) && length(input$VisuVarInt)>=1 ) 
      if(input$NormOrRaw=="norm") withProgress(message="Loading...", Plot_Visu_Phylotree(input, resDiff, dataMergeCounts()$CT_Norm, taxo_table, treeseq))
    else withProgress(message="Loading...", Plot_Visu_Phylotree(input, resDiff, dataMergeCounts()$CT_noNorm, taxo_table, treeseq))
  })
  
  
  output$PlotVisuTree <- renderTreeWeightD3({
    resDiff = ResDiffAnal()
    taxo_table = dataInput()$data$taxo
    
    
    res = NULL
    if(!is.null(resDiff$dds) && length(input$VisuVarInt)>=1){
      if(input$NormOrRaw=="norm") res = Plot_Visu_Tree(input,resDiff,dataMergeCounts()$CT_Norm,taxo_table)
      else res = Plot_Visu_Tree(input,resDiff,dataMergeCounts()$CT_noNorm,taxo_table)
    } 
    return(res)
    
  })
  
  KronaR =function(){
    resDiff = ResDiffAnal()
    taxo_table = dataInput()$data$taxo
    
    res = NULL
    if(!is.null(resDiff$dds) && length(input$VisuVarInt)>=1){
      if(input$NormOrRaw=="norm") res = Plot_Visu_Krona(input,resDiff,dataMergeCounts()$CT_Norm,taxo_table)
      else res = Plot_Visu_Krona(input,resDiff,dataMergeCounts()$CT_noNorm,taxo_table)
    }
    tempdir = tempdir()
    temp = tempfile(pattern = "file", tmpdir=tempdir, fileext =".tsv")
    write.table(res, file=temp, quote=F, sep="\t", row.names =F, col.names=F)
    Sys.chmod(tempdir, mode = "0777")
    Sys.chmod(temp, mode = "0777")
    
    return(temp)
  }
  
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
      if(input$sepdiversity=="\t") 'SHAMAN_Diversity.tsv'
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
    else if(input$PlotVisuSelect=="Heatmap") res =  d3heatmapOutput("heatmap", height = input$heightVisu+10, width=ifelse(input$modifwidthVisu,input$widthVisu,"100%"))
    else if(input$PlotVisuSelect=="Boxplot") res = plotOutput("Boxplot", height = input$heightVisu+10, width=if(input$modifwidthVisu){input$widthVisu})
    else if(input$PlotVisuSelect=="Krona"){
      if(Sys.info()["nodename"] == "ShinyPro"){
        res= tags$iframe(src=paste0("http://hub05.hosting.pasteur.fr/aghozlane/KronaRShy/?parameter=",KronaR()), height = input$heightVisu+10, width=ifelse(input$modifwidthVisu,input$widthVisu,"100%"), seamless=NA)
      }
      else{
        res= tags$iframe(src=paste0("http://127.0.0.1:5438/?parameter=",KronaR()), height = input$heightVisu+10, width=ifelse(input$modifwidthVisu,input$widthVisu,"100%"), seamless=NA)  
      }
    }
    else if(input$PlotVisuSelect=="Tree") res = treeWeightD3Output('PlotVisuTree', height = input$heightVisu+10,width=ifelse(input$modifwidthVisu,input$widthVisu,"100%"))
    else if(input$PlotVisuSelect=="Scatterplot" && !input$AddRegScatter) res = scatterD3Output("ScatterplotD3", height = input$heightVisu+10, width=ifelse(input$modifwidthVisu,input$widthVisu,"100%"))
    else if(input$PlotVisuSelect=="Scatterplot" && input$AddRegScatter) res = plotOutput("Scatterplotgg", height = input$heightVisu+10,width=if(input$modifwidthVisu){input$widthVisu})
    else if(input$PlotVisuSelect=="Diversity") res =  plotOutput("DiversityPlot", height = input$heightVisu+10, width=if(input$modifwidthVisu){input$widthVisu})
    else if(input$PlotVisuSelect=="Rarefaction") res = plotOutput("RarefactionPlot",dblclick = "RarefactionPlot_dblclick",brush = brushOpts(id = "RarefactionPlot_brush",resetOnNew = TRUE), height = input$heightVisu+10, width=if(input$modifwidthVisu){input$widthVisu})
    else if(input$PlotVisuSelect=="Phylogeny") res = PhyloTreeMetaROutput('PhyloTreeMetaR2')
    #print(res)
    return(res)
  })
  
  ## Comparison plots
  output$plotVisuComp <- renderUI({
    
    res=NULL
    if(input$PlotVisuSelectComp=="Heatmap_comp") res =  d3heatmapOutput("heatmap_comp", height = input$heightVisuComp+10, width=ifelse(input$modifwidthComp,input$widthComp,"100%"))
    if(input$PlotVisuSelectComp=="Venn") res =  d3vennROutput("VennD3", height = input$heightVisuComp+10, width=ifelse(input$modifwidthComp,input$widthComp,"100%"))
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
    
    target=values$TargetWorking
    
    if(!is.null(target)) 
    {
      namesTarget = colnames(target)[2:ncol(target)]
      selectizeInput("VisuVarInt",h6(strong("Select the variables of interest")),namesTarget, selected = namesTarget[1],multiple = TRUE)
    }
    
  })
  
  
  output$VarIntVisuScatter <- renderUI({
    
    target=values$TargetWorking
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
  
  
  # output$VarIntVisuTree <- renderUI({
  # 
  #   target=values$TargetWorking
  #   data = dataInput()$data
  #   taxo = input$TaxoSelect
  #   resDiff = ResDiffAnal()
  #   res = NULL
  # 
  #   if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0 && !is.null(taxo) && taxo!="..." && !is.null(target))
  #   {
  #     counts = dataMergeCounts()$counts
  # 
  #     Available_x = sort(rownames(counts))
  # 
  #     res = selectizeInput("TaxoTree",h6(strong(paste("Select a specific",taxo,sep=" "))),c("...",Available_x),multiple = TRUE)
  # 
  #   }
  # 
  #   return(res)
  # 
  # })
  
  #####################################################
  ##
  ##                KRONA
  ##
  #####################################################
  #output$kronar <- renderTable({
  #  data = dataInput()$data 
  #  taxo = input$TaxoSelect
  #  if(!is.null(data$counts) && !is.null(data$taxo) && nrow(data$counts)>0 && nrow(data$taxo)>0 && !is.null(taxo) && taxo!="...") 
  #  {
  #    print(counts)
  #    print(data)
  #KronaR(dat) 
  #print(data$counts)
  #krona_table=tempfile(pattern = "krona", tmpdir = tempdir(), fileext = "")
  #url=paste(krona_table, ".html", sep="")
  #system(paste("export PERL5LIB=/home/aghozlan/workspace/SHAMAN_App/KronaTools-2.6/lib:$PERL5LIB; /home/aghozlan/workspace/META10S_App/krona_bin/bin/ktImportText", krona_table))
  #system(paste("ktImportText", krona_table))
  #refs <- paste0("<a href='",  url, "' target='_blank'>krona</a>")
  
  #data.frame(refs)
  #  }
  #}, sanitize.text.function = function(x) x)
  
  
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
    target = values$TargetWorking
    labeled= values$labeled
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
    target=values$TargetWorking
    
    if(is.null(target) || taxo =="...") 
    {
      shinyjs::disable("AddFilter")
    } else {
      shinyjs::enable("AddFilter")
    }
  })
  
  
  
  
})


