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
##      LOAD FILES
##
#####################################################

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
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#   ResOutput <-eventReactive(input$RunAnalyse,{
#     
#     data = dataInput() 
# 
#   })
  
  
  
  output$SelectDataVar <- renderUI({
    
  data=dataInputTmp()
    
    selectInput("VariableSelect",p(strong("Choisissez les variables à étudier"),h6(em("Sélection multiple avec CTRL"))),
                colnames(data),
                selected=colnames(data),multiple=TRUE,size=2,selectize=FALSE)
  })
  
  
  
  output$SelectVarQuanti <- renderUI({
    
    data=dataInputTmp()
    
    Num = sapply(data,is.numeric)
    namesQuanti = colnames(data)[Num]
    
    selectInput("VariableSelectQuanti",h5(strong("Sélectionnez les variables numériques que vous souhaitez transformer en variables qualitatives"),h6(em("Sélection multiple avec CTRL"))),
                namesQuanti,
                selected=1,multiple=TRUE,size=1,selectize=FALSE)
  })
  

    output$QuantiToQuali <- renderUI({
    
      input$GoQuali
      names = isolate(input$VariableSelectQuanti)
      namesVar = paste(names,collapse = ", ")
      if(length(names)>0) HTML(paste("Les variables suivantes seront considérées comme des",'<b>',"variables qualitatives :",'</b>','<br/>',namesVar))
      else return(NULL)
       })
  
  
  dataInputType <- reactive({
    
    input$GoQuali
    
    data=dataInputTmp()
    names= colnames(data)
    
    if(!is.null(data)) 
    {
      ind=which(colnames(data)%in%isolate(input$VariableSelectQuanti))
      data[,ind]=as.data.frame(sapply(data[,ind],as.factor))      
    }
    
    return(data)
    
  })
  
  
  dataInput <- reactive({
    
    input$RefreshData
  
    data=dataInputType()
    names= colnames(data)
    rownamesTmp = rownames(data)
    
    if(!is.null(data) && !is.null(isolate(input$VariableSelect))) 
    {
      ind=which(colnames(data)%in%isolate(input$VariableSelect))
      data=as.data.frame(data[,ind])
      
      # Get the names of rows and columns
      colnames(data) = names[ind]
      rownames(data) = rownamesTmp
    }

    return(as.data.frame(data))
    
  })
  
#   observeEvent(input$RunAnalyse,{  
#     
#     ResOutput()
#   })
  
  output$DataBrutes <- renderDataTable(
    dataInput(), 
    options = list(lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
    pageLength = 10,scrollX=TRUE
  ))
  
  
  
  output$SelectUniVar <- renderUI({
    input$RunAnalyse
    
    data = dataInput()
    names = colnames(data)
 
    if(!is.null(data)) selectInput("UniVar","Variable à étudier",names)
  })
  
  
  ###########################################
  ##
  ##    Accueil exemples
  ##
  ###########################################
  

  
  ###########################################
  ##
  ##    TABLE DESCRIPTIVES STATS
  ##
  ###########################################
  
  TableQuantiOut <- reactive({
    
    input$RefreshStat
    
    data = dataInput()
    
    Num = sapply(data,is.numeric)
    dataNum = as.data.frame(data[,Num])
    
    stat = as.data.frame(signif(apply(dataNum,2,DescriptiveStat),2))
    stat = stat[as.numeric(isolate(input$IndicQuanti)),]

    return(stat)
  })
  
  
  TableQualiOut <- reactive({
    
    input$RefreshStatQuali
    
    data = dataInput()
    Num = sapply(data,is.numeric)
    dataQuali = data[,!Num]
    namesQuali = names(data)[!Num]
    
    DescriptiveStatQuali(as.data.frame(dataQuali),namesQuali,isolate(as.numeric(input$IndicQuali)))
    #as.data.frame(tmp[,isolate(as.numeric(input$IndicQuali))])
  })
  

  
  CorTab <- reactive({
    
    data = dataInput()
    Num = sapply(data,is.numeric)
    dataQuanti = as.data.frame(data[,Num])
    corel = corr.test(dataQuanti,method=input$CorelMeth)
    
    return(list(cor=signif(corel$r,2),pval=signif(corel$p,2)))
  })

    
  output$TableQuanti <- renderDataTable(
    TableQuantiOut(), 
    options = list(paging=FALSE, searching=FALSE,
                   pageLength = 12,scrollX=TRUE
    ))
  
  
  output$TableQuali <- renderDataTable(
    TableQualiOut(), 
    options = list(paging=FALSE, searching=FALSE,
                   pageLength = 10,scrollX=TRUE
    ))
  
  
  output$CorTable <- renderDataTable(
    
    CorTab()$cor,
    options = list(paging=FALSE, searching=FALSE,
                   pageLength = 10,scrollX=TRUE
    ))
  
  output$CorTableTest <- renderDataTable(
    
    CorTab()$pval,
    options = list(paging=FALSE, searching=FALSE,
                   pageLength = 10,scrollX=TRUE
    ))
  
  
  
  ###########################################
  ##
  ##   Analyse univarié
  ##
  ###########################################
  
  
  TableOutUni <- reactive({
    
    data = dataInput()
    ind = which(colnames(data)%in%input$UniVar)
    
    dataTmp = data[,ind]
    namesTmp = names(data)[ind]
    
    if(is.numeric(dataTmp))
    {
      stat = as.data.frame(signif(DescriptiveStat(dataTmp),2))
      names(stat)= namesTmp
    }
    if(!is.numeric(dataTmp))
    {
      stat = DescriptiveStatQuali(as.data.frame(dataTmp),namesTmp,seq(1,3))
    }
    stat=as.data.frame(stat)
    return(summary(dataTmp))
  })
  
  
#   output$TableUni <- renderPrint({
#     data = dataInput()
#     ind = which(colnames(data)%in%input$UniVar)
#     
#     dataTmp = data[,ind]
#     namesTmp = names(data)[ind]
#     
#     return(summary(dataTmp))
#     })
  
  
  #######################################
  ##
  ##    Info Box
  ##
  ########################################
  
  output$NumberRowBox <- renderInfoBox({
    data = dataInput()
    nbrow = nrow(data)
    infoBox(h5("Lignes"), nbrow, icon = icon("arrows-v"),color = "light-blue",fill=TRUE)
  })
  
  output$NumberColBox <- renderInfoBox({
    data = dataInput()
    nbcol = ncol(data)
    infoBox(h5("Colonnes"), nbcol, icon = icon("arrows-h"),color = "light-blue",fill=TRUE)
  })
  
  output$NumberQuantiBox <- renderInfoBox({
    data = dataInput()
    Num = sapply(data,is.numeric)
    infoBox(h5("Variables numériques"), length(which(Num)), icon = icon("line-chart"),color = "light-blue",fill=TRUE)
  })
  
  
  
  output$NumberQualiBox <- renderInfoBox({
    data = dataInput()
    Num = sapply(data,is.numeric)
    infoBox(h5("Variables qualitatives"), length(which(!Num)), icon = icon("pie-chart"),color = "light-blue",fill=TRUE)
  })
  
  
  output$RunOK <- renderInfoBox({
    input$RunAnalyse
    data = isolate(dataInput())
    if(nrow(data)>0)
    { 
      infoBox(h5(""), "L'analyse statistique a été réalisée !", icon = icon("unlock"),color = "light-blue",fill=TRUE)
    }
    else infoBox(h5(""), "Veuillez charger vos données et appuyer sur Lancer l'analyse", icon = icon("lock"),color = "light-blue",fill=TRUE)
  })
  
  output$ResTTestBox <- renderInfoBox({
    
    input$ExecuteTtestCible
    
    data = dataInput()
    ind = which(colnames(data)%in%input$UniVar)
    
    dataTmp = data[,ind]
    TestNum = is.numeric(dataTmp)
    dataTmp = data.frame(x=dataTmp)
    namesTmp = names(data)[ind]
    
    ## Numeric data
    if(TestNum)
    {
      test = t.test(dataTmp$x,mu=as.numeric((isolate(input$ValCibleTtest))))
      pval = signif(test$p.value,3)
      if(pval<=as.numeric(isolate(input$alphaTtest))/100) inf = infoBox(h4("Significatif"), paste("p-value :",pval), icon = icon("thumbs-up"),color = "green",fill=TRUE)
      if(pval>as.numeric(isolate(input$alphaTtest))/100) inf = infoBox(h4("Non significatif"), paste("p-value :",pval), icon = icon("thumbs-down"),color = "red",fill=TRUE)
    }

        if(!TestNum){
      tmp = table(dataTmp$x)
      LenTMP=length(tmp)
      test = chisq.test(tmp,p=rep(1/LenTMP,LenTMP))
      pval = signif(test$p.value,3)
      if(pval<=as.numeric(isolate(input$alphaTtest))/100) inf = infoBox(h4("Significatif"), paste("p-value :",pval), icon = icon("thumbs-up"),color = "green",fill=TRUE)
      if(pval>as.numeric(isolate(input$alphaTtest))/100) inf = infoBox(h4("Non significatif"), paste("p-value :",pval), icon = icon("thumbs-down"),color = "red",fill=TRUE)
    }
  return(inf)
    
  })
  
  
  
  output$ResTTest2sampBox <- renderInfoBox({
    
    input$ExecuteTtestMoy
    
    data = dataInput()
    
    var1 = input$VariableSelectBi1; var2 = input$VariableSelectBi2
    ind1  = which(colnames(data)%in%c(var1)); ind2  = which(colnames(data)%in%c(var2))
    
    ## Get the indexes
    ind =  c(ind1,ind2)
    
    ## Which var is numeric ?
    Num = sapply(data[,unique(ind)],is.numeric)
    
    if(length(unique(ind))==2)
    {
      if(length(which(Num))==2) 
        {
          data2 = data.frame(x=data[,ind1],y=data[,ind2])
          
          
        }

      if(length(which(Num))==1) data2 = data.frame(x=data[,ind[which(Num)]],y=data[,ind[which(!Num)]])        
    }
    
    if(length(unique(ind))==1)  data2 = data.frame(x=data[,ind],y=data[,ind])
    
    
    ## Numeric data
    if(TestNum)
    {
      test = t.test(dataTmp$x,mu=as.numeric((isolate(input$ValCibleTtest))))
      pval = signif(test$p.value,3)
      if(pval<=as.numeric(isolate(input$alphaTtest))/100) inf = infoBox(h4("Significatif"), paste("p-value :",pval), icon = icon("thumbs-up"),color = "green",fill=TRUE)
      if(pval>as.numeric(isolate(input$alphaTtest))/100) inf = infoBox(h4("Non significatif"), paste("p-value :",pval), icon = icon("thumbs-down"),color = "red",fill=TRUE)
    }
    
    if(!TestNum){
      tmp = table(dataTmp$x)
      LenTMP=length(tmp)
      test = chisq.test(tmp,p=rep(1/LenTMP,LenTMP))
      pval = signif(test$p.value,3)
      if(pval<=as.numeric(isolate(input$alphaTtest))/100) inf = infoBox(h4("Significatif"), paste("p-value :",pval), icon = icon("thumbs-up"),color = "green",fill=TRUE)
      if(pval>as.numeric(isolate(input$alphaTtest))/100) inf = infoBox(h4("Non significatif"), paste("p-value :",pval), icon = icon("thumbs-down"),color = "red",fill=TRUE)
    }
    return(inf)
    
  })
  
  
  
  output$plotuni<- renderPlot({
    
    data = dataInput()
    generateUniPlot(input,data)
  })
  
  
  output$RadioUniPlot <- renderUI({
    
    data = dataInput()
    ind = which(colnames(data)%in%input$UniVar)

    if(is.numeric(data[,ind]))
    {
      radioButtons("RadioPlotUni","Choisir une représentation",
                   choices=c("Histogramme"="hist","Boxplot"="box","Densités"="densities","Q-Q plot"="qqplot"))
    }
    else
    {
      radioButtons("RadioPlotUni","Choisir une représentation",
                   choices=c("Diagramme en barre"="BarPlot","Pie chart"="Pie"))
    }
    
    
  })
  
  
  
  
  output$CheckTestBi <- renderUI({
    
    data = dataInput()
    
    var1 = input$VariableSelectBi1; var2 = input$VariableSelectBi2
    ind1  = which(colnames(data)%in%c(var1)); ind2  = which(colnames(data)%in%c(var2))
    
    ## Get the indexes
    ind =  c(ind1,ind2)
    
    ## Which var is numeric ?
    Num = sapply(data,is.numeric)
    Num1 = Num[ind1]; Num2 = Num[ind2]
    cbox=NULL
    
    if(length(unique(ind))==2)
    {
      if(Num1 && Num2) 
      {
        cbox = checkboxGroupInput("TestBi", "Tests statistiques :",
                           c("Test d'égalité des moyennes" = "testmoy",
                             "Test corrélation" = "testcor"))
      }
      
      if((Num1 && !Num2) || (!Num1 && Num2))
      {
        cbox = checkboxGroupInput("TestBi", "Test statistique :",
                           c("ANOVA" = "TestAnova"))
      }
      
      if(!Num1 && !Num2) 
      {
        cbox =checkboxGroupInput("TestBi", "Test statistique :",
                      c("Test du Chi-2" = "TestChi2"))
      }
    }
    
    if(length(unique(ind))==1)
      {
        textInput("OneVarSelected","")
      }      
    return(cbox)
  })
  

  
  
  
  ####################################
  ##
  ##      ANALYSE BIVARIEE
  ##
  ####################################
  
  
  
  ####################################
  ##
  ##      Biplot Var quanti
  ##
  ####################################
  
  
  rangesBiplot <- reactiveValues(x = NULL, y = NULL)
  
  output$biplot_info <- renderInfoBox({
   
    if(!is.null(input$biplot_hover)){
      data = dataInput()
      
      var1 = input$VariableSelectBi1
      var2 = input$VariableSelectBi2
      
      ind  = which(colnames(data)%in%c(var1,var2))
      data2 = data.frame(x=data[,ind[1]],y=data[,ind[2]])
      
      hover=input$biplot_hover
      dist=sqrt((hover$x-data2$x)^2+(hover$y-data2$y)^2)
      
      if(min(dist)<3)
      {
        sel = which.min(dist)
        infoBox("", HTML(paste('<b>',rownames(data)[sel],'</b>',"<br/> <FONT size='2'> ",var1,":",data2$x[sel],'<br/>',var2,":",data2$y[sel], "</FONT>")), 
                icon = icon("info"),color = "light-blue")
        
      }
      else return(
        infoBox("", HTML("<em> <FONT size='2'> Passez votre souris sur les points pour plus d'information </FONT> </em>"), 
                icon = icon("info"),color = "light-blue")
      )
      
    }
    else return(
      infoBox("", HTML("<em> <FONT size='2'> Passez votre souris sur les points pour plus d'information </FONT> </em>"), 
              icon = icon("info"),color = "light-blue")
    )
    
    
  })
  
 
  
  
  output$biplot <- renderPlot({
    
    data = dataInput()
    
    generateBiPlot(input,rangesBiplot,data)
  })
  
  
  
  output$dymMenu <- renderMenu({
    
    input$RunAnalyse
    data=isolate(dataInput())
    
    if(nrow(data)>0)
    {
      sidebarMenu(
        menuItem("Gestion de données", icon = icon("table"), tabName = "datatable"),
        menuItem("Analyse statistique",
                          menuSubItem("Description générale",tabName="DescGene"),
                          menuSubItem("Analyse univariée",tabName="AnalUni"),
                          menuSubItem("Analyse bivariée",tabName="AnalBi"),
                          icon = icon("bar-chart-o"), tabName = "AnaStat")
      )
    }
  })
  
  observeEvent(input$biplot_dblclick, {
    brush <- input$biplot_brush
    if (!is.null(brush)) {
      rangesBiplot$x <- c(brush$xmin, brush$xmax)
      rangesBiplot$y <- c(brush$ymin, brush$ymax)
      
    } else {
      rangesBiplot$x <- NULL
      rangesBiplot$y <- NULL
    }
  })
  
  
  
  output$SelectDataVarBi1 <- renderUI({
    input$RunAnalyse
    data=dataInput()
    
    if(!is.null(data)){
      selectInput("VariableSelectBi1",p(strong("Choisissez les variables à étudier")),
                colnames(data),
                selected=colnames(data)[1])
      }
  })
  
  
  output$SelectDataVarBi2 <- renderUI({
    input$RunAnalyse
    data=dataInput()

    if(!is.null(data))  selectInput("VariableSelectBi2","vs",colnames(data),selected=colnames(data)[2])
  
    })
    
  
  output$RadioBiPlot <- renderUI({
    
    data = dataInput()
    
    var1 = input$VariableSelectBi1
    var2 = input$VariableSelectBi2
    
    ind  = which(colnames(data)%in%c(var1,var2))
    
    if(length(ind)==2)
    {
      radioButtons("RadioPlotBi","Choisir une représentation",
                       choices=c("Nuage de points"="Nuage","Boxplots"="box","Densités"="densities","Histogramme"="hist"))
    }
    
  })
  
  
 
  
  ## Export graph
  
  output$exportPDFbi <- downloadHandler(
    filename <- function() { paste(input$RadioPlotBi,'stat2E_plot.pdf',sep="_")},
    content <- function(file) {
      pdf(file, width = 6, height = 4)
      plot(generateBiPlot(input,rangesBiplot,dataInput()))
      dev.off()
    }
  )
  
  output$exportPNGbi <- downloadHandler(
    filename <- function() { paste(input$RadioPlotBi,'stat2E_plot.png',sep="_") },
    content <- function(file) {
      png(file, width = 600, height = 400)
      plot(generateBiPlot(input,rangesBiplot,dataInput()))
      dev.off()
    }
  )
  
  output$exportPDFuni <- downloadHandler(
    filename <- function() { paste(input$RadioPlotUni,'stat2E_plot.pdf',sep="_")},
    content <- function(file) {
      pdf(file, width = 6, height = 4)
      plot(generateUniPlot(input,dataInput()))
      dev.off()
    }
  )
  
  output$exportPNGuni <- downloadHandler(
    filename <- function() { paste(input$RadioPlotUni,'stat2E_plot.png',sep="_") },
    content <- function(file) {
      png(file, width = 600, height = 400)
      plot(generateUniPlot(input,dataInput()))
      dev.off()
    }
  )
  
  
  
})