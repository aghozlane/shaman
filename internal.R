


  GetDataFromBIOM <-function(dataBIOM)
  {
    
    counts = biom_data(dataBIOM)
    taxo = observation_metadata(dataBIOM)
    return(list(counts=counts,taxo=taxo))
  }
  
  
  GetDataFromCT <-function(dataC,dataT)
  {
    
    counts = dataC
    taxo = dataT
    return(list(counts=counts,taxo=taxo))
  }
  
  GetInteraction2 <- function(target)
  { 
    res=c()
    namesTarget = colnames(target)[2:ncol(target)]
    k=1
    for(i in 1:(length(namesTarget)-1))
    { 
      for(j in (i+1):length(namesTarget))
      { 
        res[k] = paste(namesTarget[i],":",namesTarget[j],sep="")
        k = k+1
      }
    }
    
    return(res)
  }
  



  PrintContrasts <- function (coefs, contrasts,contnames) 
  {
    contrasts = as.matrix(contrasts)
    out <-""
    
    for (i in 1:ncol(contrasts)) 
    {
      contrast <- contrasts[,i]
      contrast <- paste(ifelse(contrast > 0, "+ ", ""), contrast, sep = "")
      contrast <- gsub("( 1)|(1)", "", contrast)
      out = paste(out,paste("<b>",contnames[i], ":</b> <br/>", paste(contrast[contrast != 0], coefs[contrast != 0], collapse = " ", sep = " ")),"<br/>")
    }
    return(out)
    
  }
  
  ## Get the counts for the selected taxonomy
  GetCountsMerge <- function(dataInput,taxoSelect)
  {
    
    CT = dataInput$counts
    taxo = dataInput$taxo
    ordOTU = order(rownames(taxo))
    indOTU_annot = which(rownames(CT)%in%rownames(taxo))
    counts_annot = CT[indOTU_annot[ordOTU],]
    taxoS = taxo[ordOTU,taxoSelect]
    
    counts = aggregate(counts_annot,by=list(Taxonomy = taxoS),sum)
    
    rownames(counts)=counts[,1];counts=counts[,-1]
    return(counts)
  }


  ## Get the dds object of DESeq2
  Get_dds_object <- function(input,counts,target,design)
  {
   
    dds <- DESeqDataSetFromMatrix(countData=counts, colData=target, design=design)
    ## Size factor estimation
    dds <- estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)))
 
    dds <- estimateDispersions(dds, fitType=input$fitType)
    dds <- nbinomWaldTest(dds)
    return(list(dds = dds,counts=counts,target=target,design=design))
  }

  ## Get the design according to the input
  GetDesign <- function(input)
  {
    InterVar = input$InterestVar
    Interaction = input$Interaction2
    alltmp = c(InterVar,Interaction)
    design = as.formula(paste("~", paste0(alltmp, collapse= "+")))

    return(design)
  }
  


  ## Diagnostic Plots
  Plot_diag <- function(input,resDiff)
  {
    colors = c("dodgerblue","firebrick1","MediumVioletRed","SpringGreen")
    VarInt = input$VarInt
    dds = resDiff$dds
    counts = resDiff$counts
    target = resDiff$target
    group = as.data.frame(target[,VarInt])
    
    if(input$DiagPlot=="barplotTot") barplotTot(counts,group = group, col=colors)
    if(input$DiagPlot=="barplotNul") barPlotNul(counts, group = group, col=colors)
    if(input$DiagPlot=="densityPlot") densityPlotTot(counts, group = group, col=colors)
    if(input$DiagPlot=="MajTax") majTaxPlot(counts, group = group, col=colors)
    if(input$DiagPlot=="SERE") SEREplot(counts)
    if(input$DiagPlot=="Sfactors") diagSFactors(input,dds,frame=1) 
    if(input$DiagPlot=="SfactorsVStot") diagSFactors(input,dds,frame=2) 
    #if(input$DiagPlot=="pcaPlot") barplotTot(counts, group=target[,c(varInt1,varInt2)], col=colors)
  }



  ## barplot total
  barplotTot <- function(counts, group, cex.names = 1, col = c("lightblue","orange", "MediumVioletRed", "SpringGreen")) 
  {
    ncol1 <- ncol(group) == 1
    barplot(colSums(counts), cex.names = cex.names, main = "Total read count per sample", ylab = "Total read count", 
            ylim = c(0, max(colSums(counts)) * 1.2), density = if (ncol1) {NULL}
            else {15}, 
            angle = if (ncol1) {NULL}
            else {c(-45, 0, 45, 90)[as.integer(group[, 2])]}, col = col[as.integer(group[, 1])], las = 2)
    legend("topright", levels(group[, 1]), fill = col[1:nlevels(group[,1])], bty = "n")
    if (!ncol1)  legend("topleft", levels(group[, 2]), density = 15,col = 1, angle = c(-45, 0, 45, 90)[1:nlevels(group[, 2])], bty = "n")
  
  }


  ## barplot Nul 
  barPlotNul <-function (counts, group, cex.names = 1, col = c("lightblue","orange", "MediumVioletRed", "SpringGreen")) 
  {
    
    percentage <- apply(counts, 2, function(x) {sum(x == 0)}) * 100/nrow(counts)
    percentage.allNull <- (nrow(counts) - nrow(removeNul(counts))) * 100/nrow(counts)
    ncol1 <- ncol(group) == 1
    
    
    barplot(percentage, las = 2, col = col[as.integer(group[,1])], 
            density = if (ncol1) {NULL}
            else {15}, 
            angle = if (ncol1) {NULL}
            else {c(-45, 0, 45, 90)[as.integer(group[, 2])]},
            cex.names = cex.names, ylab = "Proportion of null counts", 
            main = "Proportion of null counts per sample", 
            ylim = c(0, 1.2 * ifelse(max(percentage) == 0, 1, max(percentage))))
    
    abline(h = percentage.allNull, lty = 2, lwd = 2)
    legend("topright", levels(group[, 1]), fill = col[1:nlevels(group[,1])], bty = "n")
    if (!ncol1) legend("topleft", levels(group[, 2]), density = 15, col = 1, angle = c(-45, 0, 45, 90)[1:nlevels(group[, 2])], bty = "n")
  }


  ## Plot density
  densityPlotTot <-function (counts, group, col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen")) 
  {
    
    counts <- removeNulCounts(counts)
    ncol1 <- ncol(group) == 1
    plot(density(log2(counts[, 1] + 1)), las = 1, lwd = 2, main = "Density of counts distribution", 
         xlab = expression(log[2] ~ (raw ~ count + 1)), 
         ylim = c(0, max(apply(counts, 2, function(x) {max(density(log2(x + 1))$y)})) * 1.05), 
         lty = if (ncol1) {1}
         else{c(1, 2, 3, 4)[as.integer(group[, 2])[1]]}, 
         col = col[as.integer(group[, 1])[1]])
    
    for (i in 2:ncol(counts)) 
    {
      lines(density(log2(counts[, i] + 1)), col = col[as.integer(group[,1])[i]], lwd = 2, 
            lty = if (ncol1) {1}
            else {c(1, 2, 3, 4)[as.integer(group[, 2])[i]]})
    }
    legend("topright", levels(group[, 1]), lty = 1, col = col[1:nlevels(group[,1])], lwd = 2, bty = "n")
    if (!ncol1) legend("topleft", levels(group[, 2]), lty = c(1, 2, 3, 4)[1:nlevels(group[, 2])], col = 1, lwd = 2, bty = "n")
    
  }


  ## Table of maj. taxo
  majTab <- function(counts,n)
  {
    seqnames <- apply(counts, 2, function(x) {
      x <- sort(x, decreasing = TRUE)
      names(x)[1:n]
    })
    seqnames <- unique(unlist(as.character(seqnames)))
    sum <- apply(counts, 2, sum)
    counts <- counts[seqnames, ]
    sum <- matrix(sum, nrow(counts), ncol(counts), byrow = TRUE)
    p <- round(100 * counts/sum, digits = 3)
    return(p)
  }


  ## Plot maj. taxo
  majTaxPlot <-function (counts, n = 3, group, cex.names = 1, col = c("lightblue",  "orange", "MediumVioletRed", "SpringGreen")) 
  {
    p = majTab(counts,n)
    maj <- apply(p, 2, max)
    seqname <- rownames(p)[apply(p, 2, which.max)]
    ncol1 <- ncol(group) == 1
    x <- barplot(maj, col = col[as.integer(group[, 1])], main = "Proportion of reads from\nmost expressed sequence",
                 ylim = c(0, max(maj) * 1.2), cex.main = 1, 
                 cex.names = cex.names, las = 2, ylab = "Proportion of reads", 
                 density = if (ncol1) {NULL}
                 else {15}, 
                 angle = if (ncol1) {NULL}
                 else {c(-45, 0, 45, 90)[as.integer(group[, 2])]})
    
    legend("topright", levels(group[, 1]), fill = col[1:nlevels(group[,1])], bty = "n")
    if (!ncol1) legend("topleft", levels(group[, 2]), density = 15, col = 1, 
                       angle = c(-45, 0, 45, 90)[1:nlevels(group[, 2])], bty = "n")
    
    for (i in 1:length(seqname)) text(x[i], maj[i]/2, seqname[i], cex = 0.8, srt = 90, adj = 0)
  }
  

  ## plot SERE Coefs
  SEREplot<-function (counts) 
  {
    sere = SEREcoef(counts)
    hc <- hclust(as.dist(sere), method = "ward.D")
    plot(hc, las = 2, hang = -1, xlab = "SERE distance, Ward criterion",main = "Cluster dendrogram\non SERE values")
    
  }
  
  
  ## Get the SERE COEF
  SEREcoef<-function(counts)
  {
    sere <- matrix(NA, ncol = ncol(counts), nrow = ncol(counts))
    for (i in 1:ncol(counts)) {
      for (j in 1:ncol(counts)) {
        sere[i, j] <- sigfun_Pearson_meta(counts[, c(i, j)])
      }
    }
    colnames(sere) <- rownames(sere) <- colnames(counts)
    sere <- round(sere, digits = 3)
    
    return(sere) 
  }
  

  ## function for the SERE coef
  sigfun_Pearson_meta <- function(observed) {
    laneTotals <- colSums(observed)
    total <- sum(laneTotals)
    fullObserved <- observed[rowSums(observed) > 0, ]
    fullLambda <- rowSums(fullObserved)/total
    fullLhat <- fullLambda > 0
    fullExpected <- outer(fullLambda, laneTotals)
    fullKeep <- which(fullExpected > 0)
    oeFull <- (fullObserved[fullKeep] - fullExpected[fullKeep])^2/fullExpected[fullKeep]
    dfFull <- length(fullKeep) - sum(fullLhat != 0)
    sqrt(sum(oeFull)/dfFull)
  }


  ## Plots of size factors
  diagSFactors<-function (input,dds,frame=1) 
  {
    geomeans <- exp(rowMeans(log(counts(dds))))
    samples <- colnames(counts(dds))
    counts.trans <- log2(counts(dds)/geomeans)
    xmin <- min(counts.trans[is.finite(counts.trans)], na.rm = TRUE)
    xmax <- max(counts.trans[is.finite(counts.trans)], na.rm = TRUE)
    
    if(!is.na(input$NbcolSfactors)) parCols = as.numeric(input$NbcolSfactors)
    else parCols = ceiling(ncol(counts.trans)/3)
    
    parRows = ceiling(ncol(counts.trans)/parCols)

    if(frame==1)
    {
      par(mfrow=c(parRows,parCols))
      for (j in 1:ncol(dds)) {
        hist(log2(counts(dds)[, j]/geomeans), nclass = 100, 
             xlab = expression(log[2] ~ (counts/geometric ~ mean)), las = 1, xlim = c(xmin, xmax), 
             main = paste("Size factors diagnostic - Sample ",samples[j], sep = ""), col = "skyblue")
        
        abline(v = log2(sizeFactors(dds)[j]), col = "red", lwd = 1.5)
      }
    }
    
    if(frame==2)
    {
      plot(sizeFactors(dds), colSums(counts(dds)), pch = 19, las = 1, 
           ylab = "Total number of reads", xlab = "Size factors", 
           main = "Diagnostic: size factors vs total number of reads")
      abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0), lty = 2, col = "grey")
    }
  }


  ############################################################
  ##
  ##              CREATE THE CONTRAST DATABASE
  ##
  ############################################################

  
  BaseContrast <- function(input,names,namesfile)
  {  

    v_tmp = c()
    filesize = file.info(namesfile)[,"size"]
    
    for(i in 1:length(names))
    {  
      Tinput = paste("input$",names[i],sep="")
      expr=parse(text=Tinput)
      val = eval(expr) 
      v_tmp[i] = as.numeric(val)
    }
    
    if(filesize!=0)
    { 
      oldContrast = read.table(namesfile,header=TRUE)
      colnamesTmp = c(colnames(oldContrast),input$ContrastName)
      mat = cbind(oldContrast,v_tmp)
    }
    else{ colnamesTmp = input$ContrastName; mat = v_tmp}
    
    write.table(mat,namesfile,row.names=FALSE,col.names = colnamesTmp)
  }
  
  
  ## Remove nul counts
  removeNulCounts <-function (counts) 
  {
    return(counts[rowSums(counts) > 0, ])
  }



DescriptiveStat <-function(vect)
{
  
  nbmiss = length(which(is.na((vect))))
  nbval = length(vect) - nbmiss
  sum = sum(vect,na.rm=TRUE)
  moy = mean(vect,na.rm=TRUE)
  var = var(vect,na.rm=TRUE)
  sd = sd(vect,na.rm=TRUE)
  CV = sd/moy
  stat = c(nbval,nbmiss,sum, summary(vect),var,sd,CV)
  
  names(stat) = c("Nb valeurs","Nb manquants","Somme","Min",
                  "1er Quartile","Mediane",'Moyenne',"3eme Quartile","Max","Variance","Ecart-type","Coeff Variation")
    
  return(stat)
}

DescriptiveStatQuali <-function(dataQuali,namesQuali,indic)
{
  
  res = matrix(0,ncol=5)
  
  for(i in 1:ncol(dataQuali))
  {
    datatmp = dataQuali[,i]
    tabQuali = table(datatmp)
    res2 = c(names(tabQuali),"Total")
    res1 = c(namesQuali[i],rep(NA,length(res2)-1))
    res3 = c(tabQuali,sum(tabQuali))
    res4 = signif(c(tabQuali/sum(tabQuali),1)*100,2)
    res5 = c(signif(cumsum(tabQuali/sum(tabQuali))*100,2),NA)
    res = rbind(res,cbind(res1,res2,res3,res4,res5))
  }
  res=res[-1,]
  rownames(res)=res[,1]
  res=res[,-1]
  colnames(res) = c("Modalites","Effectifs","%","% cumules")
  res = res[,c(1,indic+1)]
  res=as.data.frame(res)
  
  return(res)
  
}

### Gerate 1D plots
generateUniPlot<-function(input,data)
{
  
  ind = which(colnames(data)%in%input$UniVar)
  
  dataTmp = data[,ind]
  TestNum = is.numeric(dataTmp)
  dataTmp = data.frame(x=dataTmp)
  namesTmp = names(data)[ind]
  gg = NULL
  
  ## Numeric data
  if(TestNum && !is.null(input$RadioPlotUni))
  {
    if(input$RadioPlotUni=="hist")
    {
      gg = ggplot(dataTmp,aes(x=x))  + xlab(namesTmp) + theme_bw()
      if(input$HistDens=="freq") gg = gg + geom_histogram(binwidth = input$binwidth,size=input$SizeQQplot-0.5,fill=input$ColorUniplot,alpha=input$TransAlphaUni/100, color="black",aes(y = ..density..))
      if(input$HistDens=="counts") gg = gg + geom_histogram(binwidth = input$binwidth,size=input$SizeQQplot-0.5,fill=input$ColorUniplot,alpha=input$TransAlphaUni/100, color="black")
      if(input$CheckDens) gg = gg + geom_density(size=input$SizeQQplot-0.5)
      if(input$SensGraph=="Hori") gg = gg + coord_flip()  
    }
    
    if(input$RadioPlotUni=="box")
    {
      gg = ggplot(dataTmp,aes(1,x))  + geom_boxplot(fill=input$ColorUniplot,alpha=input$TransAlphaUni/100,size=input$SizeQQplot)+xlim(c(0,2)) 
      gg = gg + ylab(namesTmp) + xlab("")+ theme_bw()
      if(input$CheckAddPointsBox) gg = gg + geom_jitter()
      if(input$SensGraph=="Hori") gg = gg + coord_flip()
      
    }
    
    if(input$RadioPlotUni=="densities")
    {
      gg = ggplot(dataTmp,aes(x=x))  + theme_bw()
      gg = gg + geom_density(fill=input$ColorUniplot,alpha=input$TransAlphaUni/100,size=input$SizeQQplot-0.5) + xlab(namesTmp)
      if(input$SensGraph=="Hori") gg = gg + coord_flip() 
      
    }
    
    if(input$RadioPlotUni=="qqplot")
    {
      gg = ggplot(dataTmp, aes(sample = x)) + geom_point(stat = "qq",color = input$ColorUniplot,alpha=input$TransAlphaUni/100,size=input$SizeQQplot)
      gg = gg + ylab(namesTmp)+ theme_bw()
    }
    
    return(gg)
  }
  
  # Quali data
  if(!TestNum && !is.null(input$RadioPlotUni))
  {
    if(input$RadioPlotUni=="BarPlot")
    {
      gg = ggplot(dataTmp, aes(x))  + xlab(namesTmp)+ theme_bw()
      gg = gg+  geom_bar(fill=input$ColorUniplot,alpha=input$TransAlphaUni/100,size=input$SizeQQplot,  color="black",width=input$widthBarPlot/100)
      if(input$SensGraph=="Hori") gg = gg + coord_flip()
      if(input$BarCircular)  gg = gg + coord_polar()
    }
    
    if(input$RadioPlotUni=="Pie")
    {
      
      count = table(dataTmp$x)
      dataTmp2 = data.frame(frac = count/sum(count), xUnique = as.factor(names(count)))
      dataTmp2 = dataTmp2[order(dataTmp2$frac.Freq), ]
      dataTmp2$ymax = cumsum(dataTmp2$frac.Freq)
      dataTmp2$ymin = c(0, head(dataTmp2$ymax, n=-1))
      
      dataTmp2$xminPie = 1-input$PieWidth/100
      gg =  ggplot(dataTmp2, aes(fill=xUnique, ymax=ymax, ymin=ymin, xmax=4, xmin=4*xminPie)) + geom_rect(alpha=input$TransAlphaUni/100) + coord_polar(theta="y") 
      gg =  gg + xlim(c(0, 4)) + scale_fill_discrete(name=namesTmp)+ theme_bw()
      
    }
    return(gg)
  }
  
  
  
}


### Gerate 2D plots
generateBiPlot<-function(input,rangesBiplot,data)
{
  
  var1 = input$VariableSelectBi1
  var2 = input$VariableSelectBi2
  
  ind1  = which(colnames(data)%in%c(var1))
  ind2  = which(colnames(data)%in%c(var2))
  ind =  c(ind1,ind2)
  
  Num = sapply(data[,unique(ind)],is.numeric)
  
  if(length(unique(ind))==2)
  {
    if(length(which(Num))==2) data2 = data.frame(x=data[,ind1],y=data[,ind2])
    if(length(which(Num))==1) data2 = data.frame(x=data[,ind[which(Num)]],y=data[,ind[which(!Num)]])        
  }
  
  if(length(unique(ind))==1)  data2 = data.frame(x=data[,ind],y=data[,ind])
  
  if(!is.null(input$RadioPlotBi))
  {
    if(input$RadioPlotBi=="Nuage") 
    {
      gg = ggplot(data2,aes(x,y))  
      gg = gg + geom_point(color=input$ColorBiplot,size=input$SizePoint,alpha=input$TransAlphaBi/100) + theme_bw()
      if(!is.null(rangesBiplot$x) && !is.null(rangesBiplot$x)) gg = gg + xlim(rangesBiplot$x) + ylim(rangesBiplot$y)
      gg = gg + xlab(var1) + ylab(var2)
      if(input$CheckLM) gg = gg + geom_smooth(method='lm')
    }    
    
    if(input$RadioPlotBi=="densities" )
    {
      if(length(which(Num))==2)
      {
        vartmp = c(data[,ind[1]],data[,ind[2]])
        y=c(rep(var1,nrow(data)),rep(var2,nrow(data)))
        data2 = data.frame(x=vartmp,y=y)
      }
      
      gg = ggplot(data2,aes(x=x,fill=y)) 
      gg = gg + xlab("Valeurs")+labs(colour = var2) + scale_fill_discrete(name="Légende")
      gg = gg + geom_density(size=input$SizePoint,alpha = input$TransAlphaBi/100) + theme_bw()
      if(input$SensGraphBi=="Hori") gg = gg + coord_flip()
    }    
    
    if(input$RadioPlotBi=="hist" )
    {
      if(length(which(Num))==2)
      {
        vartmp = c(data[,ind[1]],data[,ind[2]])
        y=c(rep(var1,nrow(data)),rep(var2,nrow(data)))
        data2 = data.frame(x=vartmp,y=y)
      }
      
      gg = ggplot(data2,aes(x=x,fill=y)) 
      gg = gg + xlab("Valeurs")+labs(colour = var2) + scale_fill_discrete(name="Légende")
      gg = gg + geom_histogram(binwidth = input$binwidthBi,size=input$SizePoint-0.5,alpha = input$TransAlphaBi/100,color="black") + theme_bw()
      if(input$SensGraphBi=="Hori") gg = gg + coord_flip()
    }    
    
    if(input$RadioPlotBi=="box" )
    {
      if(length(which(Num))==2)
      {
        vartmp = c(data[,ind[1]],data[,ind[2]])
        y=c(rep(var1,nrow(data)),rep(var2,nrow(data)))
        data2 = data.frame(x=vartmp,y=y)
      }
      labtmp = unique(data2$y)
      gg = ggplot(data2,aes(y,x)) + scale_fill_discrete(name="Légende")
      gg = gg+ geom_boxplot(aes(fill=y),alpha=input$TransAlphaBi/100,size=input$SizePoint-0.5) + xlab("")+ theme_bw()
      if(input$CheckAddPointsBoxBi) gg = gg + geom_jitter(size=input$SizePoint)
      if(input$SensGraphBi=="Hori") gg = gg + coord_flip()
      
    }    
    
    return(gg)
    
  }
  
}


