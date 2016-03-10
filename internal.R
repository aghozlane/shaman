

## Modified version of expand.grid
expand.grid2.list <- function(listInput)
{
  n = length(listInput)
  if(is.list(listInput) && n>1)
  {
    l1 = listInput[[1]]
    l2 = listInput[[2]]
    res = c()
    
    for(i in l1){
      for(j in l2){ 
        res = rbind(res,paste(i,j,sep = "-"))
      }
    }
    listInput[[1]] = res
    listInput = listInput[-2]
    if(length(listInput)>1 && is.list(listInput)) res = expand.grid2.list(listInput)
  }
  else res = listInput
  return(res)
}


## Function for the rdp format
getval <- function(annotation_rdp, interest, threshold_annot){
  annotation_rdp = unlist(strsplit(annotation_rdp,"\t"))
  annotation = c(annotation_rdp[1])
  for(level in interest){
    idlevel=which(annotation_rdp == level)
    if(length(idlevel)>0){
      if(as.numeric(annotation_rdp[idlevel+1]) >= threshold_annot){
        annotation = c(annotation, gsub("\"", "", annotation_rdp[idlevel-1]))
      }
      else annotation = c(annotation, "NA")
    }
    else annotation = c(annotation, "NA")  
  }
  return(annotation)
}

## Read rdp file
read_rdp <- function(filename, threshold_annot)
{
  interest=c("phylum", "class", "order", "family", "genus")
  conn <- file(filename,open="r")
  linn <-readLines(conn)
  tab=t(sapply(1:length(linn), function(i) getval(linn[i], interest, threshold_annot)))
  close(conn)
  
  if(!TRUE%in%duplicated(tab[,1])) rownames(tab)=tab[,1];tab=tab[,-1]
  colnames(tab) = c("Phylum","Class","Order","Family","Genus")
  
  return(tab)
}






CheckCountsTable <- function(counts)
  {
    Error = NULL
    Warning = NULL
    numTest = FALSE%in%sapply(counts,is.numeric)
    if(ncol(counts)<=1){Error = "The number of columns of the counts table must be at least 2" }
    if(nrow(counts)<=1){Error = "The number of rows of the counts table must be at least 2" }
    if(numTest){Error = "The counts table must contain only numeric values" }
    if(!numTest)
    {
      if(0%in%colSums(counts)){Error = "At least one of the column of the counts table is 0" }
      if(min(counts)<0){Error = "The counts table must contain only positive values" }
    }
    if(TRUE%in%sapply(counts,is.na)){Warning = "NA values are considered as 0 is the counts table"; counts[sapply(counts,is.na)]=0}
    
  
    return(list(Error=Error,Warning=Warning,counts=counts))
  }
  
  CheckTaxoTable <- function(taxo,counts)
  {
    Error = NULL
    Warning = NULL
    if(ncol(taxo)<1){Error = "The number of columns of the taxonomy table must be at least 1" }
    if(nrow(taxo)<=1){Error = "The number of rows if the taxonomy table must be at least 2" }
    if(TRUE%in%is.numeric(taxo)){Error = "The taxonomy table must contain only character" }

    for(i in 1:ncol(taxo))
    {
      level = levels(taxo[,i])
      nb = length(level)
      if(nb==1 && level=="NA"){ Error = "At least one column contains only NA"}
    }
    
    ## Annotated features without counts
    if(!any(rownames(taxo)%in%rownames(counts))){ Error = "Some annotated features are not in the count table"}
    
    return(list(Error=Error,Warning=Warning))
  }
  
  PercentAnnot <- function(counts,taxo)
  {
    Error = NULL  
    tmp = table(rownames(counts)%in%rownames(taxo))
    Percent = tmp["TRUE"]/sum(tmp)
    if(is.na(Percent)) Percent = 0
    if(Percent==0){Error = "Counts table and annotation do not matched" }
       
    return(list(Error=Error,Percent=Percent))
  }
  
  
  GetDataFromBIOM <-function(dataBIOM)
  {
    ## Counts table
    counts = biom_data(dataBIOM)
    counts = as.matrix(counts)
    counts = as.data.frame(counts)
    CheckCounts = CheckCountsTable(counts)
    counts = CheckCounts$counts
    
    ## Taxonomy table
    taxo = as.data.frame(observation_metadata(dataBIOM))
    CheckTaxo = CheckTaxoTable(taxo,counts)
    
    ## Pourcentage of annotation
    tmp = PercentAnnot(counts,taxo)
    
    return(list(counts=counts,taxo=taxo,CheckCounts=CheckCounts,CheckTaxo=CheckTaxo,Percent=tmp$Percent,CheckPercent=tmp$Error))
  }
  
  
  GetDataFromCT <-function(dataC,dataT)
  {
    
    ## Counts table
    counts = dataC
    CheckCounts = CheckCountsTable(counts)
    counts = CheckCounts$counts
    
    ## Taxonomy table
    taxo = as.data.frame(dataT)
    CheckTaxo = CheckTaxoTable(taxo,counts)
    
    ## Pourcentage of annotation
    tmp = PercentAnnot(counts,taxo)
    
    return(list(counts=counts,taxo=taxo,CheckCounts=CheckCounts,CheckTaxo=CheckTaxo,Percent=tmp$Percent,CheckPercent=tmp$Error))
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
  


  ## Print the contrasts
  PrintContrasts <- function (coefs, contrasts,contnames) 
  {
    contrasts = as.matrix(contrasts)
    out <-""
    
    for (i in 1:ncol(contrasts)) 
    {
      contrast <- contrasts[,i]
      contrast <- paste(ifelse(contrast > 0, "+ ", ""), contrast, sep = "")
      contrast <- gsub("( 1)|(1)", "", contrast)
      out = paste(out,paste("<b>",contnames[i], "</b> <br/>", paste(contrast[contrast != 0], coefs[contrast != 0], collapse = " ", sep = " ")),"<br/>")
    }
    return(out)
    
  }

  
  
  ## Get the counts for the selected taxonomy
  GetCountsMerge <- function(input,dataInput,taxoSelect,target,design)
  {
    ## Init
    counts= NULL
    CheckTarget = FALSE
    CT_noNorm = NULL
    normFactors = NULL
    FeatureSize = NULL

    ## Counts and taxo tables
    CT = dataInput$counts
    taxo = dataInput$taxo

    ## Select cols in the target
    labels = target[,1]
    ind = which(colnames(CT)%in%labels)
    
    ## Get the feature size for the normalisation
    Size_indcol = which(toupper(colnames(CT))%in%"SIZE")
    if(length(Size_indcol)==1) FeatureSize = CT[,Size_indcol]
    
    if(length(ind)==length(labels))
    { 
      CT = CT[,ind]
      
      ## Order CT according to the target
      CT = OrderCounts(counts=CT,labels=labels)$CountsOrder
      CT_noNorm = CT
      RowProd = sum(apply(CT_noNorm,1,prod))
      
      ## Counts normalisation
      dds <- DESeqDataSetFromMatrix(countData=CT, colData=target, design=design)
      ## Normalisation with or without 0
      if(input$AccountForNA || RowProd==0) dds = estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)),geoMeans=GeoMeansCT(CT))
      if(!input$AccountForNA && RowProd!=0) dds = estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)))
       
      normFactors = sizeFactors(dds)
      # Only interesting OTU
      merged_table = merge(CT, taxo[order(rownames(CT)),], by="row.names")
      CT = merged_table[,2: (dim(CT)[2]+1)]
      taxo = merged_table[,(dim(CT)[2]+2):dim(merged_table)[2]]
      rownames(CT) = merged_table[,1]
      rownames(taxo) = merged_table[,1]
      ordOTU = order(rownames(taxo))
      counts_annot = CT
      
#       ordOTU = order(rownames(taxo))
#       indOTU_annot = which(rownames(CT)%in%rownames(taxo))
#       counts_annot = CT[indOTU_annot[ordOTU],]
     
      if(taxoSelect=="OTU") counts = counts_annot
      else{
        taxoS = taxo[ordOTU,taxoSelect]
        counts = aggregate(counts_annot,by=list(Taxonomy = taxoS),sum)
        rownames(counts)=counts[,1];counts=counts[,-1]
      }
      
      ## Ordering the counts table according to the target labels 
      tmpOrder = OrderCounts(counts,normFactors,labels)
      counts = tmpOrder$CountsOrder
      normFactors = tmpOrder$normFactorsOrder
      CheckTarget = TRUE
    }
    return(list(counts=counts,CheckTarget=CheckTarget,normFactors=normFactors, CT_noNorm=CT_noNorm))
  }
  
  ## Get the geometric mean of the counts (0 are replaced by NA values)
  GeoMeansCT <- function(CT)
  {
    CT=as.matrix(CT)
    CT[which(CT<1)]=NA
    gm = apply(CT,1,geometric.mean,na.rm=TRUE)
    return(gm)
  }
  
  
  ## Order the counts 
  OrderCounts <- function(counts,normFactors=NULL,labels)
  {
    n = length(labels)
    CountsOrder = counts
    normFactorsOrder = normFactors
    for(i in 1:n)
    {
      ind = which(labels[i]==colnames(counts))
      CountsOrder[,i] = counts[,ind]
      if(!is.null(normFactors)) normFactorsOrder[i] = normFactors[ind]
    }
    colnames(CountsOrder) = labels
    return(list(CountsOrder=CountsOrder,normFactorsOrder = normFactorsOrder))
  }
  
  
  ## Get the dds object of DESeq2
  Get_dds_object <- function(input,counts,target,design,normFactorsOTU,CT_noNorm)
  {
    dds <- DESeqDataSetFromMatrix(countData=counts, colData=target, design=design)
    sizeFactors(dds) = normFactorsOTU
    dds <- estimateDispersions(dds, fitType=input$fitType)
    dds <- nbinomWaldTest(dds)
    countsNorm = counts(dds, normalized = TRUE)
    return(list(dds = dds,raw_counts=counts,countsNorm=countsNorm,target=target,design=design,normFactors = normFactorsOTU,CT_noNorm=CT_noNorm))
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
    
    VarInt = input$VarInt
    dds = resDiff$dds
    counts = resDiff$raw_counts
    if(input$CountsType=="norm") counts = resDiff$countsNorm
    target = resDiff$target
    normFactors = resDiff$normFactors
    CT_noNorm = resDiff$CT_noNorm
    group = as.data.frame(target[,VarInt])
    rownames(group) = rownames(target)
    res = NULL
    
    if(ncol(group)>0 && nrow(counts)>0)
    { 
      colors = rep(c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78", "#2ca02c","#98df8a","#d62728","#ff9896","#9467bd","#c5b0d5","#8c564b",
                     "#c49c94","#e377c2","#f7b6d2","#7f7f7f", "#c7c7c7","#bcbd22","#dbdb8d","#17becf","#9edae5"),ceiling(nrow(target)/20))
      
      if(input$DiagPlot=="barplotTot") res = barplotTot(input,counts,group = group, col=colors)
      if(input$DiagPlot=="barplotNul") res = barPlotNul(input,counts, group = group, col=colors)
      if(input$DiagPlot=="densityPlot") res = densityPlotTot(input,counts, group = group, col=colors)
      if(input$DiagPlot=="DispPlot") res = plotDispEsts(dds)
      if(input$DiagPlot=="MajTax") res = majTaxPlot(input,counts, group = group, col=colors)
      if(input$DiagPlot=="SfactorsVStot") res = diagSFactors(input,normFactors,resDiff$raw_counts) 
      if(input$DiagPlot=="pcaPlot") res = PCAPlot_meta(input,dds, group,  type.trans = input$TransType, col = colors)
      if(input$DiagPlot=="pcoaPlot") res = PCoAPlot_meta(input,dds, group, col = colors) 
      if(input$DiagPlot=="clustPlot") res = HCPlot(input,dds,group,type.trans=input$TransType,counts,col=colors)
    }
    
    return(res)
  }

  
  HCPlot <- function (input,dds,group,type.trans,counts,col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen")) 
  {
    
    res = NULL
    
    ## Get the counts
    if (input$DistClust == "euclidean" && type.trans == "VST") counts <- assay(varianceStabilizingTransformation(dds))
    if (input$DistClust == "euclidean" && type.trans == "rlog") counts <- assay(rlogTransformation(dds))
    
    ## Get the group of leaf
    group = apply(group,1,paste, collapse = "-")    
    nb = length(unique((group)))
    
    ## Get the dendrogram
    if(input$DistClust!="sere") dist = vegdist(t(counts), method = input$DistClust)
    if(input$DistClust=="sere") dist = as.dist(SEREcoef(counts))
    hc <- hclust(dist, method = "ward.D")
    
    dend = as.dendrogram(hc)
    
    ## Get the type of dendrogram
    type <- input$typeHculst
    
    dend <- set(dend, "labels_cex", input$cexLabelDiag)
    if(input$colorHC) labels_colors(dend)<-col[as.integer(as.factor(group))][order.dendrogram(dend)]
    if(type=="hori") 
    { 
      par(cex=input$cexTitleDiag,mar=c(6,6,4,5))
      res = plot(dend, main = "Cluster dendrogram",xlab = paste(input$DistClust,"distance, Ward criterion",sep=" "),cex=input$cexLabelDiag)
    }  
    if(type!="hori")
    {
      par(cex=input$cexTitleDiag,mar=c(6,6,4,5))
      res = circlize_dendrogram(dend, labels_track_height = 0.2, dend_track_height = .3, main = "Cluster dendrogram",xlab = paste(input$DistClust,"distance, Ward criterion",sep=" "))
    }
    return(res)
  }
  
  
  ## Color for the horizontal dendro
  colLabdendo <- function(n,group) {
    
    group = apply(group,1,paste, collapse = "-")
    
    nb = length(unique((group)))
    namesGrp = names(group)

    if (is.leaf(n)) {
      a <- attributes(n)
      labCol <- rainbow(nb)[as.integer(as.factor(group))[which(namesGrp == a$label)]]
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    return(n)
  }
  
  ## Diagnostic Plots Eigen value
  Plot_diag_Eigen <- function(input,resDiff)
  {
    colors = c("dodgerblue","firebrick1","MediumVioletRed","SpringGreen")
    VarInt = input$VarInt
    dds = resDiff$dds
    counts = resDiff$counts
    target = resDiff$target
    group = as.data.frame(target[,VarInt])
    
    ## If more than 4 levels for one factor
    maxFact =max(sapply(group,FUN=function(x) length(levels(x))))
    if(maxFact>=4) colors = rainbow(maxFact) 
    
    PCAPlot_meta(input,dds, group,  type.trans = input$TransType, col = colors, plot = "eigen") 
  }
  
  Plot_diag_pcoaEigen = function(input,resDiff)
  {
    colors = c("SpringGreen","dodgerblue","black","firebrick1")
    VarInt = input$VarInt
    dds = resDiff$dds
    target = resDiff$target
    group = as.data.frame(target[,VarInt])
    rownames(group) = rownames(target)
    PCoAPlot_meta(input,dds, group, col = colors, plot = "eigen") 
  }
  
  
  

  ## barplot total
  barplotTot <- function(input,counts, group, cex.names = 1, col = c("lightblue","orange", "MediumVioletRed", "SpringGreen")) 
  {

    ncol1 <- ncol(group) == 1
    par(cex=input$cexTitleDiag,mar=c(6,6,4,5))
    barplot(colSums(counts), cex.names = cex.names, main = "Total mapped read count per sample", ylab = "Total mapped read count", 
            ylim = c(0, max(colSums(counts)) * 1.2), density = if (ncol1) {NULL}
            else {15}, 
            angle = if (ncol1) {NULL}
            else {seq(0,160,length.out =nlevels(group[, 2]))[as.integer(group[, 2])]}, col = col[as.integer(group[, 1])], las = 2)
    legend("topright", levels(group[, 1]), fill = col[1:nlevels(group[,1])], bty = "n")
    if (!ncol1)  legend("topleft", levels(group[, 2]), density = 15,col = 1, angle = seq(0,160,length.out =nlevels(group[, 2]))[1:nlevels(group[, 2])], bty = "n")
  
  }


  ## barplot Nul 
  barPlotNul <-function (input,counts, group, cex.names = 1, col = c("lightblue","orange", "MediumVioletRed", "SpringGreen")) 
  {
    
    percentage <- apply(counts, 2, function(x) {sum(x == 0)}) * 100/nrow(counts)
    percentage.allNull <- (nrow(counts) - nrow(removeNulCounts(counts))) * 100/nrow(counts)
    ncol1 <- ncol(group) == 1
    
    par(cex=input$cexTitleDiag,mar=c(6,6,4,5))

    barplot(percentage, las = 2, col = col[as.integer(group[,1])], 
            density = if (ncol1) {NULL}
            else {15}, 
            angle = if (ncol1) {NULL}
            else {seq(0,160,length.out =nlevels(group[, 2]))[as.integer(group[, 2])]},
            cex.names = cex.names, ylab = "Proportion of null counts", 
            main = "Proportion of null counts per sample", 
            ylim = c(0, 1.2 * ifelse(max(percentage) == 0, 1, max(percentage))))
    
    abline(h = percentage.allNull, lty = 2, lwd = 2)
    legend("topright", levels(group[, 1]), fill = col[1:nlevels(group[,1])], bty = "n")
    if (!ncol1) legend("topleft", levels(group[, 2]), density = 15, col = 1, angle = seq(0,160,length.out =nlevels(group[, 2]))[1:nlevels(group[, 2])], bty = "n")
  }


  ## Plot density
  densityPlotTot <-function (input,counts, group, col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen")) 
  {
    
    counts <- removeNulCounts(counts)
    ncol1 <- ncol(group) == 1
    par(cex=input$cexTitleDiag,mar=c(8,5,4,5))
    plot(density(log2(counts[, 1] + 1)), las = 1, lwd = 2, main = "Density of counts distribution", 
         xlab = expression(log[2] ~ (raw ~ count + 1)), 
         ylim = c(0, max(apply(counts, 2, function(x) {max(density(log2(x + 1))$y)})) * 1.05), 
         lty = if (ncol1) {1}
         else{rep(seq(1:6),ceiling(nlevels(group[, 2])/6))[as.integer(group[, 2])[1]]}, 
         col = col[as.integer(group[, 1])[1]])
    
    for (i in 2:ncol(counts)) 
    {
      lines(density(log2(counts[, i] + 1)), col = col[as.integer(group[,1])[i]], lwd = 2, 
            lty = if (ncol1) {1}
            else{rep(seq(1:6),ceiling(nlevels(group[, 2])/6))[as.integer(group[, 2])[i]]})
    }
    legend("topright", levels(group[, 1]), lty = 1, col = col[1:nlevels(group[,1])], lwd = 2, bty = "n")
    if (!ncol1) legend("topleft", levels(group[, 2]), lty = rep(seq(1:6),ceiling(nlevels(group[, 2])/6))[1:nlevels(group[, 2])], col = 1, lwd = 2, bty = "n")
    
  }


  ## Table of maj. taxo
  majTab <- function(input,counts,n)
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
  majTaxPlot <-function (input,counts, n = 3, group, cex.names = 1, col = c("lightblue",  "orange", "MediumVioletRed", "SpringGreen")) 
  {
    p = majTab(input,counts,n)
    maj <- apply(p, 2, max)
    seqname <- rownames(p)[apply(p, 2, which.max)]
    ncol1 <- ncol(group) == 1
    par(cex=input$cexTitleDiag,mar=c(6,6,4,5))
    x <- barplot(maj, col = col[as.integer(group[, 1])], main = "Proportion of mapped reads from\nmost expressed sequence",
                 ylim = c(0, max(maj) * 1.2), cex.main = 1, 
                 cex.names = cex.names, las = 2, ylab = "Proportion of mapped reads", 
                 density = if (ncol1) {NULL}
                 else {15}, 
                 angle = if (ncol1) {NULL}
                 else {seq(0,160,length.out =nlevels(group[, 2]))[as.integer(group[, 2])]})
    
    legend("topright", levels(group[, 1]), fill = col[1:nlevels(group[,1])], bty = "n")
    if (!ncol1) legend("topleft", levels(group[, 2]), density = 15, col = 1, 
                       angle = seq(0,160,length.out =nlevels(group[, 2]))[1:nlevels(group[, 2])], bty = "n")
    
    for (i in 1:length(seqname)) text(x[i], maj[i]/2, seqname[i], cex=input$cexLabelDiag, srt = 90, adj = 0)
  }
  

  ## plot SERE Coefs
#   SEREplot<-function(input,counts) 
#   {
#     sere = SEREcoef(counts)
#     hc <- hclust(as.dist(sere), method = "ward.D")
#     plot(hc, las = 2, hang = -1, xlab = "SERE distance, Ward criterion",main = "Cluster dendrogram\non SERE values")
#     
#   }
  
  
  ## Get the SERE COEF
  SEREcoef<-function(counts)
  {
    counts = as.matrix(counts)
    sere <- matrix(0, ncol = ncol(counts), nrow = ncol(counts))
    for (i in 1:(ncol(counts)-1)) {
      for (j in (i+1):ncol(counts)) {
        sere[i, j] <- sigfun_Pearson_meta(counts[, c(i, j)])
      }
    }
    sere=sere+t(sere)
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
  diagSFactors<-function (input,normFactors,counts) 
  {
    geomeans <- exp(rowMeans(log(counts)))
    samples <- colnames(counts)
      par(cex=input$cexTitleDiag,mar=c(6,6,4,5))
      plot(normFactors, colSums(counts), pch = 19, las = 1,cex = ifelse(input$addLabelSFact,0,input$cexLabelDiag),
           ylab = "Total number of reads", xlab = "Size factors", 
           main = "Diagnostic: size factors vs total number of reads")
      if(input$addLabelSFact) text(normFactors,colSums(counts),labels = samples,cex=input$cexLabelDiag)
      abline(lm(colSums(counts) ~ normFactors + 0), lty = 2, col = "grey")
  }

  
  ### PCoA
  PCoAPlot_meta <-function (input,dds, group_init,col = c("SpringGreen","dodgerblue","black","firebrick1"), plot = "pcoa") 
  {
    cval=c()
    time_set = 0
    # Set of shape
    shape=c(19,17,15,18)
    
    ## Var of interest
    VarInt  = input$VarInt
    
    ## Group
    group = as.character(apply(group_init,1,paste, collapse = "-"))
    
    ## Keep only some sample 
    val = c()
    for(i in 1:length(VarInt))
    { 
      Tinput = paste("input$","Mod",VarInt[i],sep="")
      expr=parse(text=Tinput)
      ## All the modalities for all the var of interest
      val = c(val,eval(expr))
    }
    if(length(VarInt)>1) Kval = apply(expand.grid(val,val),1,paste, collapse = "-")
    else Kval = val
    ind_kept = which(as.character(group)%in%Kval)
    
    ## Get the group corresponding to the modalities
    group = group[ind_kept]
    nb = length(unique((group)))
    group = as.factor(group)
    
    if(nlevels(group)!=0)
    { 
      ## Get the norm data
      counts.norm = as.data.frame(round(counts(dds)))
      if(input$CountsType=="norm") counts.norm = as.data.frame(round(counts(dds, normalized = TRUE)))
      # was removed
      counts.norm = counts.norm[,ind_kept]
  
      ## Get the distance
      if(input$DistClust!="sere") dist.counts.norm = vegdist(t(counts.norm), method = input$DistClust)
      if(input$DistClust=="sere") dist.counts.norm = as.dist(SEREcoef(counts.norm))
      
      ## Do PCoA
      pco.counts.norm = dudi.pco(d = dist.counts.norm, scannf = FALSE,nf=3)
      
      ## Get eigen values
      eigen=(pco.counts.norm$eig/sum(pco.counts.norm$eig))*100
      
      ## xlim and ylim of the plot
      min = min(pco.counts.norm$li)
      max = max(pco.counts.norm$li)
      
      ## get condition set
      condition_set=val[which(val %in% unique(group_init$condition))]
      time_set=val[which(val %in% unique(group_init$time))]
      
      ## Colors
      if(length(col)<length(condition_set) * length(time_set))# && !input$colorgroup)
      {
        col = rainbow(length(condition_set) * length(time_set))
      }
      #else if(length(col)<length(condition_set) * length(time_set) && input$colorgroup){
      #  col = rep(col[1:length(condition_set)], length(time_set))
      #}
      if (length(time_set) == 1 && length(condition_set) <= 4){
        cval = apply(expand.grid(condition_set,time_set),1,paste, collapse = "-")
        cval = sort(cval)
      }
      
      # to reactivate
      #pco.counts.norm$li = pco.counts.norm$li[ind_kept,]
      if (plot == "pcoa"){
        par(cex=input$cexTitleDiag,mar=c(6,6,4,5))
        ## Plot axis, label and circles
        plot(pco.counts.norm$li[1:2], xlab=paste("PC1 : ",round(eigen[1],1),"%") , ylab=paste("PC2 : ",round(eigen[2],1),"%"),
             xlim=c(min+0.25*min,max+0.25*max), ylim=c(min-0.1,max+0.1), cex.axis=1, cex.lab=1,lwd=2, type="n",main='Principal Coordinates Analysis ')
        # Set different shapes
        if(input$labelPCOA == "Group"){
          if(!is.null(cval)){
            for (i in 1:length(cval)){
              points(pco.counts.norm$li[which(group==cval[i]),1:2],pch=shape[i],col=col[i], cex=input$cexpoint)
            }
            s.class(dfxy = pco.counts.norm$li, fac = group, col = col, label = levels(group),
                    add.plot = TRUE, cpoint = 0, cell=input$cexcircle, clabel=input$cexLabelDiag,  cstar = input$cexstar)
          }else s.class(dfxy = pco.counts.norm$li, fac = group, col = col, label = levels(group),
                        add.plot = TRUE, cpoint = input$cexpoint, cell=input$cexcircle, clabel=input$cexLabelDiag,  cstar = input$cexstar)
        }  
        else{
          s.label(pco.counts.norm$li, clabel = input$cexLabelDiag,boxes=FALSE, add.plot = TRUE)
          s.class(dfxy = pco.counts.norm$li, fac = group, col = col, label = levels(group), add.plot = TRUE, cpoint = 0, clabel = 0, cstar = input$cexstar, cell=input$cexcircle)
        }
      }else{
        barplot(eigen[1:7], xlab="Dimensions", ylab="Eigenvalues (%)", names.arg = 1:7, col = c(rep("black", 2), rep("grey", 5)), ylim=c(0,max(eigen)+5), cex.axis=1.2, cex.lab=1.4,cex.names=1.2)
      }
  }
  
  }
  
  ### PCA
  PCAPlot_meta <-function(input,dds, group_init, n = min(500, nrow(counts(dds))), type.trans = c("VST", "rlog"), 
                           col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen"),plot="pca") 
  {
    ## Var of interest
    VarInt  = input$VarInt
    
    group = as.character(apply(group_init,1,paste, collapse = "-"))
    group_init = group
    
    ## Keep only some sample 
    val = c()
    for(i in 1:length(VarInt))
    { 
      Tinput = paste("input$","Mod",VarInt[i],sep="")
      expr=parse(text=Tinput)
      ## All the modalities for all the var of interest
      val = c(val,eval(expr))
    }
    if(length(VarInt)>1) Kval = apply(expand.grid(val,val),1,paste, collapse = "-")
    else Kval = val
    ind_kept = which(as.character(group)%in%Kval)
    
    ## Get the group corresponding to the modalities
    group = group[ind_kept]
    nb = length(unique((group)))
    group = as.factor(group)
    
    ## To select the colors
    indgrp =as.integer(as.factor(group_init))[ind_kept]
    
    
    if(nlevels(group)!=0)
    { 
      type.trans <- type.trans[1]
      
      if (type.trans == "VST") counts.trans <- assay(varianceStabilizingTransformation(dds))
      else counts.trans <- assay(rlogTransformation(dds))
      counts.trans = counts.trans[,ind_kept]
      
      rv = apply(counts.trans, 1, var, na.rm = TRUE)
      pca = prcomp(t(counts.trans[order(rv, decreasing = TRUE),][1:n, ]))
      
      if(plot=="pca")
      { 
        prp <- pca$sdev^2 * 100/sum(pca$sdev^2)
        prp <- round(prp, 2)
        ncol1 <- ncol(group) == 1
        
        abs = range(pca$x[, 1])
        abs = abs(abs[2] - abs[1])/25
        ord = range(pca$x[, 2])
        ord = abs(ord[2] - ord[1])/25
        
        par(mfrow = c(1, 2),cex=input$cexTitleDiag,mar=c(6,6,4,5))
        plot(pca$x[, 1], pca$x[, 2], las = 1, cex = input$cexTitleDiag, col = col[indgrp], 
             pch = 16,
             xlab = paste0("PC1 (", prp[1], "%)"),
             ylab = paste0("PC2 (", prp[2], "%)"), 
             main = "Principal Component Analysis"
              )
        abline(h = 0, v = 0, lty = 2, col = "lightgray")
        text(pca$x[, 1] - ifelse(pca$x[, 1] > 0, abs, -abs), pca$x[,2] - ifelse(pca$x[, 2] > 0, ord, -ord), colnames(counts.trans), col = col[indgrp],cex=input$cexLabelDiag)
        abs = range(pca$x[, 1])
        abs = abs(abs[2] - abs[1])/25
        ord = range(pca$x[, 3])
        ord = abs(ord[2] - ord[1])/25
        plot(pca$x[, 1], pca$x[, 3], las = 1, cex = input$cexTitleDiag, col = col[indgrp], 
             pch = 16,
             xlab = paste0("PC1 (", prp[1], "%)"), 
             ylab = paste0("PC3 (", prp[3], "%)"), 
             main = "Principal Component Analysis")
        abline(h = 0, v = 0, lty = 2, col = "lightgray")
        text(pca$x[, 1] - ifelse(pca$x[, 1] > 0, abs, -abs), pca$x[,3] - ifelse(pca$x[, 3] > 0, ord, -ord), colnames(counts.trans), col = col[indgrp],cex=input$cexLabelDiag)
      }
      if(plot=="eigen"){eigen = pca$sdev[1:min(7,ncol(counts.trans))]^2; barplot(eigen, xlab="Dimensions", ylab="Eigenvalues (%)", names.arg = 1:min(7,ncol(counts.trans)), col = c(rep("black", 3), rep("grey", 4)), ylim=c(0,max(eigen)+5), cex.axis=1.2, cex.lab=1.4,cex.names=1.2)}
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

  
  ############################################################
  ##
  ##              VISUALISATION PLOTS
  ##
  ############################################################
  
  GetDataToPlot <- function(input,resDiff,VarInt,ind_taxo,aggregate=TRUE)
  {
    dds = resDiff$dds
    val = c()
    list.val = list()
    counts = as.data.frame(round(counts(dds, normalized = TRUE)))
    target = resDiff$target
    counts_tmp_combined = NULL
    prop_tmp_combined = NULL
    targetInt = NULL
    namesCounts = NULL
    levelsMod = NULL
    prop_all=NULL
    ## Select a subset within the taxonomy level (default is the 12 most abundant)
    nbKept = length(ind_taxo)
    Taxonomy = rownames(counts)
    
    if (length(VarInt)>0 && nbKept>0)
    { 
      ## Get the modalities to keep
      for(i in 1:length(VarInt))
      { 
        Tinput = paste("input$","ModVisu",VarInt[i],sep="")
        expr=parse(text=Tinput)
        ## All the modalities for all the var of interest
        val = c(val,eval(expr))
        list.val[[i]] = eval(expr)
      }
      if (!is.null(val) && !is.null(list.val))
      {

        ## Create the variable to plot
        targetInt = as.data.frame(target[,VarInt])
        rownames(targetInt)=target[,1]  
        ## Combining the Varint
        if(length(VarInt)>1){targetInt$AllVar = apply(targetInt,1,paste, collapse = "-"); targetInt$AllVar = factor(targetInt$AllVar,levels =  expand.grid2.list(list.val))}
        if(length(VarInt)<=1){targetInt$AllVar = target[,VarInt]; targetInt$AllVar = factor(targetInt$AllVar,levels = val)}
        colnames(targetInt) = c(VarInt,"AllVar")
        
        ## Keep only the selected modalities
        ind_kept = which(!is.na(targetInt$AllVar))
        targetInt = targetInt[ind_kept,]
          
        levelsMod = levels(targetInt$AllVar)
        
        ## Create the counts matrix only for the selected subset
        counts_tmp = counts[Taxonomy%in%ind_taxo,]
        counts_tmp = counts_tmp[,colnames(counts_tmp)%in%rownames(targetInt)]
        
        ## Proportions over all the taxonomies
        prop_all = t(counts)/rowSums(t(counts))
        prop_all = as.data.frame(prop_all[,Taxonomy%in%ind_taxo])
        prop_all = as.matrix(prop_all[rownames(prop_all)%in%rownames(targetInt),])
        rownames(prop_all) = targetInt$AllVar
        
        ## Be careful transposition !
        if(aggregate && nrow(counts_tmp)>0 && nrow(targetInt)>0)
        { 
          counts_tmp_combined = aggregate(t(counts_tmp),by=list(targetInt$AllVar),mean)
          rownames(counts_tmp_combined) = counts_tmp_combined$Group.1
          namesCounts = counts_tmp_combined$Group.1
          counts_tmp_combined = as.matrix(counts_tmp_combined[,-1])
        }
        if(!aggregate && nrow(counts_tmp)>0 && nrow(targetInt)>0)
        {  
          counts_tmp_combined = t(counts_tmp)
          prop_tmp_combined = counts_tmp_combined/colSums(counts_tmp)
          rownames(counts_tmp_combined) = targetInt$AllVar
          namesCounts = targetInt$AllVar
          rownames(prop_tmp_combined) = targetInt$AllVar
        }
        
        ## Ordering the counts
        if(!is.null(counts_tmp_combined))
        {
          MeanCounts = apply(counts_tmp_combined,2,mean)
          ord = order(MeanCounts,decreasing=TRUE)
          counts_tmp_combined = as.matrix(counts_tmp_combined[,ord])
          if(!aggregate) prop_tmp_combined = as.matrix(prop_tmp_combined[,ord])
          prop_all = as.matrix(prop_all[,ord])
        }
      }
    }
    
      return(list(counts = counts_tmp_combined,targetInt=targetInt,prop=prop_tmp_combined,namesCounts=namesCounts,levelsMod=levelsMod,prop_all=prop_all))
    
    
  }
  
  
  
  ###########################
  ## Plots for visualisation
  ###########################
  
  Plot_Visu_Barplot <- function(input,resDiff)
  {
    
    ## Get Input for BarPlot
    VarInt = input$VisuVarInt
    ind_taxo = input$selectTaxoPlot
    
    tmp_combined = GetDataToPlot(input,resDiff,VarInt,ind_taxo)
    counts_tmp_combined = tmp_combined$counts
    nbKept = length(ind_taxo)
    SamplesNames = tmp_combined$namesCounts
    
    if(nbKept>1) namesTax = colnames(counts_tmp_combined)
    if(nbKept==1) namesTax = ind_taxo
    
    dataNull = data.frame(x=c(1,2),y=c(1,2))
    plotd3 = nvd3Plot(x ~ y , data = dataNull, type = "multiBarChart", id = 'barplotTaxoNyll',height = input$heightVisu,width=input$widthVisu)
    gg = NULL
    
    if(!is.null(counts_tmp_combined) && nrow(counts_tmp_combined)>0 && length(VarInt)>0)
    { 
      
      ## Create the data frame for the plot function
      dataBarPlot_mat = c()
      tmp_mat = matrix(0,ncol=3,nrow=nbKept)
      tmp_counts = c()
      
        for(i in 1:(nrow(counts_tmp_combined)))
        {
          ## Taxo
          tmp_mat[1:nbKept,1] = namesTax
          
          ## Counts
  
          tmpProp = counts_tmp_combined[i,]
          if(input$CountsOrProp=="prop")
          { 
            tmpProp = round(tmpProp/sum(tmpProp),3)
            tmpProp = as.numeric(tmpProp/sum(tmpProp) * 100)
          }
          tmp_counts = c(tmp_counts,tmpProp)      
          
          ## Meta data
          tmp_mat[1:nbKept,3] = as.character(rep(SamplesNames[i],nbKept))
          
          ## Conbined the sample
          dataBarPlot_mat = rbind(dataBarPlot_mat,tmp_mat)
        }
        
        
        ## Add numeric vector to the dataframe
        dataBarPlot_mat = as.data.frame(dataBarPlot_mat)
        
        colnames(dataBarPlot_mat) = c("Taxonomy","Proportions","AllVar")
        dataBarPlot_mat[,2] = tmp_counts
        if(input$SensPlotVisu == "Vertical") Sens = "multiBarChart"
        if(input$SensPlotVisu == "Horizontal") Sens = "multiBarHorizontalChart"
      
        plotd3 <- nvd3Plot(Proportions ~ AllVar | Taxonomy, data = dataBarPlot_mat, type = Sens, id = 'barplotTaxo',height = input$heightVisu,width=input$widthVisu)
        plotd3$chart(stacked = TRUE)
      
        ##################################
        ## Same plot in ggplot2 for export
        ##################################
      
        tax.colors=rep(c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78", "#2ca02c","#98df8a","#d62728","#ff9896","#9467bd","#c5b0d5","#8c564b",
                         "#c49c94","#e377c2","#f7b6d2","#7f7f7f", "#c7c7c7","#bcbd22","#dbdb8d","#17becf","#9edae5"),ceiling(nbKept/20))
        
        dataBarPlot_mat$Taxonomy = factor(dataBarPlot_mat$Taxonomy,levels = namesTax)
      
        gg= ggplot(dataBarPlot_mat, aes(x=AllVar, y=Proportions, fill=Taxonomy)) 
        gg= gg + geom_bar(stat="identity")
        gg= gg + theme_bw()+ scale_fill_manual(values=tax.colors)
        gg = gg +theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank()) 
        if(input$CountsOrProp=="prop") gg = gg+labs(y="Relative abundance (%)",x="")
        if(input$CountsOrProp=="counts") gg = gg+labs(y="Abundance",x="")
        if(input$SensPlotVisu == "Horizontal") gg = gg + coord_flip()
    } 
    return(list(plotd3=plotd3,gg=gg))
  }
  
  
  
######################################################
##
##            HEATMAP
##
######################################################
  
  
  Plot_Visu_Heatmap <- function(input,resDiff,export=FALSE){
  
  VarInt = input$VisuVarInt
  ind_taxo = input$selectTaxoPlot
  
  counts_tmp_combined = GetDataToPlot(input,resDiff,VarInt,ind_taxo)$counts
  
  if(!is.null(counts_tmp_combined) && nrow(counts_tmp_combined)>0)
  { 
    ## Transform to log2
    counts_tmp_combined = log2(GetDataToPlot(input,resDiff,VarInt,ind_taxo)$counts+1)
   
    col <- switch(input$colors,
                  "green-blue"=colorRampPalette(brewer.pal(9,"GnBu"))(200),
                  "blue-white-red"=colorRampPalette(rev(brewer.pal(9, "RdBu")))(200),
                  "purple-white-orange"=colorRampPalette(rev(brewer.pal(9, "PuOr")))(200),
                  "red-yellow-green"= colorRampPalette(rev(brewer.pal(9,"RdYlGn")))(200))
    
    ## Transpose matrix if Horizontal
    if(input$SensPlotVisu=="Horizontal") counts_tmp_combined = t(as.matrix(counts_tmp_combined))
    
    if(!export) plot = d3heatmap(counts_tmp_combined, dendrogram = "none", Rowv = NA, Colv = NA, na.rm = TRUE,width = input$widthVisu, height = input$heightVisu, show_grid = FALSE, colors = col, scale = input$scaleHeatmap,cexRow = 0.6)
    
    if(export) plot = heatmap.2(counts_tmp_combined, dendrogram = "none", Rowv = NA, Colv = NA, na.rm = TRUE, density.info="none", margins=c(12,8),trace="none",srtCol=45,col = col, scale = input$scaleHeatmap,cexRow = 0.6)
    return(plot)
  }

  
  }

  ######################################################
  ##
  ##            BOXPLOT
  ##
  ######################################################
  
  
  Plot_Visu_Boxplot <- function(input,resDiff,alpha=0.7){
    
    gg = NULL
    
    ## Colors
    colors = rep(c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78", "#2ca02c","#98df8a","#d62728","#ff9896","#9467bd","#c5b0d5","#8c564b",
                   "#c49c94","#e377c2","#f7b6d2","#7f7f7f", "#c7c7c7","#bcbd22","#dbdb8d","#17becf","#9edae5"),ceiling(nrow(resDiff$target)/20))
    
    ## Get Input for BoxPlot
    VarInt = input$VisuVarInt
    ind_taxo = input$selectTaxoPlot
    
    
    tmp_merge = GetDataToPlot(input,resDiff,VarInt,ind_taxo,aggregate=FALSE)
    counts_tmp_combined = tmp_merge$counts
    levelsMod = tmp_merge$levelsMod
    nbKept = length(ind_taxo)
    
    if(!is.null(counts_tmp_combined) && nrow(counts_tmp_combined)>0 && !is.null(levelsMod))
    { 
    
      if(input$typeDataBox == "Relative") counts_tmp_combined = tmp_merge$prop_all
      else counts_tmp_combined = log2(counts_tmp_combined+1)
      if(nbKept==1) colnames(counts_tmp_combined)=ind_taxo

      ## Create the data frame for the plot function
      dataBarPlot_mat = c()
      tmp_mat = matrix(0,ncol=3,nrow=nbKept)
      tmp_counts = c()
      
      for(i in 1:(nrow(counts_tmp_combined)))
      {
        ## Taxo
        tmp_mat[1:nbKept,1] = colnames(counts_tmp_combined)
        
        ## Counts        
        tmpProp = counts_tmp_combined[i,]
        tmp_counts = c(tmp_counts,tmpProp)      
        
        ## Meta data
        tmp_mat[1:nbKept,3] = as.character(rownames(counts_tmp_combined)[i],nbKept)
        
        ## Conbined the sample
        dataBarPlot_mat = rbind(dataBarPlot_mat,tmp_mat)
      }
    
      
      dataBarPlot_mat = as.data.frame(dataBarPlot_mat)
      
      colnames(dataBarPlot_mat) = c("Taxonomy","Value","Samples")
      dataBarPlot_mat[,2] = tmp_counts

      if(is.null(input$BoxColorBy) || length(VarInt)<=1){ dataBarPlot_mat$Colors = dataBarPlot_mat$Samples}
      if(!is.null(input$BoxColorBy) && length(VarInt)>1)
        { 
        tmp = strsplit(as.character(dataBarPlot_mat$Samples),"-")
        ind = which(VarInt%in%input$BoxColorBy)
        dataBarPlot_mat$Colors = rapply(tmp, function(x) paste(x[ind],collapse ="-"), how = "unlist")
      }
      # print(dataBarPlot_mat$Colors)
      dataBarPlot_mat$Samples = factor(dataBarPlot_mat$Samples,levels=levelsMod)
      
      gg = ggplot(dataBarPlot_mat,aes(x=Samples,y=Value,fill=Colors))  + geom_boxplot(alpha=alpha) +theme_bw()
      gg = gg  +theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),panel.background = element_blank(),
                      panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) 
      gg = gg + ylab(paste(input$typeDataBox, "abundance")) +scale_fill_manual(values = colors) + guides(fill=FALSE)
      if(input$CheckAddPointsBox) gg = gg + geom_point(position=position_jitterdodge(dodge.width=0.9))
      if(input$SensPlotVisu=="Horizontal") gg = gg + coord_flip()
      if(nbKept>1) gg = gg + facet_wrap(~ Taxonomy,scales = input$ScaleBoxplot)
    }
    
   return(gg)
    
  }
  
  
  
  ######################################################
  ##
  ##            Scatter plot
  ##
  ######################################################
  
  
  Plot_Visu_Scatterplot<- function(input,resDiff,export=FALSE,lmEst = FALSE,CorEst=FALSE){
    
    plot = NULL
    regCoef = NULL
    Rsq = NULL
    cor.est = NULL
    cor.pvalue = NULL
    dds = resDiff$dds
    counts = as.data.frame(round(counts(dds, normalized = TRUE)))
    target = as.data.frame(resDiff$target)
    data = cbind(target,log2(t(counts)+1))
    
    ## Get Input for ScatterPlot
    Xvar = input$Xscatter
    Yvar = input$Yscatter
    ColBy = input$ColorBy
    PchBy = input$PchBy
    PointSize = input$PointSize
    
    x_var = if (is.null(Xvar)) NULL else data[,Xvar]
    y_var = if (is.null(Yvar)) NULL else data[,Yvar]
    col_var = if (ColBy== "None" || is.null(ColBy)) NULL else data[,ColBy]
    symbol_var = if (PchBy == "None" || is.null(PchBy)) NULL else data[,PchBy]
    size_var = if (PointSize == "None" || is.null(PointSize))  NULL else data[,PointSize]
    
    if(!export && !input$AddRegScatter && !lmEst && !CorEst){
      plot = scatterD3(x = x_var,
                            y = y_var,
                            lab = rownames(data),
                            xlab = Xvar,
                            ylab = Yvar,
                            col_var = col_var,
                            col_lab = ColBy,
                            symbol_var = symbol_var,
                            symbol_lab = PchBy,
                            size_var = size_var,
                            size_lab = PointSize,
                            key_var = rownames(data),
                            height = input$heightVisu,
                            point_opacity = 0.7,
                            labels_size = input$SizeLabelScatter,
                            transitions = TRUE)
      return(plot)
    }

    if(export || input$AddRegScatter){
      if(!lmEst && !CorEst){
        col_var = if (ColBy== "None" || is.null(ColBy)) 1 else data[,ColBy]
        symbol_var = if (PchBy == "None" || is.null(PchBy)) factor(rep(1,nrow(data))) else data[,PchBy]
        size_var = if (PointSize == "None" || is.null(PointSize))  1 else data[,PointSize]
        
        plot = ggplot(data, aes(x = x_var, y = y_var)) + geom_point(aes(color=col_var,size =size_var,shape = symbol_var),alpha=0.7) +theme_bw()
        if(input$SizeLabelScatter!=0) plot = plot + geom_text(aes(label=rownames(data),color=col_var,size=input$SizeLabelScatter/20),vjust = 0,nudge_y =0.05)
        plot = plot + xlab(Xvar) + ylab(Yvar)
        if(input$AddRegScatter) plot = plot + geom_smooth(method="lm")
      return(plot)
      }
    }
    if(lmEst && !CorEst)
    {
      res = lm(y_var~x_var)
      sumRes = summary(res)
      regCoef = sumRes$coefficients 
      rownames(regCoef) = c("Intercept",Xvar)
      Rsq = sumRes$r.squared
      return(list(regCoef=regCoef,Rsq = Rsq))
    }
    if(CorEst)
    {
      print(head(data))
      
      cor.est = cor(as.matrix(data),method = input$CorMeth)
      #cor.pvalue = cor.test(data,method = input$CorMeth)
      return(list(cor.est=cor.est))
    }
  }
  
  
  
  ######################################################
  ##
  ##            GLOBAL VIEW
  ##
  ######################################################
  
  
  Plot_Visu_Diversity <- function(input,resDiff,type="point"){
    
    gg = NULL
    dds = resDiff$dds
    counts = round(counts(dds, normalized = TRUE))
    #target = resDiff$target
    
    ## Get Input for the plot
    VarInt = input$VisuVarInt
    #VarIntBoxDiv = input$VarBoxDiv 
    ind_taxo = rownames(counts)
    
    tmp = GetDataToPlot(input,resDiff,VarInt,ind_taxo,aggregate=FALSE)
    counts_tmp_combined = tmp$counts
    targetInt = tmp$targetInt

    if(nrow(counts_tmp_combined)>0 && !is.null(counts_tmp_combined) && !is.null(targetInt))
    { 
      print(TaxoNumber(counts_tmp_combined))
      alpha <- tapply(TaxoNumber(counts_tmp_combined), targetInt$AllVar, mean)
      gamma <- TaxoNumber(counts_tmp_combined, targetInt$AllVar)
      beta = gamma/alpha - 1
      nb = length(alpha)
#       dataTmp = data.frame(value=c(alpha,beta,gamma),
#                            diversity = c(rep("Alpha",nb),rep("Beta",nb),rep("Gamma",nb)),
#                            Var = as.character(rep(names(alpha),3)),
#                            X = as.character(rep(targetInt[,VarIntBoxDiv],3)))
      dataTmp = data.frame(value=c(alpha,beta,gamma),
                     diversity = c(rep("Alpha",nb),rep("Beta",nb),rep("Gamma",nb)),
                     Var = as.character(rep(names(alpha),3)))
     
      ## Merge targetInt et dataTmp par rapport à Var
#       VectX = c()
#       for(i in 1:nb)
#       {
#         ## If duplicated, take only one row
#         tmpX = which(targetInt$AllVar%in%names(alpha)[i])[1]
#         VectX = c(VectX,targetInt[tmpX,VarIntBoxDiv]) 
#       }
#       print(VectX)
#       dataTmp$X =  as.character(rep(VectX,3))
                               
      dataTmp = dataTmp[dataTmp$diversity%in%input$WhichDiv,]

      if(type=="point")
      { 
        gg = ggplot(dataTmp, aes(x=Var, y=value, color=diversity)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        gg = gg + geom_point(size=input$sizePointGlobal) 
        if(input$SensPlotVisu=="Horizontal") gg = gg + coord_flip()
        if(input$SplitVisuGlobal==TRUE) gg = gg + facet_wrap(~ diversity)
      }
#       if(type=="box")
#       { 
#         gg = ggplot(dataTmp,aes(x=X,y=value,fill=diversity))  + geom_boxplot(alpha=0.7) + theme_bw()  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#         gg = gg + geom_point(size=input$sizePointGlobal) 
#         gg = gg + geom_point(position=position_jitterdodge(dodge.width=0.9))
#         if(input$SensPlotVisuGlobal=="Horizontal") gg = gg + coord_flip()
#         if(input$SplitVisuGlobal==TRUE) gg = gg + facet_wrap(~ diversity) 
#       }
      
#       nvd3Plot(value ~ Var | diversity, data = dataTmp, id = 'Scachart', type = 'lineChart',height = 1000,width=1000)
#       p1$xAxis(axisLabel = 'Variable of interest')
    }
    return(gg)
    
  }

  
  
  ## Rarefaction
  Plot_Visu_Rarefaction <- function(input,resDiff,xlim,ylim,ylab="Species"){
    
    PlotRare = NULL
    dds = resDiff$dds
    
    ## Taxo in columns
    counts = t(round(counts(dds, normalized = TRUE)))
    
    if(nrow(counts)>0 && !is.null(counts))
    { 
      max <- max(rowSums(counts))
      raremax <- min(rowSums(counts))
      #PlotRare = rarefaction_curve(counts, step = 10, taxo = "Taxonomy level") 
      options(warn=-1)
      PlotRare = rarecurve(counts, step = max(1,ceiling(max/60)),sample=raremax, col = "blue", cex = 0.9,xlim=xlim,ylim=ylim, ylab=ylab) 
      options(warn=0)
    }
    return(PlotRare)
  }
  
  
  
  ## Get the non-zero taxo by sample  
  TaxoNumber <-  function (x, groups, mar = 1) 
  {
    if (!missing(groups)) 
    {
      if (length(groups) == 1) groups = rep(groups, nrow(x))
      x = aggregate(x, list(groups), max)
      rownames(x) = x[, 1]; x = x[, -1]
    }
    if (length(dim(x)) > 1) res = apply(x > 0, mar, sum)
    else res = sum(x > 0)
  }
  
  
  rarefaction_curve <- function (x, step = 1, taxo ="Species") 
  {
    
    tot = rowSums(x)
    S = TaxoNumber(x)
    if (any(S <= 0)) {
      x <- x[S > 0, , drop = FALSE]
      tot <- tot[S > 0]
      S <- S[S > 0]
    }
    nr <- nrow(x)

    out <- lapply(seq_len(nr), function(i) {
      n <- seq(1, tot[i], by = step)
      if (n[length(n)] != tot[i]) n <- c(n, tot[i])
      drop(rarefy(x[i, ], n))
    })
    
    
    df = data.frame()
    
    for(i in 1: length(out))
    {
      dftmp = data.frame(x=attr(out[[i]], "Subsample"),y=out[[i]],samples=rep(rownames(x)[i],length(out[[i]])))
      df = rbind(df,dftmp)
    }
    
    Nmax = sapply(out, function(x) max(attr(x, "Subsample")))
    Smax = sapply(out, max)
    
#     plot =  nvd3Plot(y ~ x | samples, data = df, id = 'chart', type = 'lineChart',height=600)
#     plot$xAxis(axisLabel = 'Sample size')
    
    plot =  ggplot(df,aes=c(x=x,y=y,  group=samples, colour=samples)) + geom_line()+xlab('Sample size') 
    plot = plot + theme_bw() + theme(legend.position="bottom")

    return(plot)
  }
  
  
  TableDiff_print <- function(input,BaseContrast,resDiff, info = NULL) 
  {
    VarInt = input$VarInt
    dds = resDiff$dds
    counts = resDiff$counts
    target = resDiff$target
    group = ""
    if(length(VarInt)>1) for(i in VarInt){ group = paste(group,target[,i], sep="-") }
    else group = target[,VarInt]
    conds = unique(group)
    
    result = list()
    alpha = input$AlphaVal
    cooksCutoff = ifelse(input$CooksCutOff!='Auto',ifelse(input$CooksCutOff!=Inf,input$CutOffVal,Inf),TRUE)
    result[[input$ContrastList_table]] <- results(dds,contrast=BaseContrast[,input$ContrastList_table],pAdjustMethod=input$AdjMeth,
                                                  cooksCutoff=cooksCutoff,
                                                  independentFiltering=input$IndFiltering,alpha=alpha)
    
    #names(result) <- gsub("_", " ", names(result))
    if (is.null(info))  info <- data.frame(Id = rownames(result[[1]]))
    else names(info)[1] <- "Id"
    if (any(duplicated(info[, 1]))) stop("Duplicated IDs in the annotations")
    counts <- data.frame(Id = rownames(counts(dds)), counts(dds), round(counts(dds, normalized = TRUE)))
    colnames(counts) <- c("Id", colnames(counts(dds)), paste0("norm.", colnames(counts(dds))))
    bm <- data.frame(Id = rownames(result[[1]]), baseMean = round(result[[1]][,"baseMean"], 2))
    
    base <- merge(merge(info, counts, by = "Id", all.y = TRUE), bm, by = "Id")
    tmp <- base[, paste("norm", colnames(counts(dds)), sep = ".")]
    for (cond in conds) {
      base[, cond] <- round(apply(as.data.frame(tmp[, group == cond]), 1, mean), 0)
    }
    complete.complete <- base
    complete <- vector("list", length(result))
    names(complete) <- names(result)
    for (name in names(result)) {
      complete.name <- base
      
      ### ??????????
      conds.supp <- setdiff(conds, gsub("\\(|\\)", "", unlist(strsplit(name, " vs "))))
      if (length(conds.supp) > 0) {
        complete.name <- complete.name[, -which(names(complete.name) %in% conds.supp)]
        samples.supp <- colnames(counts(dds))[group %in% conds.supp]
        col.supp <- c(samples.supp, paste("norm", samples.supp, sep = "."))
        complete.name <- complete.name[, -which(names(complete.name) %in% col.supp)]
      }
      ### ??????????

      res.name <- data.frame(Id = rownames(result[[name]]), 
                             FC = round(2^(result[[name]][, "log2FoldChange"]), 3),
                             log2FoldChange = round(result[[name]][, "log2FoldChange"], 3),
                             pvalue = result[[name]][, "pvalue"], 
                             padj = result[[name]][, "padj"])
      complete.name <- merge(complete.name, res.name, by = "Id")
      mcols.add <- data.frame(Id = rownames(counts(dds)), 
                              dispGeneEst = mcols(dds)$dispGeneEst, dispFit = mcols(dds)$dispFit, 
                              dispMAP = mcols(dds)$dispMAP, dispersion = mcols(dds)$dispersion, 
                              betaConv = mcols(dds)$betaConv, maxCooks = mcols(dds)$maxCooks)
      if (is.null(cooksCutoff)){
        m <- nrow(attr(dds, "modelMatrix"))
        p <- ncol(attr(dds, "modelMatrix"))
        cooksCutoff <- qf(0.99, p, m - p)
      }
      mcols.add$outlier <- ifelse(mcols(dds)$maxCooks > cooksCutoff, "Yes", "No")
      complete.name <- merge(complete.name, mcols.add, by = "Id")
      complete[[name]] <- complete.name
      up.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange>=0), ]
      up.name <- up.name[order(up.name$padj), ]
      down.name <- complete.name[which(complete.name$padj<=alpha & complete.name$betaConv & complete.name$log2FoldChange<=0), ]
      down.name <- down.name[order(down.name$padj), ]

      name <- gsub(" ", "", name)
      keep <- c("Id","baseMean","FC","log2FoldChange","padj")
#       complete.complete[, paste(name, keep, sep = ".")] <- complete.name[, keep]
    }
    #return(list(complete=complete.name,up=up.name,down=down.name))
    return(list(complete=complete.name[,keep],up=up.name[,keep],down=down.name[,keep]))
  }
  
  
  Get_log2FC_padj <-function(input,BaseContrast,resDiff, info = NULL)
  {
    log2FC = NULL
    padj = NULL
    VarInt = input$VarInt
    dds = resDiff$dds
    counts = resDiff$counts
    target = resDiff$target
    SelContrast = colnames(BaseContrast)
    nbCont = length(SelContrast)
    result = list()
    alpha = input$AlphaVal
    cooksCutoff = ifelse(input$CooksCutOff!='Auto',ifelse(input$CooksCutOff!=Inf,input$CutOffVal,Inf),TRUE)

      for(i in 1:nbCont)
      { 
        cont = as.character(SelContrast[i])
        result[[cont]] <- results(dds,contrast=BaseContrast[,cont],pAdjustMethod=input$AdjMeth,
                                  cooksCutoff=cooksCutoff,
                                  independentFiltering=input$IndFiltering,alpha=alpha)
      }
      log2FC = as.matrix(round(result[[SelContrast[1]]][, "log2FoldChange"], 3))
      padj = as.matrix(round(result[[SelContrast[1]]][, "padj"], 3))
      if(nbCont>1)
      {
        for(i in 2:nbCont)
        {
          log2FC = cbind(log2FC,round(result[[SelContrast[i]]][, "log2FoldChange"], 3))
          padj = cbind(padj,round(result[[SelContrast[i]]][, "padj"], 7))
        }
        colnames(log2FC) = names(result)
        colnames(padj) = names(result)
      }
      rownames(log2FC) = rownames(result[[SelContrast[1]]])
      rownames(padj) = rownames(result[[SelContrast[1]]])

    return(list(log2FC=as.data.frame(log2FC),padj=padj))
  }
  
  
  Plot_Visu_Heatmap_FC <- function(input,BaseContrast,resDiff,export=FALSE){
    
    res = NULL
    SelContrast = input$ContrastList_table_FC
    log2FC = Get_log2FC_padj(input,BaseContrast,resDiff, info = NULL)$log2FC

    
    if(!is.null(log2FC) && length(SelContrast)>=2)
    { 
      cont = which(colnames(log2FC)%in%SelContrast)
      log2FC = log2FC[,SelContrast] 
      ind_taxo = input$selectTaxoPlot
      ind = rownames(log2FC)%in%ind_taxo
      log2FC = log2FC[ind,]
      
      
      col <- c(colorRampPalette(c("royalblue4","royalblue3","royalblue2","royalblue1","white"))(n = 100),colorRampPalette(c("white",  "firebrick1", "firebrick2", "firebrick3", "firebrick4"))(n = 100))
      ## Transpose matrix if Horizontal
      if(input$SensPlotVisu=="Horizontal") log2FC = t(as.matrix(log2FC))
      
      if(!export) res = d3heatmap(log2FC, dendrogram = "row", Rowv = TRUE, Colv = NA, na.rm = TRUE, width = input$widthVisu, height = input$heightVisu, show_grid = FALSE, colors = col, scale = input$scaleHeatmap,cexRow = input$LabelSizeHeatmap,cexCol =input$LabelSizeHeatmap, offsetCol=input$LabelColOffsetHeatmap,offsetRow=input$LabelRowOffsetHeatmap)
    
      
      if(export) res = heatmap.2(log2FC, dendrogram = "row", Rowv = TRUE, Colv = NA, na.rm = TRUE, width = input$widthVisu, height = input$heightVisu, show_grid = FALSE, colors = col, scale = input$scaleHeatmap,cexRow = input$LabelSizeHeatmap,cexCol =input$LabelSizeHeatmap, offsetCol=input$LabelColOffsetHeatmap,offsetRow=input$LabelRowOffsetHeatmap)

    }
    return(res)
  }


  