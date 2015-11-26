


  GetDataFromBIOM <-function(dataBIOM)
  {
    
    counts = biom_data(dataBIOM)
    counts = as.matrix(counts)
    counts = as.data.frame(counts)
    taxo = as.data.frame(observation_metadata(dataBIOM))

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
      out = paste(out,paste("<b>",contnames[i], ":</b> <br/>", paste(contrast[contrast != 0], coefs[contrast != 0], collapse = " ", sep = " ")),"<br/>")
    }
    return(out)
    
  }

  
  
  ## Get the counts for the selected taxonomy
  GetCountsMerge <- function(input,dataInput,taxoSelect,target,design)
  {
    counts= NULL
    CheckTarget = FALSE
    
    ## Counts and taxo tables
    CT = dataInput$counts
    taxo = dataInput$taxo
        
    ## Select cols in the target
    labels = target[,1]
    ind = which(colnames(CT)%in%labels)
    
    
    if(length(ind)==length(labels))
    { 
      CT = CT[,ind]
      
      ## Order CT according to the target
      CT = OrderCounts(CT,labels)
#       ind0 = which(rowSums(CT)==0)
#       if(length(ind0)>0) CT = CT[-ind0,]
      
      ## Counts normalisation
      dds <- DESeqDataSetFromMatrix(countData=CT, colData=target, design=design)
      dds <- estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)))

      CT = as.data.frame(round(counts(dds, normalized = TRUE)))
      ordOTU = order(rownames(taxo))
      indOTU_annot = which(rownames(CT)%in%rownames(taxo))
      counts_annot = CT[indOTU_annot[ordOTU],]
      
      if(taxoSelect=="OTU") counts = counts_annot
      else{
      taxoS = taxo[ordOTU,taxoSelect]
      counts = aggregate(counts_annot,by=list(Taxonomy = taxoS),sum)
      rownames(counts)=counts[,1];counts=counts[,-1]
      }
      
      ## Ordering the counts table according to the target labels 
      counts = OrderCounts(counts,labels)
      CheckTarget = TRUE
    }
    return(list(counts=counts,CheckTarget=CheckTarget))
  }

  ## Order the counts 
  OrderCounts <- function(counts,labels)
  {
    n = length(labels)
    CountsOrder = counts

    for(i in 1:n)
    {
      
      ind = which(labels[i]==colnames(counts))
      CountsOrder[,i] = counts[,ind]
    }
    colnames(CountsOrder) = labels
    return(CountsOrder)
  }
  
  
  ## Get the dds object of DESeq2
  Get_dds_object <- function(input,counts,target,design)
  {
    
    dds <- DESeqDataSetFromMatrix(countData=counts, colData=target, design=design)
    normFactors = rep(1,nrow(target))
    ## Size factor estimation
    #dds <- estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)))
    #normalizationFactors(dds) <- normFactors
    sizeFactors(dds)<- normFactors
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
    #colors = c("dodgerblue","firebrick1","MediumVioletRed","SpringGreen")
    colors = c("SpringGreen","dodgerblue","black","firebrick1")
    VarInt = input$VarInt
    dds = resDiff$dds
    counts = resDiff$counts
    target = resDiff$target
    group = as.data.frame(target[,VarInt])
    rownames(group) = rownames(target)
    
    ## If more than 4 levels for one factor
    if(length(VarInt)>1)  maxFact =max(sapply(group,FUN=function(x) length(levels(x))))
    else maxFact = length(levels(group))
    if(maxFact>=4) colors = rainbow(maxFact) 
    
    if(input$DiagPlot=="barplotTot") barplotTot(input,counts,group = group, col=colors)
    if(input$DiagPlot=="barplotNul") barPlotNul(input,counts, group = group, col=colors)
    if(input$DiagPlot=="densityPlot") densityPlotTot(input,counts, group = group, col=colors)
    if(input$DiagPlot=="MajTax") majTaxPlot(input,counts, group = group, col=colors)
    if(input$DiagPlot=="SERE") SEREplot(input,counts)
    if(input$DiagPlot=="Sfactors") diagSFactors(input,dds,frame=1) 
    if(input$DiagPlot=="SfactorsVStot") diagSFactors(input,dds,frame=2) 
    if(input$DiagPlot=="pcaPlot") PCAPlot_meta(input,dds, group,  type.trans = input$TransType, col = colors)
    if(input$DiagPlot=="pcoaPlot") PCoAPlot_meta(input,dds, group) 
    if(input$DiagPlot=="clustPlot") HCPlot(input,dds,group,type.trans=input$TransType)
  }

  
#   HCPlot <- function (input,dds,group,type.trans,col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen")) 
#   {
#     counts = as.data.frame(round(counts(dds, normalized = TRUE)))
#     if (type.trans == "VST") counts.trans <- assay(varianceStabilizingTransformation(dds))
#     if (type.trans == "rlog") counts.trans <- assay(rlogTransformation(dds))
#     
#     hc <- hclust(dist(t(counts.trans)), method = "ward.D")
#     
#     type <- switch(input$typeHculst,
#                   "radial"="radial",
#                   "fan"="fan",
#                   "triangle"="cladogram",,
#                   "hori"= "hori",
#                   "verti"=NULL)
#     
#     par(cex=input$cexLabelDiag,mar=c(12,5,8,5))
#     if(input$colorHC && type=="hori") 
#     {
#       hc = dendrapply(as.dendrogram(hc),colLabdendo,group) 
#       plot(hc, xlab = "Euclidean distance, Ward criterion", main = "Cluster dendrogram")
#     }
#     
#     if(!input$colorHC && type=="hori") 
#     {
#       plot(hc, xlab = "Euclidean distance, Ward criterion", main = "Cluster dendrogram",hang=-1)
#     }
#     
#     if(type!="hori") 
#     { 
#       group = apply(group,1,paste, collapse = "-")
#       nb = length(unique(group))
#       plot(as.phylo(hc), type= type,label.offset = 1, tip.color = ifelse(input$colorHC, rainbow(nb)[as.integer(as.factor(group))], rep(1,nb)))
#     }
#     dev.off() 
#   }
  
  HCPlot <- function (input,dds,group,type.trans,col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen")) 
  {
    
    ## Get the counts
    counts = as.data.frame(round(counts(dds, normalized = TRUE)))
    if (type.trans == "VST") counts.trans <- assay(varianceStabilizingTransformation(dds))
    if (type.trans == "rlog") counts.trans <- assay(rlogTransformation(dds))
    
    ## Get the group of leaf
    group = apply(group,1,paste, collapse = "-")    
    nb = length(unique((group)))
    
    ## Get the dendrogram
    hc <- hclust(dist(t(counts.trans)), method = "ward.D")
    dend = as.dendrogram(hc)
    
    ## Get the type of dendrogram
    type <- switch(input$typeHculst,
                   "fan"="fan",
                   "hori"= "hori")
    
    dend <- set(dend, "labels_cex", input$cexLabelDiag)
    if(input$colorHC) labels_colors(dend)<-rainbow(nb)[as.integer(as.factor(group))][order.dendrogram(dend)]
    
    if(type=="hori") 
    { 
      par(mar = c(8,4,4,2))
      plot(dend, main = "Cluster dendrogram")
    }  
    if(type!="hori")
    {
      par(mar = c(0.3,2,0.3,2))
      circlize_dendrogram(dend, labels_track_height = 0.2, dend_track_height = .3, main = "Cluster dendrogram")
    }
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
    par(cex=input$cexLabelDiag,mar=c(12,5,4,5))
    barplot(colSums(counts), cex.names = cex.names, main = "Total mapped read count per sample", ylab = "Total mapped read count", 
            ylim = c(0, max(colSums(counts)) * 1.2), density = if (ncol1) {NULL}
            else {15}, 
            angle = if (ncol1) {NULL}
            else {c(-45, 0, 45, 90)[as.integer(group[, 2])]}, col = col[as.integer(group[, 1])], las = 2)
    legend("topright", levels(group[, 1]), fill = col[1:nlevels(group[,1])], bty = "n")
    if (!ncol1)  legend("topleft", levels(group[, 2]), density = 15,col = 1, angle = c(-45, 0, 45, 90)[1:nlevels(group[, 2])], bty = "n")
  
  }


  ## barplot Nul 
  barPlotNul <-function (input,counts, group, cex.names = 1, col = c("lightblue","orange", "MediumVioletRed", "SpringGreen")) 
  {
    
    percentage <- apply(counts, 2, function(x) {sum(x == 0)}) * 100/nrow(counts)
    percentage.allNull <- (nrow(counts) - nrow(removeNulCounts(counts))) * 100/nrow(counts)
    ncol1 <- ncol(group) == 1
    
    par(cex=input$cexLabelDiag,mar=c(12,5,4,5))

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
  densityPlotTot <-function (input,counts, group, col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen")) 
  {
    
    counts <- removeNulCounts(counts)
    ncol1 <- ncol(group) == 1
    par(cex=input$cexLabelDiag,mar=c(8,5,4,5))
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

    x <- barplot(maj, col = col[as.integer(group[, 1])], main = "Proportion of mapped reads from\nmost expressed sequence",
                 ylim = c(0, max(maj) * 1.2), cex.main = 1, 
                 cex.names = cex.names, las = 2, ylab = "Proportion of mapped reads", 
                 density = if (ncol1) {NULL}
                 else {15}, 
                 angle = if (ncol1) {NULL}
                 else {c(-45, 0, 45, 90)[as.integer(group[, 2])]})
    
    legend("topright", levels(group[, 1]), fill = col[1:nlevels(group[,1])], bty = "n")
    if (!ncol1) legend("topleft", levels(group[, 2]), density = 15, col = 1, 
                       angle = c(-45, 0, 45, 90)[1:nlevels(group[, 2])], bty = "n")
    
    for (i in 1:length(seqname)) text(x[i], maj[i]/2, seqname[i], cex=input$cexLabelDiag, srt = 90, adj = 0)
  }
  

  ## plot SERE Coefs
  SEREplot<-function(input,counts) 
  {
    sere = SEREcoef(counts)
    print(sere)
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
    print("OK1")
    laneTotals <- colSums(observed)
    total <- sum(laneTotals)
    print("OK2")
    fullObserved <- observed[rowSums(observed) > 0, ]
    fullLambda <- rowSums(fullObserved)/total
    fullLhat <- fullLambda > 0
    print("OK3")
    fullExpected <- outer(fullLambda, laneTotals)
    fullKeep <- which(fullExpected > 0)
    print(fullKeep)
    print(fullExpected)
    oeFull <- (fullObserved[fullKeep] - fullExpected[fullKeep])^2/fullExpected[fullKeep]
    print(oeFull)
    dfFull <- length(fullKeep) - sum(fullLhat != 0)
    print(dfFull)
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

  
  ### PCoA
  PCoAPlot_meta <-function (input,dds, group_init,col = c("SpringGreen","dodgerblue","black","firebrick1"), plot = "pcoa") 
  {
    cval=c()
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
    
    ## Get the norm data
    counts.norm = as.data.frame(round(counts(dds, normalized = TRUE)))
    
    # was removed
    counts.norm = counts.norm[,ind_kept]
    
    ## Get the distance
    dist.counts.norm = vegdist(t(counts.norm), method = input$DistPCOA)
    
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
    print(condition_set)
    print(time_set)
    if (length(time_set) == 1 && length(condition_set) <= 4){
      cval = apply(expand.grid(condition_set,time_set),1,paste, collapse = "-")
      cval = sort(cval)
    }
    print(col)
    # to reactivate
    #pco.counts.norm$li = pco.counts.norm$li[ind_kept,]
    if (plot == "pcoa"){
      ## Plot axis, label and circles
      plot(pco.counts.norm$li[1:2], xlab=paste("PC1 : ",round(eigen[1],1),"%") , ylab=paste("PC2 : ",round(eigen[2],1),"%"),
           xlim=c(min+0.25*min,max+0.25*max), ylim=c(min-0.1,max+0.1), cex.axis=1, cex.lab=1,lwd=2, type="n")
      # Set different shapes
      if(input$labelPCOA == "Group"){
        print(cval)
        print(length(cval))
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
  
  ### PCA
  PCAPlot_meta <-function (input,dds, group, n = min(500, nrow(counts(dds))), type.trans = c("VST", "rlog"), 
                           col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen"),plot="pca") 
  {
    type.trans <- type.trans[1]
    
    if (type.trans == "VST") counts.trans <- assay(varianceStabilizingTransformation(dds))
    else counts.trans <- assay(rlogTransformation(dds))
    
    rv = apply(counts.trans, 1, var, na.rm = TRUE)
    pca = prcomp(t(counts.trans[order(rv, decreasing = TRUE),][1:n, ]))
    
   
    
    
    
    if(plot=="pca")
    { 
      prp <- pca$sdev^2 * 100/sum(pca$sdev^2)
      prp <- round(prp, 2)
      ncol1 <- ncol(group) == 1
      
      par(mfrow = c(1, 2))
      
      abs = range(pca$x[, 1])
      abs = abs(abs[2] - abs[1])/25
      ord = range(pca$x[, 2])
      ord = abs(ord[2] - ord[1])/25
      
      par(mar=c(8,5,4,5))
      plot(pca$x[, 1], pca$x[, 2], las = 1, cex = 2, col = col[as.integer(group[,1])], 
           pch = if (ncol1) {16}
           else {c(16:18, 25)[as.integer(group[, 2])]},
           xlab = paste0("PC1 (", prp[1], "%)"),
           ylab = paste0("PC2 (", prp[2], "%)"), 
           main = "Principal Component Analysis",
            )
      abline(h = 0, v = 0, lty = 2, col = "lightgray")
      text(pca$x[, 1] - ifelse(pca$x[, 1] > 0, abs, -abs), pca$x[,2] - ifelse(pca$x[, 2] > 0, ord, -ord), colnames(counts.trans), col = col[as.integer(group[, 1])])
      abs = range(pca$x[, 1])
      abs = abs(abs[2] - abs[1])/25
      ord = range(pca$x[, 3])
      ord = abs(ord[2] - ord[1])/25
      plot(pca$x[, 1], pca$x[, 3], las = 1, cex = 2, col = col[as.integer(group[, 1])], 
           pch = if (ncol1) {16}
           else {c(16:18, 25)[as.integer(group[, 2])]}, 
           xlab = paste0("PC1 (", prp[1], "%)"), 
           ylab = paste0("PC3 (", prp[3], "%)"), 
           main = "Principal Component Analysis")
      abline(h = 0, v = 0, lty = 2, col = "lightgray")
      text(pca$x[, 1] - ifelse(pca$x[, 1] > 0, abs, -abs), pca$x[,3] - ifelse(pca$x[, 3] > 0, ord, -ord), colnames(counts.trans), col = col[as.integer(group[, 1])],cex=input$cexLabelDiag)
    }
    
    if(plot=="eigen") barplot(pca$sdev^2, main = "Eigen values of the PCA", names.arg = 1:length(pca$sdev), xlab = "Axes")
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
  
  GetDataToPlot <- function(resDiff,VarInt,ind_taxo,aggregate=TRUE)
  {
    dds = resDiff$dds
    counts = as.data.frame(round(counts(dds, normalized = TRUE)))
    target = resDiff$target
    counts_tmp_combined = NULL
    prop_tmp_combined = NULL
    targetInt = NULL
    ## Select a subset within the taxonomy level (default is the 12 most abundant)
    nbKept = length(ind_taxo)
    Taxonomy = rownames(counts)
    
    if (length(VarInt)>0 && nbKept>0)
    { 
      ## Create the variable to plot
      targetInt = as.data.frame(target[,VarInt])
      rownames(targetInt)=target[,1]  
      if(length(VarInt)>1) targetInt$AllVar = apply(targetInt,1,paste, collapse = "-")
      if(length(VarInt)<=1)  targetInt$AllVar = target[,VarInt]
      colnames(targetInt) = c(VarInt,"AllVar")
      ## Create the counts matrix only for the selected subset
      counts_tmp = counts[Taxonomy%in%ind_taxo,]

      ## Be careful transposition !
      if(aggregate)
      { 
        counts_tmp_combined = aggregate(t(counts_tmp),by=list(targetInt$AllVar),sum)
        rownames(counts_tmp_combined) = counts_tmp_combined$Group.1
        counts_tmp_combined = as.matrix(counts_tmp_combined[,-1])
      }
      if(!aggregate)
      {  
        counts_tmp_combined = t(counts_tmp)
        prop_tmp_combined = counts_tmp_combined/colSums(counts)
        rownames(counts_tmp_combined) = targetInt$AllVar
        rownames(prop_tmp_combined) = targetInt$AllVar
      }
      ## Ordering the counts
      MeanCounts = apply(counts_tmp_combined,2,mean)
      ord = order(MeanCounts,decreasing=TRUE)
      counts_tmp_combined = as.matrix(counts_tmp_combined[,ord])
      if(!aggregate) prop_tmp_combined = as.matrix(prop_tmp_combined[,ord])
    }
    
      return(list(counts = counts_tmp_combined,targetInt=targetInt,prop=prop_tmp_combined))
    
    
  }
  
  
  
  ###########################
  ## Plots for visualisation
  ###########################
  
  Plot_Visu_Barplot <- function(input,resDiff)
  {

    ## Get Input for BarPlot
    VarInt = input$VisuVarIntBP
    ind_taxo = input$selectTaxoPlotBP
    
    counts_tmp_combined = GetDataToPlot(resDiff,VarInt,ind_taxo)$counts
    nbKept = length(ind_taxo)
    
    if(!is.null(counts_tmp_combined) && nrow(counts_tmp_combined)>0)
    { 
      counts_tmp_combined = GetDataToPlot(resDiff,VarInt,ind_taxo)$counts
      Taxonomy = rownames(counts_tmp_combined)
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
          if(input$CountsOrProp=="prop")
          { 
            tmpProp = round(tmpProp/sum(tmpProp),3)
            tmpProp = as.numeric(tmpProp/sum(tmpProp) * 100)
          }
          tmp_counts = c(tmp_counts,tmpProp)      
          
          ## Meta data
          tmp_mat[1:nbKept,3] = as.character(rep(rownames(counts_tmp_combined)[i],nbKept))
          
          ## Conbined the sample
          dataBarPlot_mat = rbind(dataBarPlot_mat,tmp_mat)
        }
        
        
        ## Add numeric vector to the dataframe
        dataBarPlot_mat = as.data.frame(dataBarPlot_mat)
        
        colnames(dataBarPlot_mat) = c("Taxonomy","Proportions","AllVar")
        dataBarPlot_mat[,2] = tmp_counts
  
        plotd3 <- nvd3Plot(Proportions ~ AllVar | Taxonomy, data = dataBarPlot_mat, type = input$SensPlotVisuBP, id = 'barplotTaxo',height = input$heightVisu,width=input$widthVisu)
        plotd3$chart(stacked = TRUE)
    } 
    else{ 
      ## Pb affichage quand data NULL
      dataNull = data.frame(x=c(1,2),y=c(1,2))
      plotd3 = nvd3Plot(x ~ y , data = dataNull, type = "multiBarChart", id = 'barplotTaxoNyll',height = input$heightVisu,width=input$widthVisu)
    }
    return(plotd3)
  }
  
  
  
######################################################
##
##            HEATMAP
##
######################################################
  
  
  Plot_Visu_Heatmap <- function(input,resDiff){
  
  VarInt = input$VisuVarIntHM
  ind_taxo = input$selectTaxoPlotHM
  
  counts_tmp_combined = GetDataToPlot(resDiff,VarInt,ind_taxo)$counts
  
  if(!is.null(counts_tmp_combined) && nrow(counts_tmp_combined)>0)
  { 
    ## Transform to log2
    counts_tmp_combined = log2(GetDataToPlot(resDiff,VarInt,ind_taxo)$counts+1)
   
    col <- switch(input$colors,
                  "green-blue"=colorRampPalette(brewer.pal(9,"GnBu"))(200),
                  "blue-white-red"=colorRampPalette(rev(brewer.pal(9, "RdBu")))(200),
                  "purple-white-orange"=colorRampPalette(rev(brewer.pal(9, "PuOr")))(200),
                  "red-yellow-green"= colorRampPalette(rev(brewer.pal(9,"RdYlGn")))(200))
    
    ## Transpose matrix if Horizontal
    if(input$SensPlotVisuHM=="Horizontal") counts_tmp_combined = t(as.matrix(counts_tmp_combined))
         #print(counts_tmp_combined)
    return(heatmap.2(counts_tmp_combined, dendrogram = "none", Rowv = NA, Colv = NA, na.rm = TRUE, density.info="none", margins=c(12,8),trace="none",srtCol=45,
                    col = col, scale = input$scaleHeatmap,cexRow = 0.6))
#     return(d3heatmap(counts_tmp_combined, dendrogram = "none", Rowv = NA, Colv = NA, na.rm = TRUE, 
#                      width = 1500, height = 1000, show_grid = FALSE, colors = col, scale = input$scaleHeatmap,
#                      cexRow = 0.6))
  }

  
  }

  ######################################################
  ##
  ##            BOXPLOT
  ##
  ######################################################
  
  
  Plot_Visu_Boxplot <- function(input,resDiff){
    
    gg = NULL
    ## Get Input for BoxPlot
    VarInt = input$VisuVarIntBoxP
    ind_taxo = input$selectTaxoPlotBoxP
    
    tmp_merge = GetDataToPlot(resDiff,VarInt,ind_taxo,aggregate=FALSE)
    counts_tmp_combined = tmp_merge$counts

    nbKept = length(ind_taxo)
    
    if(!is.null(counts_tmp_combined) && nrow(counts_tmp_combined)>0)
    { 
      Taxonomy = rownames(counts_tmp_combined)
    
      if(input$typeDataBox == "Relative") 
      { 
        counts_tmp_combined = tmp_merge$prop
      }
      if(input$typeDataBox == "Log2") counts_tmp_combined = log2(counts_tmp_combined+1)
      
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
      
      gg = ggplot(dataBarPlot_mat,aes(x=Samples,y=Value,fill=Samples))  + geom_boxplot(alpha=0.7) + theme_bw()  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
      gg = gg + ylab(input$typeDataBox)
      if(input$CheckAddPointsBox) gg = gg + geom_point(position=position_jitterdodge(dodge.width=0.9))
      if(input$SensPlotVisuBoxP=="Horizontal") gg = gg + coord_flip()
      if(nbKept>1) gg = gg + facet_wrap(~ Taxonomy)
    }
    
   return(gg)
    
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
    VarInt = input$VisuVarIntDiv
    VarIntBoxDiv = input$VarBoxDiv 
    ind_taxo = rownames(counts)
    
    tmp = GetDataToPlot(resDiff,VarInt,ind_taxo,aggregate=FALSE)
    counts_tmp_combined = tmp$counts
    targetInt = tmp$targetInt

    if(nrow(counts_tmp_combined)>0 && !is.null(counts_tmp_combined) && !is.null(targetInt))
    { 
      alpha <- tapply(TaxoNumber(counts_tmp_combined), targetInt$AllVar, mean)
      gamma <- TaxoNumber(counts_tmp_combined, targetInt$AllVar)
      beta = gamma/alpha - 1
      nb = length(alpha)
      dataTmp = data.frame(value=c(alpha,beta,gamma),
                           diversity = c(rep("Alpha",nb),rep("Beta",nb),rep("Gamma",nb)),
                           Var = as.character(rep(names(alpha),3)), X = as.character(rep(targetInt[,VarIntBoxDiv],3)))
     
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
        if(input$SensPlotVisuGlobal=="Horizontal") gg = gg + coord_flip()
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
  
  
  
  
  
  #   ## Get tables from contrasts
#   GetTableFromContrast  <- function (object, contrast, name, lfcThreshold = 0, altHypothesis = c("greaterAbs","lessAbs", "greater", "less"), listValues = c(1, -1), cooksCutoff, 
#                             independentFiltering = TRUE, alpha = 0.1, filter, theta, pAdjustMethod = "BH", format = c("DataFrame", "GRanges", "GRangesList"), test, addMLE = FALSE, parallel = FALSE, 
#                             BPPARAM = bpparam()) 
#   {
#     name = resultsNames(object)[length(resultsNames(object))]
#     
#     if (!"results" %in% mcols(mcols(object))$type) {
#       stop("cannot find results columns in object, first call DESeq, nbinomWaldTest, or nbinomLRT")
#     }
#     if (missing(test)) {
#       test <- attr(object, "test")
#     }
#     else if (test == "Wald" & attr(object, "test") == "LRT") {
#       object <- makeWaldTest(object)
#     }
#     else if (test == "LRT" & attr(object, "test") == "Wald") {
#       stop("the LRT requires the user run nbinomLRT or DESeq(dds,test='LRT')")
#     }
#     format <- match.arg(format, choices = c("DataFrame", "GRanges", 
#                                             "GRangesList"))
#     if (addMLE & !attr(object, "betaPrior")) {
#       stop("addMLE=TRUE is only for when a beta prior was used. otherwise, the log2 fold changes are already MLE")
#     }
#     if (format == "GRanges" & is(rowData(object), "GRangesList")) {
#       if (any(elementLengths(rowData(object)) == 0)) {
#         stop("rowData is GRangesList and one or more GRanges have length 0. Use format='DataFrame' or 'GRangesList'")
#       }
#     }
#     hasIntercept <- attr(terms(design(object)), "intercept") == 
#       1
#     isExpanded <- attr(object, "modelMatrixType") == "expanded"
#     termsOrder <- attr(terms.formula(design(object)), "order")
#     if ((test == "Wald") & (isExpanded | !hasIntercept) & missing(contrast) &  all(termsOrder < 2)) 
#       {
#       if (missing(name)) 
#       {
#         designVars <- all.vars(design(object))
#         lastVarName <- designVars[length(designVars)]
#         lastVar <- colData(object)[[lastVarName]]
#         if (is.factor(lastVar)) {
#           nlvls <- nlevels(lastVar)
#           contrast <- c(lastVarName, levels(lastVar)[nlvls], levels(lastVar)[1])
#         }
#       }
#     }
#     if (missing(name)) {
#       #name <- lastCoefName(object)
#       name=""
#     }
#     altHypothesis <- match.arg(altHypothesis, choices = c("greaterAbs", "lessAbs", "greater", "less"))
#     stopifnot(lfcThreshold >= 0)
#     stopifnot(length(lfcThreshold) == 1)
#     stopifnot(length(altHypothesis) == 1)
#     stopifnot(length(alpha) == 1)
#     stopifnot(length(pAdjustMethod) == 1)
#     stopifnot(length(listValues) == 2 & is.numeric(listValues))
#     stopifnot(listValues[1] > 0 & listValues[2] < 0)
#     if (length(name) != 1 | !is.character(name)) {
#       stop("the argument 'name' should be a character vector of length 1")
#     }
#     if (lfcThreshold == 0 & altHypothesis == "lessAbs") {
#       stop("when testing altHypothesis='lessAbs', set the argument lfcThreshold to a positive value")
#     }
#     print(names(mcols(object)))
#     WaldResults <- paste0("WaldPvalue_", name) %in% names(mcols(object))
#     print("OKG1.3")
#     LRTResults <- "LRTPvalue" %in% names(mcols(object))
#     print("OKG1.4")
#     if (!(WaldResults | LRTResults)) {
#       stop("cannot find appropriate results in the DESeqDataSet.\npossibly nbinomWaldTest or nbinomLRT has not yet been run.")
#     }
#     if (!missing(contrast)) {
#       if (is.character(contrast) & length(contrast) != 3) {
#         stop("'contrast', as a character vector of length 3, should have the form:\ncontrast = c('factorName','numeratorLevel','denominatorLevel'),\nsee the manual page of ?results for more information")
#       }
#       if (is.list(contrast) & length(contrast) == 1) {
#         contrast <- list(contrast[[1]], character())
#       }
#       if (is.list(contrast) & length(contrast) != 2) {
#         stop("'contrast', as a list, should have length 2,\nsee the manual page of ?results for more information")
#       }
#       if (is.list(contrast) & !(is.character(contrast[[1]]) & 
#                                   is.character(contrast[[2]]))) {
#         stop("'contrast', as a list of length 2, should have character vectors as elements,\nsee the manual page of ?results for more information")
#       }
#       print("OKG2")
#       res <- if (!parallel) {cleanContrast(object, contrast, expanded = isExpanded, listValues = listValues, test = test)}
#       else if (parallel) {
#         nworkers <- BPPARAM$workers
#         idx <- factor(sort(rep(seq_len(nworkers), length = nrow(object))))
#         do.call(rbind, bplapply(levels(idx), function(l) {
#           cleanContrast(object[idx == l, , drop = FALSE], 
#                         contrast, expanded = isExpanded, listValues = listValues, 
#                         test = test)
#         }, BPPARAM = BPPARAM))
#       }
#     }
#     else {
#       log2FoldChange <- getCoef(object, name)
#       lfcSE <- getCoefSE(object, name)
#       stat <- getStat(object, test, name)
#       pvalue <- getPvalue(object, test, name)
#       res <- cbind(mcols(object)["baseMean"], log2FoldChange, 
#                    lfcSE, stat, pvalue)
#       names(res) <- c("baseMean", "log2FoldChange", "lfcSE", 
#                       "stat", "pvalue")
#     }
#     print("OKG4")
#     rownames(res) <- rownames(object)
#     if (addMLE) {
#       if (!missing(contrast)) {
#         if (is.numeric(contrast)) 
#           stop("addMLE only implemented for: contrast=c('condition','B','A')")
#         if (is.list(contrast)) 
#           stop("addMLE only implemented for: contrast=c('condition','B','A')")
#         res <- cbind(res, mleContrast(object, contrast))
#       }
#       else {
#         mleName <- paste0("MLE_", name)
#         mleNames <- names(mcols(object))[grep("MLE_", names(mcols(object)))]
#         if (!mleName %in% mleNames) 
#           stop("MLE_ plus 'name' was not found as a column in mcols(dds)")
#         mleColumn <- mcols(object)[mleName]
#         names(mleColumn) <- "lfcMLE"
#         mcols(mleColumn)$description <- paste("log2 fold change (MLE):", 
#                                               name)
#         res <- cbind(res, mleColumn)
#       }
#       res <- res[, c("baseMean", "log2FoldChange", "lfcMLE", 
#                      "lfcSE", "stat", "pvalue")]
#     }
#     print("OKG5")
#     if (!(lfcThreshold == 0 & altHypothesis == "greaterAbs")) {
#       if (test == "LRT") {
#         warning("tests of log fold change above or below a theshold are Wald tests.\nLikelihood ratio test p-values are overwritten")
#       }
#       if (altHypothesis == "greaterAbs") {
#         newStat <- sign(res$log2FoldChange) * pmax(0, (abs(res$log2FoldChange) - 
#                                                          lfcThreshold))/res$lfcSE
#         newPvalue <- pmin(1, 2 * pnorm(abs(res$log2FoldChange), 
#                                        mean = lfcThreshold, sd = res$lfcSE, lower.tail = FALSE))
#       }
#       else if (altHypothesis == "lessAbs") {
#         if (attr(object, "betaPrior")) {
#           stop("testing altHypothesis='lessAbs' requires setting the DESeq() argument betaPrior=FALSE")
#         }
#         newStatAbove <- pmax(0, lfcThreshold - res$log2FoldChange)/res$lfcSE
#         pvalueAbove <- pnorm(res$log2FoldChange, mean = lfcThreshold, 
#                              sd = res$lfcSE, lower.tail = TRUE)
#         newStatBelow <- pmax(0, res$log2FoldChange + lfcThreshold)/res$lfcSE
#         pvalueBelow <- pnorm(res$log2FoldChange, mean = -lfcThreshold, 
#                              sd = res$lfcSE, lower.tail = FALSE)
#         newStat <- pmin(newStatAbove, newStatBelow)
#         newPvalue <- pmax(pvalueAbove, pvalueBelow)
#       }
#       else if (altHypothesis == "greater") {
#         newStat <- pmax(0, res$log2FoldChange - lfcThreshold)/res$lfcSE
#         newPvalue <- pnorm(res$log2FoldChange, mean = lfcThreshold, 
#                            sd = res$lfcSE, lower.tail = FALSE)
#       }
#       else if (altHypothesis == "less") {
#         newStat <- pmax(0, lfcThreshold - res$log2FoldChange)/res$lfcSE
#         newPvalue <- pnorm(res$log2FoldChange, mean = -lfcThreshold, 
#                            sd = res$lfcSE, lower.tail = TRUE)
#       }
#       print("OKG6")
#       res$stat <- newStat
#       res$pvalue <- newPvalue
#     }
#     m <- nrow(attr(object, "dispModelMatrix"))
#     p <- ncol(attr(object, "dispModelMatrix"))
#     if (m > p) {
#       defaultCutoff <- qf(0.99, p, m - p)
#       if (missing(cooksCutoff)) {
#         cooksCutoff <- defaultCutoff
#       }
#       stopifnot(length(cooksCutoff) == 1)
#       if (is.logical(cooksCutoff) & cooksCutoff) {
#         cooksCutoff <- defaultCutoff
#       }
#     }
#     else {
#       cooksCutoff <- FALSE
#     }
#     print("OKG7")
#     performCooksCutoff <- (is.numeric(cooksCutoff) | cooksCutoff)
#     if ((m > p) & performCooksCutoff) {
#       cooksOutlier <- mcols(object)$maxCooks > cooksCutoff
#       res$pvalue[cooksOutlier] <- NA
#     }
#     if (sum(mcols(object)$replace, na.rm = TRUE) > 0) {
#       nowZero <- which(mcols(object)$replace & mcols(object)$baseMean == 
#                          0)
#       res$log2FoldChange[nowZero] <- 0
#       if (addMLE) {
#         res$lfcMLE[nowZero] <- 0
#       }
#       res$lfcSE[nowZero] <- 0
#       res$stat[nowZero] <- 0
#       res$pvalue[nowZero] <- 1
#     }
#     print("OKG8")
#     if (independentFiltering) {
#       if (missing(filter)) {
#         filter <- res$baseMean
#       }
#       if (missing(theta)) {
#         lowerQuantile <- mean(filter == 0)
#         if (lowerQuantile < 0.95) 
#           upperQuantile <- 0.95
#         else upperQuantile <- 1
#         theta <- seq(lowerQuantile, upperQuantile, length = 20)
#       }
#       stopifnot(length(theta) > 1)
#       stopifnot(length(filter) == nrow(object))
#       filtPadj <- filtered_p(filter = filter, test = res$pvalue, 
#                              theta = theta, method = pAdjustMethod)
#       numRej <- colSums(filtPadj < alpha, na.rm = TRUE)
#       j <- which.max(numRej)
#       res$padj <- filtPadj[, j, drop = TRUE]
#       cutoffs <- quantile(filter, theta)
#       filterThreshold <- cutoffs[j]
#       filterNumRej <- data.frame(theta = theta, numRej = numRej)
#     }
#     else {
#       res$padj <- p.adjust(res$pvalue, method = pAdjustMethod)
#     }
#     mcols(res)$type[names(res) == "padj"] <- "results"
#     mcols(res)$description[names(res) == "padj"] <- paste(pAdjustMethod, "adjusted p-values")
#     
#     print("OKG9")
#     deseqRes <- DESeqResults(res)
#     if (independentFiltering) {
#       attr(deseqRes, "filterThreshold") <- filterThreshold
#       attr(deseqRes, "filterNumRej") <- filterNumRej
#     }
#     if (format == "DataFrame") {
#       return(deseqRes)
#     }
#     else if (format == "GRangesList") {
#       if (class(rowData(object)) == "GRanges") 
#         message("rowData is GRanges")
#       out <- rowData(object)
#       mcols(out) <- deseqRes
#       return(out)
#     }
#     else if (format == "GRanges") {
#       if (class(rowData(object)) == "GRangesList") {
#         message("rowData is GRangesList, unlisting the ranges")
#         out <- unlist(range(rowData(object)))
#         mcols(out) <- deseqRes
#         return(out)
#       }
#       else {
#         out <- rowData(object)
#         mcols(out) <- deseqRes
#         return(out)
#       }
#     }
#   }
#   
#   
#   TableDiff_result <- function(input,BaseContrast,resDiff)
#   {
#     VarInt = input$VarInt
#     dds = resDiff$dds
#     counts = resDiff$counts
#     target = resDiff$target
#     group = as.data.frame(target[,VarInt])
#     result = list()
#     
#     result[[input$ContrastList_table]] <- GetTableFromContrast(dds,contrast=BaseContrast[,input$ContrastList_table],pAdjustMethod=input$AdjMeth,
#                                cooksCutoff=ifelse(!is.null(input$CooksCutOff),input$CutOffVal,TRUE),
#                                independentFiltering=input$IndFiltering,alpha=input$AlphaVal)
#     return(result)
#   }
#   
  
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
      if (is.null(cooksCutoff)) {
        m <- nrow(attr(dds, "modelMatrix"))
        p <- ncol(attr(dds, "modelMatrix"))
        cooksCutoff <- qf(0.99, p, m - p)
      }
      mcols.add$outlier <- ifelse(mcols(dds)$maxCooks > cooksCutoff, "Yes", "No")
      complete.name <- merge(complete.name, mcols.add, by = "Id")
      complete[[name]] <- complete.name
      up.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange >= 0), ]
      up.name <- up.name[order(up.name$padj), ]
      down.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange <= 0), ]
      down.name <- down.name[order(down.name$padj), ]
      name <- gsub(" ", "", name)
#       keep <- c("FC", "log2FoldChange", "padj")
#       complete.complete[, paste(name, keep, sep = ".")] <- complete.name[, keep]
    }
    #return(list(complete=complete.name,up=up.name,down=down.name))
    return(list(complete=complete.name[,c("Id","baseMean","FC","log2FoldChange","padj")],up=up.name[,c("Id","baseMean","FC","log2FoldChange","padj")],down=down.name[,c("Id","baseMean","FC","log2FoldChange","padj")]))
  }
  
  
  