#@ This file contains all the functions for the 
#@ diagnostic plots of SHAMAN


##############################################################
##
##        Main function for the diagnostic plots
##
##############################################################

Plot_diag <- function(input,resDiff,tree,getTable=FALSE)
{
  
  VarInt = input$VarInt
  dds = resDiff$dds
  counts = resDiff$raw_counts
  if(input$CountsType=="Normalized") counts = resDiff$countsNorm
  target = resDiff$target
  normFactors = resDiff$normFactors
  
  ## Counts at the OTU level
  CT = resDiff$CT_noNorm
  if(input$CountsType=="Normalized") CT = resDiff$CT_Norm
  
  group = as.data.frame(target[,VarInt])
  rownames(group) = rownames(target)
  res = NULL
  
  if(ncol(group)>0 && nrow(counts)>0 && !getTable)
  { 
    colors = rep(c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78", "#2ca02c","#98df8a","#d62728","#ff9896","#9467bd","#c5b0d5","#8c564b",
                   "#c49c94","#e377c2","#f7b6d2","#7f7f7f", "#c7c7c7","#bcbd22","#dbdb8d","#17becf","#9edae5"),ceiling(nrow(target)/20))
    
    if(input$DiagPlot=="barplotTot") res = barplotTot(input,counts,group = group, col=colors)
    if(input$DiagPlot=="barplotNul") res = barPlotNul(input,counts, group = group, col=colors)
    if(input$DiagPlot=="densityPlot") res = densityPlotTot(input,counts, group = group, col=colors)
    if(input$DiagPlot=="boxplotNorm") res = boxplotNorm(input,CT,group = group, col=colors)
    if(input$DiagPlot=="DispPlot") res = plotDispEsts(dds)
    if(input$DiagPlot=="MajTax") res = majTaxPlot(input,counts, group = group, col=colors)
    if(input$DiagPlot=="SfactorsVStot") res = diagSFactors(input,normFactors,resDiff$raw_counts) 
    if(input$DiagPlot=="pcaPlot") res = PCAPlot_meta(input,dds, group,  type.trans = input$TransType, col = colors)
    if(input$DiagPlot=="pcoaPlot") res = PCoAPlot_meta(input,dds, group,CT,tree, col = colors) 
    if(input$DiagPlot=="nmdsPlot") res = NMDSPlot(input, dds, group,CT,tree, col = colors) 
    if(input$DiagPlot=="clustPlot") res = HCPlot(input,dds,group,type.trans=input$TransType,counts,CT,tree,col=colors)
  }
  if(getTable && input$DiagPlot=="pcaPlot") res = Get_pca_table(input,dds, group,  type.trans = input$TransType)
  if(getTable && input$DiagPlot=="pcoaPlot") res = Get_pcoa_table(input,dds, group,CT,tree)$table
  if(getTable && input$DiagPlot=="nmdsPlot") res = Get_nmds_table(input,dds, group,CT,tree)$table
  return(res)
}




##############################################################
##
##        Permanova test
##
##############################################################

Perma_test_Diag <- function(input,resDiff,tree)
{
  
  VarInt = input$VarInt
  dds = resDiff$dds
  
  ## Counts at the OTU level
  CT = resDiff$CT_noNorm
  if(input$CountsType=="Normalized") CT = resDiff$CT_Norm
  
  target = resDiff$target
  normFactors = resDiff$normFactors
  
  group = as.data.frame(target[,VarInt])
  rownames(group) = rownames(target)
  res = NULL
  
  if(ncol(group)>0 && !is.null(dds))
  { 
    if(input$DiagPlot=="pcoaPlot") res = Get_pcoa_table(input,dds, group,CT,tree)$test$aov.tab$`Pr(>F)`[1]
    if(input$DiagPlot=="nmdsPlot") res = Get_nmds_table(input,dds, group,CT,tree)$test$aov.tab$`Pr(>F)`[1]
  }

  return(res)
}


##############################################################
##
##          Plot functions 
##
##############################################################


## Hierarchical clustering
HCPlot <- function (input,dds,group,type.trans=NULL,counts=NULL,CT,tree,col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen")) 
{
  
  res = NULL
  if(!is.null(dds) && !is.null(counts) && !is.null(type.trans) && !is.null(input$DistClust)){
    ## Get the counts

    if (input$DistClust == "euclidean" && type.trans == "VST") counts <- assay(varianceStabilizingTransformation(dds))
    if (input$DistClust == "euclidean" && type.trans == "rlog") counts <- assay(rlogTransformation(dds))
    
    ## Get the group of leaf
    group = apply(group,1,paste, collapse = "-")    
    nb = length(unique((group)))
    
    ## Get the dendrogram
    if(input$DistClust=="sere") dist.counts.norm = as.dist(SEREcoef(counts))
    else if(input$DistClust=="Unifrac"){
      tmp = UniFracDist(CT,tree)
      if(is.null(tree) || is.null(tmp)) dist.counts.norm = NULL
      else
      {
        dist.counts.norm = switch(input$DistClustUnifrac,
                                  "WU" = as.dist(tmp[, , "d_1"]), 
                                  "UWU" = as.dist(tmp[, , "d_UW"]),
                                  "VAWU" = as.dist(tmp[, , "d_VAW"]))
      }
      
    }
    else if(input$DistClust  %in% getDistMethods()){
      dist = as.dist(distance(t(sweep(counts,2,colSums(counts),`/`)), method=input$DistClust))
      dist[is.na(dist)]=0.0
      dist.counts.norm = dist
    }
    else  dist.counts.norm = vegdist(t(counts), method = input$DistClust)
  
    if(!is.null(dist.counts.norm))
    {
  
      hc <- hclust(dist.counts.norm, method = "ward.D")
      
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
      else
      {
        par(cex=input$cexTitleDiag,mar=c(6,6,4,5))
        res = circlize_dendrogram(dend, labels_track_height = 0.2, dend_track_height = .3, main = "Cluster dendrogram",xlab = paste(input$DistClust,"distance, Ward criterion",sep=" "))
      }
    }
  }
  return(res)
}

## PCA Eigen value
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


## PCOA Eigen value
Plot_diag_pcoaEigen = function(input,resDiff,tree)
{
  colors = c("SpringGreen","dodgerblue","black","firebrick1")
  VarInt = input$VarInt
  dds = resDiff$dds
  target = resDiff$target
  group = as.data.frame(target[,VarInt])
  rownames(group) = rownames(target)
  CT = resDiff$CT_noNorm
  if(input$CountsType=="Normalized") CT = resDiff$CT_Norm
  PCoAPlot_meta(input,dds, group,CT=CT,tree, col = colors, plot = "eigen") 
}


## NMDS stress plot
Plot_diag_nmdsStress = function(input,resDiff,tree)
{
  VarInt = input$VarInt
  dds = resDiff$dds
  target = resDiff$target
  group = as.data.frame(target[,VarInt])
  rownames(group) = rownames(target)
  CT = resDiff$CT_noNorm
  if(input$CountsType=="Normalized") CT = resDiff$CT_Norm
  NMDSPlot(input,dds, group,CT=CT,tree=tree, plot = "stress") 
}




## Boxplot for the counts normalized/no normalized
boxplotNorm <- function(input,CT, group, col = c("lightblue","orange", "MediumVioletRed", "SpringGreen")) 
{
  
  ncol1 <- ncol(group) == 1
  par(cex=input$cexTitleDiag,mar=c(12,6,4,5))
  if(input$RemoveNullValue) CT[CT==0] = NA
  
  ### Boxplots of the counts
  my.boxplot(log2(CT+1), las = 2, pol.col = col[as.integer(group[,1])],
             pol.density = if (ncol1) {NULL}
             else {15}, 
             pol.angle = if (ncol1) {NULL}
             else {seq(0,160,length.out =nlevels(group[, 2]))[as.integer(group[, 2])]},
             main = paste(input$CountsType, "counts distribution"), ylab = expression(log[2] ~ ( count + 1)))
  legend("topright", levels(group[, 1]), fill = col[1:nlevels(group[,1])], bty = "n")
  if (!ncol1)  legend("topleft", levels(group[, 2]), density = 15,col = 1, angle = seq(0,160,length.out =nlevels(group[, 2]))[1:nlevels(group[, 2])], bty = "n")
  
}

## barplot total
barplotTot <- function(input,counts, group, cex.names = 1, col = c("lightblue","orange", "MediumVioletRed", "SpringGreen")) 
{
  
  ncol1 <- ncol(group) == 1
  par(cex=input$cexTitleDiag,mar=c(12,6,4,5))
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
  
  par(cex=input$cexTitleDiag,mar=c(12,6,4,5))
  
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
  par(cex=input$cexTitleDiag,mar=c(12,5,4,5))
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
  par(cex=input$cexTitleDiag,mar=c(12,6,4,5))
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
PCoAPlot_meta <-function (input, dds, group_init, CT,tree,col = c("SpringGreen","dodgerblue","black","firebrick1"), plot = "pcoa") 
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
  
  if(nlevels(group)!=0 && !is.null(input$PCaxe1) && !is.null(input$PCaxe2))
  { 
    ## Get the norm data
    counts.norm = as.data.frame(round(counts(dds)))
    if(input$CountsType=="Normalized") counts.norm = as.data.frame(round(counts(dds, normalized = TRUE)))
    # was removed
    counts.norm = counts.norm[,ind_kept]
    # print(head(counts.norm))
    ## Get the distance
    if(input$DistClust=="sere") dist.counts.norm = as.dist(SEREcoef(counts.norm))
    else if(input$DistClust=="Unifrac"){
      #tmp = UniFracDist(CT,tree)
      tmp = UniFracDist(counts.norm,tree)
      if(is.null(tree) || is.null(tmp)) dist.counts.norm = NULL
      else
      {
        dist.counts.norm = switch(input$DistClustUnifrac,
                                  "WU" = as.dist(tmp[, , "d_1"]), 
                                  "UWU" = as.dist(tmp[, , "d_UW"]),
                                  "VAWU" = as.dist(tmp[, , "d_VAW"]))
      }
      
    }
    else if(input$DistClust  %in% getDistMethods()){
      dist = as.dist(distance(t(sweep(counts.norm,2,colSums(counts.norm),`/`)), method=input$DistClust))
      dist[is.na(dist)]=0.0
      dist.counts.norm = dist
    }
    else  dist.counts.norm = vegdist(t(counts.norm), method = input$DistClust)
    #"additive_symm"
    if(!is.null(dist.counts.norm))
    {
      ## Do PCoA
      pco.counts.norm = dudi.pco(d = dist.counts.norm, scannf = FALSE,nf=ncol(counts.norm))
      
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
        v_axes = c(as.numeric(gsub("PC","",input$PCaxe1)),as.numeric(gsub("PC","",input$PCaxe2)))
        
        plot(pco.counts.norm$li[v_axes], 
             xlab=paste(input$PCaxe1, ": ",round(eigen[v_axes[1]],1),"%") , 
             ylab=paste(input$PCaxe2, ": ",round(eigen[v_axes[2]],1),"%"),
             xlim=c(min+0.25*min,max+0.25*max), ylim=c(min-0.1,max+0.1), 
             cex.axis=1, cex.lab=1,lwd=2, type="n")
        
        title(main='Principal Coordinates Analysis ',cex.main=1.5)
        ## Add a subtitle
        par(font.main=3)
        if(input$DistClust=="Unifrac") title(main=paste("\n","\n",input$DistClustUnifrac,"distance",sep=" "),cex.main=1)
        else title(main=paste("\n","\n",input$DistClust,"distance",sep=" "),cex.main=1)
        
        # Set different shapes
        if(input$labelPCOA == "Group"){
          if(!is.null(cval)){
            for (i in 1:length(cval)){
              points(pco.counts.norm$li[which(group==cval[i]),v_axes],pch=shape[i],col=col[i], cex=input$cexpoint)
            }
            s.class(dfxy = pco.counts.norm$li[v_axes], fac = group, col = col, label = levels(group),
                    add.plot = TRUE, cpoint = 0, cell=input$cexcircle, clabel=input$cexLabelDiag,  cstar = input$cexstar)
          }else s.class(dfxy = pco.counts.norm$li[v_axes], fac = group, col = col, label = levels(group),
                        add.plot = TRUE, cpoint = input$cexpoint, cell=input$cexcircle, clabel=input$cexLabelDiag,  cstar = input$cexstar)
        }  
        else{
          s.label(pco.counts.norm$li, clabel = input$cexLabelDiag,boxes=FALSE, add.plot = TRUE)
          s.class(dfxy = pco.counts.norm$li, fac = group, col = col, label = levels(group), add.plot = TRUE, cpoint = 0, clabel = 0, cstar = input$cexstar, cell=input$cexcircle)
        }
      }else{
        v_axes = c(as.numeric(gsub("PC","",input$PCaxe1)),as.numeric(gsub("PC","",input$PCaxe2)))
        nbBar = max(7,max(v_axes))
        col = rep("grey",nbBar)
        col[v_axes] = "black"
        barplot(eigen[1:nbBar], xlab="Dimensions", ylab="Eigenvalues (%)", names.arg = 1:nbBar, col = col, ylim=c(0,max(eigen)+5), cex.axis=1.2, cex.lab=1.4,cex.names=1.2)
      }
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
  
  # save(val,Kval,dds,group_init,type.trans,VarInt,ind_kept,file="testLDA")
  ## Get the group corresponding to the modalities
  group = group[ind_kept]
  nb = length(unique((group)))
  group = as.factor(group)
  
  ## To select the colors
  indgrp =as.integer(as.factor(group_init))[ind_kept]
  
  
  if(nlevels(group)!=0 && !is.null(input$PCaxe1) && !is.null(input$PCaxe2))
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

      abs = range(pca$x[, as.numeric(gsub("PC","",input$PCaxe1))])
      abs = abs(abs[2] - abs[1])/25
      ord = range(pca$x[, as.numeric(gsub("PC","",input$PCaxe2))])
      ord = abs(ord[2] - ord[1])/25
      
      
      plot(pca$x[,as.numeric(gsub("PC","",input$PCaxe1))],pca$x[,as.numeric(gsub("PC","",input$PCaxe2))], las = 1, cex = input$cexTitleDiag, col = col[indgrp], 
           pch = 16,
           xlab = paste0(input$PCaxe1," (", prp[as.numeric(gsub("PC","",input$PCaxe1))], "%)"),
           ylab = paste0(input$PCaxe2," (", prp[as.numeric(gsub("PC","",input$PCaxe2))], "%)"),
           main = "Principal Component Analysis"
      )
      abline(h = 0, v = 0, lty = 2, col = "lightgray")
      text(pca$x[, as.numeric(gsub("PC","",input$PCaxe1))] - ifelse(pca$x[, as.numeric(gsub("PC","",input$PCaxe1))] > 0, abs, -abs), pca$x[,as.numeric(gsub("PC","",input$PCaxe2))] - ifelse(pca$x[,as.numeric(gsub("PC","",input$PCaxe2))] > 0, ord, -ord), colnames(counts.trans), col = col[indgrp],cex=input$cexLabelDiag)
      
    }
    if(plot=="eigen"){
      nbBar = max(7,max(c(as.numeric(gsub("PC","",input$PCaxe1)),as.numeric(gsub("PC","",input$PCaxe2)))))
      col = rep("grey",nbBar)
      eigen = pca$sdev[1:nbBar]^2
      col[c(as.numeric(gsub("PC","",input$PCaxe1)),as.numeric(gsub("PC","",input$PCaxe2)))] = "black"
      barplot(eigen, xlab="Dimensions", ylab="Eigenvalues (%)", names.arg = 1:nbBar, col = col, ylim=c(0,max(eigen)+5), cex.axis=1.2, cex.lab=1.4,cex.names=1.2)}
    
  }
}



##############################################################
##
##          Useful functions
##
##############################################################

## Remove nul counts
removeNulCounts <-function (counts) 
{
  return(counts[rowSums(counts) > 0, ])
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

## Create your own boxplots to get hatched boxplots
my.boxplot <- function(x, pol.col = 1, pol.density = NULL, pol.angle = 45,
                       bxp.pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5), ...){
  res <- boxplot(x, pars = bxp.pars, ...) # que boxplot se d�merde avec ses arguments
  n <- ncol(res$stats) # nombre de boxplots
  density <- if(is.null(pol.density)){NULL}else{rep(pol.density, length = n)}
  angle <- if(is.null(pol.angle)){NULL}else{rep(pol.angle, length = n)}
  col <- if(is.null(pol.col)){NULL}else{rep(pol.col, length = n)}
  # Ajout des textures
  ex <- bxp.pars$boxwex/2 
  for(i in 1:n){
    polygon(c(i - ex, i - ex, i + ex, i + ex),
            c(res$stats[2, i], res$stats[4, i], res$stats[4, i], res$stats[2, i]),
            density = density[i], angle = angle[i], col = col[i])
    segments(i-ex,res$stats[3,i],i+ex,res$stats[3,i],lwd=3,col="black",lend=1)
  }
}


### Get PCOA table (useful to get the number of axes)
Get_pcoa_table <-function (input, dds, group_init,CT,tree) 
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
  permanova_test = NULL
  
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
    if(input$CountsType=="Normalized") counts.norm = as.data.frame(round(counts(dds, normalized = TRUE)))
    # was removed
    counts.norm = counts.norm[,ind_kept]
    
    ## Get the distance
    if(input$DistClust=="sere") dist.counts.norm = as.dist(SEREcoef(counts.norm))
    else if(input$DistClust=="Unifrac") {
      tmp = UniFracDist(counts.norm,tree)
      if(is.null(tree) || is.null(tmp)) dist.counts.norm = NULL
      else{
          dist.counts.norm = switch(input$DistClustUnifrac,
                                "WU" = as.dist(tmp[, , "d_1"]), 
                                "UWU" = as.dist(tmp[, , "d_UW"]),
                                "VAWU" = as.dist(tmp[, , "d_VAW"]))
      }
    }
    else if(input$DistClust  %in% getDistMethods()){
      dist = as.dist(distance(t(sweep(counts.norm,2,colSums(counts.norm),`/`)), method=input$DistClust))
      dist[is.na(dist)]=0.0
      dist.counts.norm = dist
    }
    else {
      dist.counts.norm = vegdist(t(counts.norm), method = input$DistClust)
    }
    
    if(!is.null(dist.counts.norm)){ 
      ## To get always the same result 
      set.seed(666)
      permanova_test = adonis(dist.counts.norm~group)
  
      ## Do PCoA
      pco.counts.norm = dudi.pco(d = dist.counts.norm, scannf = FALSE,nf=ncol(counts.norm))
      return(list(table = pco.counts.norm$li,test=permanova_test))
    } else return(list(table = NULL,test=NULL))
    
  }
}

### Get PCA table (useful to get the number of axes)
Get_pca_table <-function(input,dds, group_init, n = min(500, nrow(counts(dds))), type.trans = c("VST", "rlog")) 
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
    return(pca$x)
  }
}


Get_nmds_table <-function(input,dds, group_init,CT,tree) 
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
    ## Get the norm data
    counts.norm = as.data.frame(round(counts(dds)))
    if(input$CountsType=="Normalized") counts.norm = as.data.frame(round(counts(dds, normalized = TRUE)))
    # was removed
    counts.norm = counts.norm[,ind_kept]
    
    ## Get the distance
    if(input$DistClust!="sere" && input$DistClust!="Unifrac") dist.counts.norm = vegdist(t(counts.norm), method = input$DistClust)
    if(input$DistClust=="sere") dist.counts.norm = as.dist(SEREcoef(counts.norm))
    if(input$DistClust=="Unifrac") {
      tmp = UniFracDist(CT,tree)
      if(is.null(tree) || is.null(tmp)) dist.counts.norm = NULL
      if(!is.null(tree)) {dist.counts.norm = switch(input$DistClustUnifrac,
                                                    "WU" = as.dist(tmp[, , "d_1"]), 
                                                    "UWU" = as.dist(tmp[, , "d_UW"]),
                                                    "VAWU" = as.dist(tmp[, , "d_VAW"])
      )
      }
    
    }
    
    if(!is.null(dist.counts.norm)){ 
      
      ## To get always the same result 
      set.seed(666)
      permanova_test = adonis(dist.counts.norm~group)
      
      ## Do NMDS
      # nmds.counts.norm = metaMDS(dist.counts.norm,k=min(round((nrow(counts.norm)-1)/2-1),round(ncol(counts.norm)/2)), trymax = 1)
      nmds.counts.norm = metaMDS(dist.counts.norm,k=3, trymax = 1)
      proj = nmds.counts.norm$points
      return(list(table = proj,test=permanova_test,nmds=nmds.counts.norm))
      
    } else return(list(table = NULL,test=NULL))
    
  }
}





### NMDS
NMDSPlot <-function (input, dds, group_init, CT,tree,col = c("SpringGreen","dodgerblue","black","firebrick1"),plot="nmds") 
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
  
  
  if(nlevels(group)!=0 && !is.null(input$PCaxe1) && !is.null(input$PCaxe2))
  { 
    
    ## Get the norm data
    counts.norm = as.data.frame(round(counts(dds)))
    if(input$CountsType=="Normalized") counts.norm = as.data.frame(round(counts(dds, normalized = TRUE)))
    # was removed
    counts.norm = counts.norm[,ind_kept]
    
    ## Get the distance
    if(input$DistClust!="sere" && input$DistClust!="Unifrac") dist.counts.norm = vegdist(t(counts.norm), method = input$DistClust)
    if(input$DistClust=="sere") dist.counts.norm = as.dist(SEREcoef(counts.norm))
    if(input$DistClust=="Unifrac") {
      tmp = UniFracDist(CT,tree)
      if(is.null(tree) || is.null(tmp)) dist.counts.norm = NULL
      if(!is.null(tree)) {dist.counts.norm = switch(input$DistClustUnifrac,
                                                    "WU" = as.dist(tmp[, , "d_1"]), 
                                                    "UWU" = as.dist(tmp[, , "d_UW"]),
                                                    "VAWU" = as.dist(tmp[, , "d_VAW"])
      )
      }
    
    }

    if(!is.null(dist.counts.norm)){ 
      
      ## Do NMDS
      # nmds.counts.norm = metaMDS(dist.counts.norm,k=min(round((nrow(counts.norm)-1)/2-1),round(ncol(counts.norm)/2)), trymax = 25)
      nmds.counts.norm = metaMDS(dist.counts.norm,k=3, trymax = 25)
      if(plot =="nmds") {

        proj = nmds.counts.norm$points
        
        ## xlim and ylim of the plot
        min = min(proj); max = max(proj)
        abs = range(proj[, as.numeric(gsub("PC","",input$PCaxe1))])
        abs = abs(abs[2] - abs[1])/25
        ord = range(proj[, as.numeric(gsub("PC","",input$PCaxe2))])
        ord = abs(ord[2] - ord[1])/25
        
        plot(proj[,as.numeric(gsub("PC","",input$PCaxe1))],proj[,as.numeric(gsub("PC","",input$PCaxe2))], las = 1, cex = input$cexTitleDiag, col = col[indgrp], 
             pch = 16,
             xlab = paste("MDS",gsub("PC","",input$PCaxe1)),
             ylab = paste("MDS",gsub("PC","",input$PCaxe2))
        )
        title(main= "Non-metric multidimensional scaling",cex.main=1.5)
        par(font.main=3)
        if(input$DistClust=="Unifrac") title(main=paste("\n","\n",input$DistClustUnifrac,"distance",sep=" "),cex.main=1)
        else title(main=paste("\n","\n",input$DistClust,"distance",sep=" "),cex.main=1)
        abline(h = 0, v = 0, lty = 2, col = "lightgray")
        text(proj[,as.numeric(gsub("PC","",input$PCaxe1))] - ifelse(proj[,as.numeric(gsub("PC","",input$PCaxe1))] > 0, abs, -abs), proj[,as.numeric(gsub("PC","",input$PCaxe2))] - ifelse(proj[,as.numeric(gsub("PC","",input$PCaxe2))] > 0, ord, -ord), colnames(counts.norm), col = col[indgrp],cex=input$cexLabelDiag)
      }
      if(plot=="stress") stressplot(nmds.counts.norm,pch=16,p.col="black")
    }
  }
  
}



### Unifrac distance 
UniFracDist<-function(counts,tree)
{
  counts_otu = round(t(as.matrix(counts)),0)
  # unifracs = vegdist(counts_otu) ## Default value
  unifracs = NULL
  if(!is.null(counts_otu) && !is.null(tree)){
    unifracs <- GUniFrac(counts_otu, tree, alpha=c(0, 0.5, 1))$unifracs
  }
  return(unifracs) 
}


