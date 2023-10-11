#@ This file contains all the functions for the 
#@ diagnostic plots of SHAMAN


##############################################################
##
##        Main function for the diagnostic plots
##
##############################################################

Plot_diag <- function(input,resDiff,tree = NULL,getTable=FALSE, pca_calcul = NULL)
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
      # colors = switch(input$colorsdiag,
      #                 "shaman-palette1"=rep(c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78", "#2ca02c","#98df8a","#d62728","#ff9896","#9467bd","#c5b0d5","#8c564b",
      #                "#c49c94","#e377c2","#f7b6d2","#7f7f7f", "#c7c7c7","#bcbd22","#dbdb8d","#17becf","#9edae5"),ceiling(nrow(target)/20)),
      #                "shaman-palette2"=rep(c("#4c4c4c", "#bdbdbd", "#ff4cff", "#bc62f4", "#1613ef", "#1391ef", "#ef1340"), ceiling(nrow(target)/7))
      # )
      colors = switch(input$colorsdiag,
                      "retro palette" = rep(c(
                        "#048789", "#503D2E", "#D44D27", "#E2A72E", "#EFEBC8",
                                 "#107E7D", "#7A3D3D", "#F15A22", "#F7D488", "#F4A259",
                                 "#005F6B", "#4D314A", "#BF6B63", "#FF8C42", "#FF3C38"
                      ), ceiling(nrow(target)/15)),
                      "easter palette" = rep(c(
                        "#DED4FF", "#A6E7FF", "#C7FFCD", "#FFF7AD", "#FFDEBA",
                                 "#D1F2A5", "#EFFAB4", "#FFC48C", "#FF9F80", "#FF6B6B",
                                 "#69D2E7", "#A7DBD8", "#E0E4CC", "#F38630", "#FA6900"
                      ), ceiling(nrow(target)/15)),
                      "warm palette" = rep(c(
                        "#FFCC0D", "#FF7326", "#FF194D", "#BF2669", "#702A8C",
                                 "#FFD700", "#FF8C00", "#FF4500", "#FF6347", "#FF1493",
                                 "#FF69B4", "#FFA07A", "#FF7F50", "#CD5C5C", "#8B0000"
                      ), ceiling(nrow(target)/15)),
                      "basic palette (1)" = rep(c(
                        "#f44336", "#e81e63", "#9c27b0", "#673ab7", "#3f51b5",
                                 "#2196f3", "#03a9f4", "#00bcd4", "#009688", "#4caf50",
                                 "#8bc34a", "#cddc39", "#ffeb3b", "#ffc107", "#ff9800",
                                 "#ff5722", "#795548", "#9e9e9e", "#607d8b", "#000000"
                      ), ceiling(nrow(target)/20)),
                      "basic palette (2)" = rep(c(
                        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
                                 '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
                                 '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
                                 '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080',
                                 '#ffffff', '#000000'
                      ), ceiling(nrow(target)/20)),
                      "basic palette (3)" = rep(c(
                        '#C0362C', '#FF8642', '#F4DCB5', '#816C5B', '#C3B7AC',
                                 '#E7E3D7', '#668D3C', '#B1DDA1', '#E5F3CF', '#0097AC',
                                 '#3CD6E6', '#97EAF4', '#007996', '#06C2F4', '#FAD8FA'
                      ), ceiling(nrow(target)/15)),
                      "basic palette (4)" = rep(c(
                        "#6929c4", "#1192e8", "#005d5d", "#9f1853", "#fa4d56",
                                 "#570408", "#198038", "#002d9c", "#ee538b", "#b28600",
                                 "#009d9a", "#012749", "#8a3800", "#a56eff"
                      ), ceiling(nrow(target)/14))
      )
      
      
      
      if(input$DiagPlot=="barplotTot") res = barplotTot(input,counts,group = group, col=colors)
      if(input$DiagPlot=="barplotNul") res = barPlotNul(input,counts, group = group, col=colors)
      if(input$DiagPlot=="densityPlot") res = densityPlotTot(input,counts, group = group, col=colors)
      if(input$DiagPlot=="boxplotNorm") res = boxplotNorm(input,CT,group = group, col=colors)
      if(input$DiagPlot=="DispPlot") res = plotDispEsts(dds)
      if(input$DiagPlot=="MajTax") res = majTaxPlot(input,counts, group = group, col=colors)
      if(input$DiagPlot=="SfactorsVStot") res = diagSFactors(input,normFactors,resDiff$raw_counts) 
      if(input$DiagPlot=="pcaPlot") res = PCAPlot_meta(input,dds, group,  type.trans = input$TransType, col = colors, pca_res = pca_calcul)
      if(input$DiagPlot=="pcoaPlot") res = PCoAPlot_meta(input,dds, group,CT,tree, col = colors, resDiff = resDiff)$plot 
      if(input$DiagPlot=="nmdsPlot") res = NMDSPlot(input, dds, group,CT,tree, col = colors) 
      if(input$DiagPlot=="clustPlot") res = HCPlot(input,dds,group,type.trans=input$TransType,counts,CT,tree,col=colors)
    }
    if(getTable && input$DiagPlot=="pcaPlot") res = Get_pca_table(input,dds, group,  type.trans = input$TransType)$table
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
    if(input$DiagPlot=="pcaPlot") res = Get_pca_table(input,dds, group)$test$aov.tab$`Pr(>F)`[1]
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
      
      #dend <- set(dend, "labels_cex", input$cexLabelDiag)
      if(input$colorHC) labels_colors(dend)<-col[as.integer(as.factor(group))][order.dendrogram(dend)]
      if(type=="hori") 
      { 
        par(cex=input$cexTitleDiag,mar=c(6,6,4,5))
        res = plot(dend, main = "Cluster dendrogram",xlab = paste(input$DistClust,"distance, Ward criterion",sep=" "),cex=input$cexLabelDiag)
      }  
      else
      {
        par(cex=input$cexTitleDiag,mar=c(6,6,4,5))
        res = circlize_dendrogram(dend, labels_track_height = 0.2, dend_track_height = .3, main = "Cluster dendrogram",xlab = paste(input$DistClust,"distance, Ward criterion",sep=" "), 
                                  labels_cex = input$cexLabelDiag)
      }
    }
  }
  return(res)
}

## PCA Eigen value
Plot_diag_Eigen <- function(input,resDiff, pca)
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
  
  res = PCAPlot_meta(input,dds, group,  type.trans = input$TransType, col = colors, plot = "eigen", pca_res = pca)
  return(res)
}

Plot_diag_Contrib <- function(input, resDiff, pca)
{
  colors = switch(input$colorsdiag,
                      "retro palette" = rep(c(
                        "#048789", "#503D2E",  "#FF3C38"
                      )),
                      "easter palette" = rep(c(
                        "#DED4FF", "#A6E7FF", "#C7FFCD"
                      )),
                      "warm palette" = rep(c(
                        "#FF6347", "#FF194D", "#702A8C"
                                 ))
              )
  VarInt = input$VarInt
  dds = resDiff$dds
  counts = resDiff$counts
  target = resDiff$target
  group = as.data.frame(target[,VarInt])
  
  ## If more than 4 levels for one factor
  maxFact =max(sapply(group,FUN=function(x) length(levels(x))))
  if(maxFact>=4) colors = rainbow(maxFact)
  
  res = PCAPlot_meta(input,dds, group,  type.trans = input$TransType, col = colors, plot = "contrib", pca_res = pca)
  return(res)
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
  PCoAPlot_meta(input,dds, group,CT=CT,tree, col = colors, plot = "eigen", resDiff = resDiff) 
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

compute_arrows <-  function(given_pcoa, trait_df) {
  
  n <- nrow(trait_df)
  points.stand <- scale(given_pcoa$vectors)
  
  # Compute covariance of variables with all axes
  S <- cov(trait_df, points.stand)
  
  # Select only positive eigenvalues
  pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))]
  
  # Standardize value of covariance (see Legendre & Legendre 1998)
  U <- S %*% diag((pos_eigen/(n - 1))^(-0.5))
  colnames(U) <- colnames(given_pcoa$vectors)
  
  # Add values of covariances inside object
  given_pcoa$U <- U
  
  return(given_pcoa)
}


### PCoA
PCoAPlot_meta <-function (input, dds, group_init, CT,tree,col = c("SpringGreen","dodgerblue","black","firebrick1"), plot = "pcoa", pcoa_theme = pca_themes(), resDiff) 
{
  dist.counts.norm = NULL
  target = resDiff$target
  pp = NULL
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
  group = factor(group, levels = Kval)
  
  if(nlevels(group)!=0 && !is.null(input$PCaxe1) && !is.null(input$PCaxe2))
  { 
    ## Get the norm data
    counts.norm = as.data.frame(round(counts(dds)))
    if(input$CountsType=="Normalized") counts.norm = as.data.frame(round(counts(dds, normalized = TRUE)))
    # was removed
    counts.norm = counts.norm[,ind_kept]
    # print(head(counts.norm))
    ## Get the distance
    if(input$DiagPlot == 'pcaPlot'){
      res = Get_pca_table(input, dds, group_init)$table
      return(list(plot = NULL, dataDiv = res))
    }
    else{
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
      
      # pcoa_df <- data.frame(pco.counts.norm$li) 
      # 
      # 
      # pcoa_df <- pcoa_df %>% 
      #   tibble::as_tibble(rownames="samples")
      # 
      # eigen_df <- data.frame(
      #   Dimensions = 1:length(pco.counts.norm$eig),
      #   PercentageExplained = (pco.counts.norm$eig/sum(pco.counts.norm$eig))*100 #Percentage explained 
      # )
      # print(pco.counts.norm$li)
      # print(pco.counts.norm$co)
      #Code d'Elise
      # results_PCoA <- ape::pcoa(as.dist(dist.counts.norm))
      # print(results_PCoA)
      # 
      # annot_PCoA <- tibble::tibble(Sample.name = rownames(results_PCoA$vectors)) %>% dplyr::left_join(target, by = c("Sample.name" = "SampleID")) #val is the vector containing the samples
      # results_PCoA <-  compute_arrows(results_PCoA, t(counts.norm))
      # 
      # print("After compute arrows")
      # print(results_PCoA)
      # 
      # DA_results_globale <- Get_pcoa_table(input, dds, group_init,CT,tree)$test 
      # print(DA_results_globale)
      # 
      # variables <-  as.data.frame(0.5*(results_PCoA$U/max(abs(results_PCoA$U[,1:2])))) %>% #scale(center=F)  %>%
      #   tibble::rownames_to_column("Sample.name") %>%
      #   dplyr::mutate(signif="Sample.name" %in% dplyr::filter(DA_results_globale, pvalue_adjusted<0.05)$Id) %>%
      #   dplyr::select(signif, everything())
      # 
      # print(variables)
      # 
      # 
      # to_plot <-  as.data.frame(results_PCoA$vectors) %>% tibble::rownames_to_column("Sample.name") %>% 
      #   dplyr::left_join(annot_PCoA)
      # 
      # print(to_plot)
      # 
      # means <- to_plot %>%  dplyr::group_by("Sample.name")
      #   dplyr::summarise(dplyr::across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
      
      # A1_axis <- paste0("A", as.character(gsub("PC", "", input$PCaxe1))) #A1
      # A2_axis <- paste0("A", as.character(gsub("PC", "", input$PCaxe2))) #A2
      # PC1_int_value <- as.numeric(gsub("PC", "", input$PCaxe1)) #int
      # PC2_int_value <- as.numeric(gsub("PC", "", input$PCaxe2))
      # PCo1_axis <- paste0("PCoA", as.character(gsub("PC", "", input$PCaxe1))) #PCoA1
      # PCo2_axis <- paste0("PCoA", as.character(gsub("PC", "", input$PCaxe2)))
      
    }
      if (plot == "pcoa"){
        # pp <- ggplot(data = pcoa_df, aes_string(x = A1_axis, y = A2_axis)) + 
        #   geom_hline(yintercept = 0, linetype = 2) +
        #   geom_vline(xintercept = 0, linetype = 2) +
        #   geom_point() +
        #   labs(
        #     x = paste0(PCo1_axis, " (", round(eigen_df$PercentageExplained[PC1_int_value], digits = 2), "%)"), 
        #     y = paste0(PCo2_axis, " (", round(eigen_df$PercentageExplained[PC2_int_value], digits = 2), "%)"),
        #     title = "Principal Coordinate Analysis",
        #     subtitle = paste0(input$DistClust, " distance")
        #   ) +          
        #   theme_linedraw() + pcoa_theme 
        
        #print(pcoa_df)
        par(cex=input$cexTitleDiag,mar=c(6,6,4,5))
        # ## Plot axis, label and circles
        v_axes = c(paste0("A", as.character(gsub("PC","",input$PCaxe1))),paste0("A", as.character(gsub("PC","",input$PCaxe2))))        
        pp = plot(pco.counts.norm$li[v_axes],
             xlab=paste(input$PCaxe1, ": ",round(eigen[v_axes[1]],1),"%") ,
             ylab=paste(input$PCaxe2, ": ",round(eigen[v_axes[2]],1),"%"),
             xlim=c(min+0.25*min,max+0.25*max), ylim=c(min-0.1,max+0.1),
             cex.axis=1, cex.lab=1,lwd=2, type="n")
        title(main='Principal Coordinates Analysis ',cex.main=1.5)
        # # pp = ggplot(pcoa_df, aes(x = v_axes[1], y =v_axes[2])) + geom_point() +
        # #   xlab(paste(input$PCaxe1, ": ",round(eigen[v_axes[1]],1),"%")) + 
        # #   ylab(paste(input$PCaxe2, ": ",round(eigen[v_axes[2]],1),"%")) +
        # #   ggtitle("Principal Coordinates Analysis") + pcoa_theme +theme_linedraw()
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
        pp = barplot(eigen[1:nbBar], xlab="Dimensions", ylab="Eigenvalues (%)", names.arg = 1:nbBar, col = col, ylim=c(0,max(eigen)+5), cex.axis=1.2, cex.lab=1.4,cex.names=1.2)
        # eigen_df$ChosenComponents <- ifelse(eigen_df$Dimensions %in% c(as.numeric(gsub("PC", "", input$PCaxe1)), 
        #                                                                as.numeric(gsub("PC", "", input$PCaxe2))), "Chosen", "Not Chosen")
        # 
        # pp <- ggplot(data = eigen_df, aes(x = Dimensions, y = PercentageExplained, fill = ChosenComponents)) +
        #   geom_bar(stat = "identity", width = 0.7) +
        #   labs(
        #     x = "Dimensions",
        #     y = "Percentage Explained Variance",
        #     title = "Scree Plot"
        #   ) +
        #   scale_fill_manual(values = c("Chosen" = "#1D4D68", "Not Chosen" = "lightgray")) +
        #   ylim(0, max(eigen_df$PercentageExplained)) +
        #   geom_text_repel(aes(label = round(PercentageExplained, 2)), vjust = -0.5) +
        #   theme_linedraw() +
        #   theme(
        #     axis.text.x = element_text(angle = 45, hjust = 1),
        #     legend.position = "none"
        #   ) + pcoa_theme +
        #   coord_cartesian(clip = "off")
        
        }
    }
  }
  return(list(plot=pp,dataDiv = dist.counts.norm))
}


### PCA

pca_themes <-  function(input = NULL) {
  
  if(!is.null(input))
  {
  pca_theme <- theme(
    plot.title = element_text(family = "Helvetica", face = "bold", size = 15 *input$cexTitlePlotDiag, hjust = 0.5),
    plot.subtitle = element_text(family = "Helvetica", size = 10 * input$cexSubtitleDiag, hjust = 0.5),
    legend.title = element_text(colour = "black", size =10*input$cexLegendDiag, face = "bold.italic", family = "Helvetica"),
    legend.text = element_text(face = "italic", size =10*input$cexLegendDiag *2/3, colour = "black", family = "Helvetica"),
    axis.title = element_text(family = "Helvetica", size = input$cexAxisLabelDiag *10, colour = "black"),
    axis.text = element_text(family = "Courier", colour = "black", size = input$cexScaleLabelDiag*10)
  )
  }
  else{
    pca_theme <- theme(
      plot.title = element_text(family = "Helvetica", face = "bold", size = (15), hjust = 0.5),
      plot.subtitle = element_text(family = "Helvetica", size = (10), hjust = 0.5),
      legend.title = element_text(colour = "black", face = "bold.italic", family = "Helvetica"),
      legend.text = element_text(face = "italic", colour = "black", family = "Helvetica"),
      axis.title = element_text(family = "Helvetica", size = 10, colour = "black"),
      axis.text = element_text(family = "Courier", colour = "black", size = 10)
    )
  }
  
  return(pca_theme = pca_theme)
}



PCAPlot_meta <-function(input,dds, group_init, n = min(500, nrow(counts(dds))), type.trans = c("VST", "rlog"), 
                        col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen"), plot="pca", pca_res, pca_theme = pca_themes()) 
{
        if(!is.null(pca_res)  && !is.null(input$PCaxe1) && !is.null(input$PCaxe2))
        {
          row.names(pca_res$var$coord) <- substr(rownames(pca_res$var$coord), 1, 20)
          row.names(pca_res$ind$coord) <- substr(rownames(pca_res$ind$coord), 1, 20)
          Groups = pca_res$group
          #PCaxes
          PC1 <- pca_res$ind$coord[, as.numeric(gsub("PC", "", input$PCaxe1))]
          PC2 <- pca_res$ind$coord[, as.numeric(gsub("PC", "", input$PCaxe2))]
          PC1_axis = as.numeric(gsub("PC", "", input$PCaxe1)) #int
          PC2_axis = as.numeric(gsub("PC", "", input$PCaxe2))
          Dim_1 = paste0("Dim.", as.character(PC1_axis))
          Dim_2 = paste0("Dim.", as.character(PC2_axis))
          contribThreshold <- round(length(rownames(pca_res$var$contrib)) - round(input$varSlider*length(rownames(pca_res$var$contrib))))
          if(contribThreshold < 1)
            contribThreshold = contribThreshold +1
          if (plot == "pca") {
            pca_theme = pca_themes(input)
            #Useful variable to center the plot
            x_range <- max(abs(min(PC1)), abs(max(PC1)))
            y_range <- max(abs(min(PC2)), abs(max(PC2)))
            x_lim <- c(as.numeric(-x_range), as.numeric(x_range))
            y_lim <- c(as.numeric(-y_range), as.numeric(y_range))
            var_x_max <- max(abs(pca_res$ind$coord[, min(PC1_axis, PC2_axis)]))
            var_y_max <- max(abs(pca_res$ind$coord[, max(PC1_axis, PC2_axis)]))
            x_range <- max(var_x_max, x_range)
            y_range <- max(y_range, var_y_max)
            # Add some margin to the limits for better visibility
            margin_factor <- 0.1  # You can adjust this factor to suit your needs
            x_margin <- var_x_max * input$cexTitleDiag * margin_factor
            y_margin <- var_y_max * input$cexTitleDiag * margin_factor
            
            # Update the limits with margin
            x_limit <- c(-var_x_max * input$cexTitleDiag * 0.75 - x_margin, var_x_max * input$cexTitleDiag * 0.75 + x_margin)
            y_limit <- c(-var_y_max * input$cexTitleDiag * 0.75 - y_margin, var_y_max * input$cexTitleDiag * 0.75 + y_margin)
            pca_plot <- NULL
            if (input$radioPCA == 1) {
              pca_plot <- fviz_pca_ind(pca_res, axes = c(PC1_axis, PC2_axis), pointsize = input$cexLabelDiag, label = "none", geom = "point") +
                geom_smooth(method = "lm") +
                geom_point(aes(x = PC1, y = PC2, colour = Groups, fill = Groups), size = input$cexLabelDiag) +
                scale_colour_manual(values = col, aesthetics = c("colour", "fill")) +
                coord_cartesian(ylim = c(-y_range * input$cexTitleDiag * 0.75, y_range * input$cexTitleDiag * 0.75)) +
                theme_linedraw() + pca_theme + labs(subtitle = "Euclidean distance")  
            } else if (input$radioPCA == 3) {
              pca_plot <- fviz_pca_var(pca_res, axes = c(PC1_axis, PC2_axis), col.var = "contrib", label = "none", col.circle = "black", select.var = list (contrib = contribThreshold)) +
                labs(color = "Contribution") +
                theme_linedraw() + pca_theme +
                scale_color_gradient(low = "#1565C0", high = "#b92b27") + xlim(c(-1, 1)) +ylim(c(-1, 1))
            }
            if(input$checkCenteringPCA && (input$radioPCA!=3))
            {
              pca_plot <- pca_plot + xlim(x_lim) +ylim(y_lim)
            }
            
            #Checking for the label
            if(input$checkLabelIndiv && (input$radioPCA !=2))
            {
              if (input$radioPCA == 1) {
                pca_plot <- fviz_pca_ind(pca_res, labelsize = input$cexPointsLabelDiag*4, axes = c(PC1_axis, PC2_axis), pointsize = input$cexLabelDiag, repel = TRUE, geom = c("point", "text")) +
                  geom_smooth(method="lm") +
                  geom_point(aes(x = PC1, y = PC2, colour = Groups, fill = Groups), size = input$cexLabelDiag) +
                  geom_label(label.size = input$cexPointsLabelDiag) +
                  geom_text_repel(check_overlap = TRUE) +
                  scale_colour_manual(values = col, aesthetics = c("colour", "fill")) +
                  coord_cartesian(ylim=c(-y_range * input$cexTitleDiag*0.75, y_range * input$cexTitleDiag *0.75)) +
                  theme_linedraw() + pca_theme + labs(subtitle = "Euclidean distance")
              } else if(input$radioPCA == 3)
              {
                pca_plot <- fviz_pca_var(pca_res, labelsize = input$cexPointsLabelDiag*4, axes = c(PC1_axis, PC2_axis), col.var = "contrib", repel = TRUE, geom = c("arrow", "text"), col.circle = "black", select.var = list(contrib = contribThreshold)) +
                  labs(color = "Contribution") +
                  scale_color_gradient(low = "#1565C0", high = "#b92b27") +
                  theme_linedraw() + pca_theme + xlim(c(-1, 1)) +ylim(c(-1, 1))
              }
              if(input$checkCenteringPCA && (input$radioPCA!=3))
              {
                pca_plot <- pca_plot + coord_cartesian(xlim = c(-var_x_max * input$cexTitleDiag*0.75, var_x_max * input$cexTitleDiag *0.75), ylim=c(-var_y_max * input$cexTitleDiag*0.75, var_y_max * input$cexTitleDiag *0.75))
              }
            }
            if(input$radioPCA == 2){
              if(!is.null(input$checkLabelBiplot))
              {
                if(input$checkLabelBiplot[1] == 2)
                {
                  pca_plot <- fviz_pca_biplot(pca_res, labelsize = input$cexPointsLabelDiag*4, axes = c(PC1_axis, PC2_axis), pointsize = input$cexLabelDiag, col.var = "#1d4d68", repel = TRUE, geom.ind = "point", geom.var = c("arrow", "text"), label = "var", select.ind = list(contrib = input$indivSlider), select.var = list(contrib = contribThreshold)) + 
                     geom_smooth(method="lm") +
                     geom_point(aes(x = PC1, y = PC2, colour = Groups, fill = Groups), size = input$cexLabelDiag) +
                     geom_text_repel(check_overlap = TRUE) +
                     scale_colour_manual(values = col, aesthetics = c("colour", "fill")) +
                     coord_cartesian(ylim=c(-y_range * input$cexTitleDiag*0.1, y_range * input$cexTitleDiag*0.1)) +
                     theme_linedraw() + pca_theme + labs(subtitle = "Euclidean distance")
                } else if(input$checkLabelBiplot[1] == 1)
                {
                  pca_plot <- fviz_pca_biplot(pca_res, labelsize = input$cexPointsLabelDiag*4, axes = c(PC1_axis, PC2_axis), pointsize = input$cexLabelDiag, col.var = "#1d4d68", repel = TRUE, geom.ind = c("point", "text"), label = "ind", select.ind = list(contrib = input$indivSlider), select.var = list(contrib = contribThreshold)) +
                    geom_smooth(method="lm") +
                    geom_point(aes(x = PC1, y = PC2, colour = Groups, fill = Groups), size = input$cexLabelDiag) +
                    geom_text_repel(check_overlap = TRUE) +
                    scale_colour_manual(values = col, aesthetics = c("colour", "fill")) +
                    coord_cartesian(ylim=c(-y_range * input$cexTitleDiag*0.1, y_range * input$cexTitleDiag*0.1)) +
                    theme_linedraw() + pca_theme +labs(subtitle = "Euclidean distance") 
                }
                if(length(input$checkLabelBiplot) == 2) #if both of label buttons are pressed, the input vector has a length of 2
                {
                  pca_plot <- fviz_pca_biplot(pca_res, labelsize = input$cexPointsLabelDiag*4, axes = c(PC1_axis, PC2_axis), pointsize = input$cexLabelDiag, col.var = "#1d4d68", repel = TRUE, geom.ind = c("point", "text"), geom.var = c("arrow", "text"), select.ind = list(contrib = input$indivSlider), select.var = list(contrib = contribThreshold)) +
                    geom_smooth(method="lm") +
                    geom_point(aes(x = PC1, y = PC2, colour = Groups, fill = Groups), size = input$cexLabelDiag) +
                    geom_text_repel(check_overlap = TRUE) +
                    scale_colour_manual(values = col, aesthetics = c("colour", "fill")) +
                    coord_cartesian(ylim=c(-y_range * input$cexTitleDiag*0.1, y_range * input$cexTitleDiag*0.1)) +
                    theme_linedraw() + pca_theme +labs(subtitle = "Euclidean distance") 
                }
              }
              if(is.null(input$checkLabelBiplot)) #None of the label buttons are pressed
              {
                pca_plot <- fviz_pca_biplot(pca_res, labelsize = input$cexPointsLabelDiag*4, axes = c(PC1_axis, PC2_axis), pointsize = input$cexLabelDiag, col.var = "#1d4d68", repel = TRUE, geom.ind = "point", geom.var = "arrow", label = "none", select.ind = list(contrib = input$indivSlider), select.var = list(contrib = contribThreshold)) +
                  geom_point(aes(x = PC1, y = PC2, colour = Groups, fill = Groups), size = input$cexLabelDiag) +
                  geom_text(check_overlap = TRUE) +
                  scale_colour_manual(values = col, aesthetics = c("colour", "fill")) +
                  geom_smooth(method="lm") +
                  coord_cartesian(ylim=c(-y_range * input$cexTitleDiag*0.1, y_range * input$cexTitleDiag*0.1)) +
                  theme_linedraw() + pca_theme +labs(subtitle = "Euclidean distance")
              }
              if(input$checkCenteringPCA)
              {
                pca_plot <- pca_plot + coord_cartesian(xlim = x_limit, ylim = y_limit)              }
            }
            return(pca_plot)
          }
          
          if(plot == "eigen")
          {
            #Find the chosen PCA axis in order to color the chosen bars in blue
            pca_res$eigen_df$ChosenComponents <- ifelse(pca_res$eigen_df$Dimensions %in% c(as.numeric(gsub("PC", "", input$PCaxe1)), 
                                                                                           as.numeric(gsub("PC", "", input$PCaxe2))), "Chosen", "Not Chosen")
            eigen_plot <- ggplot(data = pca_res$eigen_df, aes(x = Dimensions, y = PercentageExplained, fill = ChosenComponents)) +
                       geom_bar(stat = "identity", width = 0.7) +
                       labs(
                         x = "Dimensions",
                         y = "Percentage Explained Variance",
                         title = "Scree Plot"
                       ) +
                       scale_fill_manual(values = c("Chosen" = "#1D4D68", "Not Chosen" = "lightgray")) +
                       ylim(0, max(pca_res$eigen_df$PercentageExplained)) +
                       geom_text_repel(aes(label = round(PercentageExplained, 2)), vjust = -0.5) +
                       theme_linedraw() +
                       theme(
                         axis.text.x = element_text(angle = 45, hjust = 1),
                         legend.position = "none"
                       ) +
                       pca_theme +
                       coord_cartesian(clip = "off")
            return(eigen_plot)
          }
          
          if(plot == "contrib")
          {
            Dim_1 = paste0("Dim.", as.character(PC1_axis))
            Dim_2 = paste0("Dim.", as.character(PC2_axis))
            min_axis <- min(PC1_axis, PC2_axis)
            
            contrib_df_filtered <- pca_res$contrib_df %>%
              dplyr::filter(Dimension %in% paste0("Dim.", c(PC1_axis, PC2_axis))) %>%
              dplyr::arrange(Variable, Dimension) %>%
              tidyr::pivot_wider(names_from = Dimension, values_from = Contribution) %>%
              dplyr::arrange(desc(!!sym(paste0("Dim.", min_axis)))) %>%
              dplyr::slice(1:contribThreshold)
            

            summed_contrib_df <- contrib_df_filtered %>%
              dplyr::mutate(summed_dimensions = rowSums(contrib_df_filtered[, c(Dim_1, Dim_2)])) %>%
              dplyr::select(Variable, summed_dimensions)
            if(input$sumContrib)
            {
              melted_df <- reshape2::melt(summed_contrib_df, id.vars = "Variable")
              
              # Plot
              contrib_plot <- ggplot(melted_df, aes(x = reorder(Variable, -value), y = value, fill = variable)) +
                geom_bar(stat = "identity", position = "stack", alpha = 0.3) +
                labs(title = "Sum of contributions",
                     y = "Total contribution",
                     x = input$TaxoSelect) +
                theme_linedraw() +
                pca_theme +
                theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
                        legend.title = element_blank()) +
                scale_fill_manual(values = "#9f1853") + 
                guides(fill = FALSE)
              
            }
            else{
            contrib_plot <- ggplot(contrib_df_filtered, aes(x = reorder(Variable, -!!sym(Dim_1)), y = !!sym(Dim_1), fill = "Dimension 1")) +
              geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.3) +
              geom_bar(aes(y = !!sym(Dim_2), fill = "Dimension 2"),
                       stat = "identity", position = position_dodge(width = 0.8), alpha = 0.3) +
              labs(title = paste0("Contribution to Dimensions ", as.character(PC1_axis), " and ", as.character(PC2_axis)), x = input$TaxoSelect, y = "Contribution") +
              theme_linedraw() +
              pca_theme +
              theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
                    legend.title = element_blank()) +
              scale_fill_manual(values = c("Dimension 1" = "#1192e8", "Dimension 2" = "#198038"),
                                name = "Dimension",
                                breaks = c("Dimension 1", "Dimension 2"),
                                labels = c(paste("Dimension", as.character(PC1_axis)), paste("Dimension", as.character(PC2_axis)))) +
              guides(fill = guide_legend(title = NULL))
            }
            
            # contrib_plot <- ggplot(contrib_df_filtered, aes(x = Dimension, y = Contribution, fill = Variable, width = CumulativeProportion)) +
            #   geom_bar(stat = "identity") +
            #   coord_flip() +
            #   theme_minimal() +
            #   labs(x = "Dimension", y = "Contribution") +
            #   scale_x_discrete(expand = c(0.1, 0)) + theme_linedraw() + pca_theme
            
          }
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
  dist.counts.norm = NULL
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
  
  dist.counts.trans = NULL
  
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
    
    ## To get always the same result 
    set.seed(666)
    #permanova_test = adonis(dist.counts.norm~group)
    type.trans <- type.trans[1]
    
    if (type.trans == "VST") counts.trans <- SummarizedExperiment::assay(varianceStabilizingTransformation(dds))
    else counts.trans <- SummarizedExperiment::assay(rlogTransformation(dds))
    counts.trans = counts.trans[,ind_kept]
    
    rv = apply(counts.trans, 1, var, na.rm = TRUE)
    dat <- t(counts.trans[order(rv, decreasing = TRUE),][1:n, ]) %>% data.frame
    pca <- FactoMineR::PCA(dat, ncp = 12, scale.unit = TRUE, graph = FALSE)
    dist.counts.trans <- factoextra::get_dist(dat, method = "euclidean")
    #pca = prcomp(t(counts.trans[order(rv, decreasing = TRUE),][1:n, ]))
    permanova_test = vegan::adonis(dist.counts.trans~group)
    pcaTable <- pca$ind$coord
    return(list(table = pcaTable, test = permanova_test, dist_val = dist.counts.trans))
  }
  else
    return(list(table = NULL, test = NULL))
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
      #print(proj[, 1:ncol(proj)])
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


