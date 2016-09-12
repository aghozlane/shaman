#@ This file contains all the functions needed to
#@ to load, check and transform the data

## Add news to the home page
addNews <- function(date ="",title="",text="")
{
  res=list()
  res$r1 = paste("<b><font size='+1'>",date,"</font></b>", " - ", "<b><font size='+1'>",title,"</font></b><br/>")
  res$r2 = paste("<p><font color='grey'>",text,"</font></p><hr/>")
  
  return(HTML(unlist(res)))
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





## Check the format of the counts table
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

## Check the format of the taxonomy table
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

## Get the percentage of annotated OTU
PercentAnnot <- function(counts,taxo)
{
  Error = NULL  
  tmp = table(rownames(counts)%in%rownames(taxo))
  Percent = tmp["TRUE"]/sum(tmp)
  if(is.na(Percent)) Percent = 0
  if(Percent==0){Error = "Counts table and annotation do not matched" }
  
  return(list(Error=Error,Percent=Percent))
}


## Get the counts and the taxo tables from the BIOM format file.
GetDataFromBIOM <-function(dataBIOM)
{
  ## Counts table
  counts = biom_data(dataBIOM)
  counts = as.matrix(counts)
  ## Change of - to . is risky
  colnames(counts) = gsub("-",".",colnames(counts))
  counts = as.data.frame(counts)
  CheckCounts = CheckCountsTable(counts)
  counts = CheckCounts$counts
  
  ## Taxonomy table
  taxo = as.data.frame(observation_metadata(dataBIOM))
  OTUnames = rownames(taxo)
  ## Modif taxo table (remove p__,... and change the colnames)
  taxo = as.data.frame(sapply(taxo,gsub,pattern="^.*__",replacement=""))
  colnames(taxo) = c("Kingdom", "Phylum","Class","Order","Family","Genus","Species")
  rownames(taxo) = OTUnames
  ## Remove empty row
  taxo[taxo==""] = NA
  taxo[taxo=="Unassigned"] = NA
  taxo=taxo[rowSums(is.na(taxo))!=dim(taxo)[2], ]
  
  CheckTaxo = CheckTaxoTable(taxo,counts)
  
  ## Pourcentage of annotation
  tmp = PercentAnnot(counts,taxo)
  
  return(list(counts=counts,taxo=taxo,CheckCounts=CheckCounts,CheckTaxo=CheckTaxo,Percent=tmp$Percent,CheckPercent=tmp$Error))
}

## Check the data
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
  
  ## Get the normalization variable (normalization can be done according to this variable)
  VarNorm = input$VarNorm
  
  if(length(ind)==length(labels))
  { 
    if(input$TypeTable == "MGS"){
      ## Get the feature size for the normalisation
      Size_indcol = which(toupper(colnames(CT))%in%"SIZE")
      if(length(Size_indcol)==1) FeatureSize = CT[,Size_indcol]
      else print("Size parameter is missing in the count matrix")
      # Consider only counts
      CT = CT[,ind]
      # Divide by gene length
      CT = CT / FeatureSize * 1000
      # Convert matrix as integer
      CT_int=t(apply(CT,1,as.integer))
      rownames(CT_int)=rownames(CT)
      colnames(CT_int)=colnames(CT)
      CT=CT_int
    } else CT = CT[,ind]
    
    save.image("test.RData")
    ## Order CT according to the target
    CT = OrderCounts(counts=CT,labels=labels)$CountsOrder
    CT_noNorm = CT
    RowProd = sum(apply(CT_noNorm,1,prod))
    
    ## Create the dds object
    dds <- DESeqDataSetFromMatrix(countData=CT, colData=target, design=design)
    
    if(is.null(VarNorm)){
      ## Counts normalisation
      ## Normalisation with or without 0
      if(input$AccountForNA || RowProd==0) dds = estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)),geoMeans=GeoMeansCT(CT))
      if(!input$AccountForNA && RowProd!=0) dds = estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)))
      normFactors = sizeFactors(dds)
      
    } else{
      group = as.data.frame(target[,VarNorm])
      group = apply(group,1,paste, collapse = "-")
      normFactors = c()
      mod = unique(group)
      ## At least 2 samples are needed for the normalization
      if(min(table(group))>1){
        for(i in unique(group))
        {
          indgrp = which(group==i) 
          CT_tmp = CT[,indgrp]
          CT_tmp = removeNulCounts(CT_tmp) 
          target_tmp = data.frame(labels = rownames(target)[indgrp])
          dds_tmp <- DESeqDataSetFromMatrix(countData=CT_tmp, colData=target_tmp, design=~labels)
          if(input$AccountForNA) dds_tmp = estimateSizeFactors(dds_tmp,locfunc=eval(as.name(input$locfunc)),geoMeans=GeoMeansCT(CT_tmp))
          if(!input$AccountForNA) dds_tmp = estimateSizeFactors(dds_tmp,locfunc=eval(as.name(input$locfunc)))
          normFactors[indgrp] = sizeFactors(dds_tmp)
        }
      } else{
        if(input$AccountForNA || RowProd==0) dds = estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)),geoMeans=GeoMeansCT(CT))
        if(!input$AccountForNA && RowProd!=0) dds = estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)))
        normFactors = sizeFactors(dds)
      }
      
      sizeFactors(dds) = normFactors
    }
    
    ## Keep normalized OTU table
    CT_Norm = counts(dds, normalized=TRUE)
    save(CT_Norm,dds,CT,taxo,taxoSelect,file="test.RData")
    # Only interesting OTU
    # merged_table = merge(CT, taxo[order(rownames(CT)),], by="row.names")
    merged_table = merge(CT, taxo, by="row.names")
    CT = merged_table[,2: (dim(CT)[2]+1)]
    taxo = merged_table[,(dim(CT)[2]+2):dim(merged_table)[2]]
    rownames(CT) = merged_table[,1]
    rownames(taxo) = merged_table[,1]
    #ordOTU = order(rownames(taxo))
    counts_annot = CT
    #       ordOTU = order(rownames(taxo))
    #       indOTU_annot = which(rownames(CT)%in%rownames(taxo))
    #       counts_annot = CT[indOTU_annot[ordOTU],]
    ## Aggregate matrix
    if(taxoSelect=="OTU/Gene") counts = counts_annot
    else{
      if(input$TypeTable == "MGS"){
        taxoS = taxo[,input$TypeTable]
        counts = aggregate(counts_annot,by=list(Taxonomy = taxoS),mean)
        rownames(counts)=counts[,1]
        counts=counts[,-1]
        counts_int=t(apply(counts,1,as.integer))
        rownames(counts_int)=rownames(counts)
        colnames(counts_int)=colnames(counts)
        counts=counts_int
      }
      if(taxoSelect != "MGS"){
        #taxoS = taxo[ordOTU,taxoSelect]
        taxoS = taxo[,taxoSelect]
        counts = aggregate(counts_annot,by=list(Taxonomy = taxoS),sum)
        rownames(counts)=counts[,1];counts=counts[,-1]
      }
    }
    
    ## Ordering the counts table according to the target labels 
    tmpOrder = OrderCounts(counts,normFactors,labels)
    counts = tmpOrder$CountsOrder
    normFactors = tmpOrder$normFactorsOrder
    CheckTarget = TRUE
  }
  return(list(counts=counts,CheckTarget=CheckTarget,normFactors=normFactors, CT_noNorm=CT_noNorm, CT_Norm =CT_Norm))
  #return(list(counts=counts,target=target[ind,],labeled=labeled,normFactors=normFactors, CT_noNorm=CT_noNorm))
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

