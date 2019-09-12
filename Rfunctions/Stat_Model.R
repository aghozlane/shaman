#@ This file contains all the functions for the 
#@ statistical modelling (dds object, interaction, ...),
#@ and the functions concerning the contrasts



##############################################################
##
##          Statistical modelling
##
##############################################################


## Get the possible interaction for the statistical model
GetInteraction2 <- function(target,VarInt)
{ 
  res=c()
  
  if(length(VarInt)>=2)
  {
    k=1
    for(i in 1:(length(VarInt)-1))
    { 
      for(j in (i+1):length(VarInt))
      { 
        res[k] = paste(VarInt[i],":",VarInt[j],sep="")
        k = k+1
      }
    }
  }
  return(res)
}


## Get the dds object of DESeq2
Get_dds_object <- function(input,counts,target,design,normFactorsOTU,CT_noNorm,CT_Norm,tree=NULL)
{
  dds <- DESeq2shaman::DESeqDataSetFromMatrix(countData=counts, colData=target, design=design,ignoreRank = TRUE)
  sizeFactors(dds) = normFactorsOTU
  ddstmp = NULL
  
#   dds <- estimateDispersions(dds, fitType=input$fitType)
#   if(as.numeric(R.Version()$major)>=3 && as.numeric(R.Version()$minor) >=1.3){
#     dds <- nbinomWaldTest(dds)
#   }else{
#     dds <- nbinomWaldTest(dds,modelMatrixType = "expanded")
#   }
#   countsNorm = counts(dds, normalized = TRUE)

  # if(as.numeric(R.Version()$major)>=3 && as.numeric(R.Version()$minor) >=1.3){
  #   if(nrow(counts)*nrow(target)>=50000)
  #   {
  #     try(DESeq(dds,fitType=input$fitType,parallel = TRUE,minReplicatesForReplace = Inf) ->ddstmp,silent = TRUE)
  #     if(!is.null(ddstmp)) dds = ddstmp
  #   }
  #   if(nrow(counts)*nrow(target)<50000 || is.null(ddstmp))
  #   {
  #     dds <- estimateDispersions(dds, fitType=input$fitType)
  #     dds <- nbinomWaldTest(dds)
  #   }
  # }else{
    # dds <- estimateDispersions(dds, fitType=input$fitType)
    # dds <- nbinomWaldTest(dds,modelMatrixType = "expanded")
   if(nrow(counts)*nrow(target)>=50000)
    {
     dds = DESeq2shaman::DESeq(dds,fitType=input$fitType,parallel = TRUE,minReplicatesForReplace = Inf)
        # try(DESeq(dds,fitType=input$fitType,parallel = TRUE,minReplicatesForReplace = Inf) ->ddstmp,silent = TRUE)
        # if(!is.null(ddstmp)) dds = ddstmp
    }
   if(nrow(counts)*nrow(target)<50000 || is.null(dds))
    {
        dds <- DESeq2shaman::estimateDispersions(dds, fitType=input$fitType)
        dds <- DESeq2shaman::nbinomWaldTest(dds)
    }
    # dds <-  DESeq(dds,fitType=input$fitType,modelMatrixType = "expanded",parallel = TRUE)
  # }
  countsNorm = DESeq2shaman::counts(dds, normalized = TRUE)
  
  
  #save(dds,file="dds.RData")
  return(list(dds = dds,raw_counts=counts,countsNorm=countsNorm,target=target,design=design,normFactors = normFactorsOTU,CT_noNorm=CT_noNorm,CT_Norm=CT_Norm))
}


## Get the design according to the input
GetDesign <- function(input,target)
{
  design = NULL
  Interaction = NULL
  InterVar = input$InterestVar
  if(!any(!InterVar%in%colnames(target)))
  {
    if(length(InterVar)>1) Interaction = input$Interaction2
    alltmp = c(InterVar,Interaction)
    if(length(alltmp)>0) design = as.formula(paste("~", paste0(alltmp, collapse= "+")))
  }
  return(design)
}



##############################################################
##
##          Contrasts
##
##############################################################

## Create the contrast base
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

## Create the contrast base with the buttons
BaseContrastEasy <- function(input,names,namesfile,target)
{  
  
  v_tmp = rep(0,length(names))
  filesize = file.info(namesfile)[,"size"]
  F1 = NULL
  nameContrast = ""
  
  ## Get the selected modalities
  M1 = input$Select1_contrast
  M2 = input$Select2_contrast

  if(length(input$Interaction2)>0) F1 = input$Select3_contrast
  ## Get the name of the parameter corresponding to the modalities
  InterVar = input$InterestVar
  if(length(InterVar)>1) {Sel_Var = InterVar[which(unlist(lapply(target[,InterVar],FUN = function(x){M1%in%x})))]} else Sel_Var = InterVar
  names1dds = paste(Sel_Var,M1,sep="")
  names2dds = paste(Sel_Var,M2,sep="")

  ## fill the vector
  ind1 = which(names%in%names1dds)
  ind2 = which(names%in%names2dds)
  
  if(length(ind1)>0) v_tmp[ind1] = 1
  if(length(ind2)>0) v_tmp[ind2] = -1
  
  
  nameContrast = paste(M1,"_vs_",M2,sep="")
  
  if(F1!="All" && !is.null(F1)){
    Sel_Var_For = InterVar[which(unlist(lapply(target[,InterVar],FUN = function(x){F1%in%x})))]
    ## Depends on the interation
    namesfor1 = paste(Sel_Var,M1,".",Sel_Var_For,F1,sep="")
    namesfor1.1 = paste(Sel_Var_For,F1,".",Sel_Var,M1,sep="")
    namesfor2 = paste(Sel_Var,M2,".",Sel_Var_For,F1,sep="")
    namesfor2.1 = paste(Sel_Var_For,F1,".",Sel_Var,M2,sep="")
    ind1.for = which(names%in%c(namesfor1,namesfor1.1))
    ind2.for = which(names%in%c(namesfor2,namesfor2.1))
    if(length(ind1.for)>0) v_tmp[ind1.for] = 1
    if(length(ind2.for)>0) v_tmp[ind2.for] = -1
    nameContrast = paste(nameContrast,"_for_",F1,sep="")
  }
  
  if(filesize!=0)
  { 
    oldContrast = read.table(namesfile,header=TRUE)
    colnamesTmp = c(colnames(oldContrast),nameContrast)
    mat = cbind(oldContrast,v_tmp)
  }
  else{ colnamesTmp = nameContrast; mat = v_tmp}
  
  write.table(mat,namesfile,row.names=FALSE,col.names = colnamesTmp)
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

## Check the rank of the matrix for the statistical model.
CheckMatrixRank <- function(design,target)
{
  testRank = FALSE

  if(!is.null(design) && length(as.character(design))) 
  {
    tmp = strsplit(as.character(design)[2],"[ ]*[+][ ]*")
    testDesign = any(unlist(tmp)%in%colnames(target))
  }
  
  if(nrow(target)>=1 && ncol(target)>=2 && nrow(target)>1 && !is.null(design) && testDesign){
    modelMatrix <- model.matrix(design, data = as.data.frame(target))
    rk = qr(modelMatrix)$rank
    testRank = (ncol(modelMatrix)==rk)
  }
    return(testRank)
}


