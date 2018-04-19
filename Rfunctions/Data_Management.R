#@ This file contains all the functions needed to
#@ to load, check, filter and transform the data 

create_forked_task <- function(expr) {
  makeReactiveBinding("state")
  state <- factor("running",
                  levels = c("running", "success", "error", "cancel"),
                  ordered = TRUE
  )
  
  result <- NULL
  
  # Launch the task in a forked process. This always returns
  # immediately, and we get back a handle we can use to monitor
  # or kill the job.
  task_handle <- parallel::mcparallel({
    force(expr)
  })
  
  # Poll every 100 milliseconds until the job completes
  o <- observe({
    res <- parallel::mccollect(task_handle, wait = FALSE)
    if (is.null(res)) {
      invalidateLater(100)
    } else {
      o$destroy()
      if (!is.list(res) || length(res) != 1 || !inherits(res[[1]], "try-error")) {
        state <<- "success"
        result <<- res[[1]]
      } else {
        state <<- "error"
        result <<- attr(res[[1]], "condition", exact = TRUE)
      }
    }
  })
  
  list(
    completed = function() {
      state != "running"
    },
    result = function() {
      if (state == "running") {
        # If running, abort the current context silently.
        # We've taken a reactive dependency on "state" so if
        # the state changes the context will invalidate.
        req(FALSE)
      } else if (state == "success") {
        return(result)
      } else if (state == "error") {
        stop(result)
      } else if (state == "cancel") {
        validate(need(FALSE, "The operation was cancelled"))
      }
    },
    cancel = function() {
      if (state == "running") {
        state <<- "cancel"
        o$destroy()
        tools::pskill(task_handle$pid, tools::SIGTERM)
        tools::pskill(-task_handle$pid, tools::SIGTERM)
        parallel::mccollect(task_handle, wait = FALSE)
      }
    }
  )
}


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
CheckCountsTable <- function(counts, MGSTable=FALSE)
{
  Error = NULL
  Warning = NULL
  
  if(is.null(counts) && is.null(Error)){Error = "There is no counts table" }
  

  if(ncol(counts)<=1 && is.null(Error)){Error = "The number of columns of the counts table must be at least 2" }
  if(nrow(counts)<=1 && is.null(Error)){Error = "The number of rows of the counts table must be at least 2" }
  
  if(is.null(Error)) 
  {
    numTest = FALSE%in%sapply(counts,is.numeric)
    if(numTest) Error = "The counts table must contain only numeric values" 
    if(!numTest)
    {
      if(0%in%colSums(counts)){Error = "At least one of the column of the counts table is 0" }
      if(min(counts)<0){Error = "The counts table must contain only positive values" }
      if(MGSTable && length(which(toupper(colnames(counts))%in%"SIZE")) != 1){Error="The counts table must contain a column named SIZE providing the length of each gene"}
    }
  }
  if(TRUE%in%sapply(counts,is.na) && is.null(Error)){Warning = "NA values are considered as 0 is the counts table"; counts[sapply(counts,is.na)]=0}
  
  return(list(Error=Error,Warning=Warning,counts=counts))
}

## Check the format of the taxonomy table
CheckTaxoTable <- function(taxo,counts, MGSTable=FALSE, taxoCreated=FALSE)
{
  Error = NULL
  Warning = NULL
  if(taxoCreated){Warning = "No taxonomy table has been uploaded, the analysis can only be done at the OTU/gene level"}
  if(ncol(taxo)<1 && is.null(Error)){Error = "The number of columns of the taxonomy table must be at least 1" }
  else if(nrow(taxo)<=1 && is.null(Error)){Error = "The number of rows if the taxonomy table must be at least 2" }
  if(TRUE%in%is.numeric(taxo) && is.null(Error) ){Error = "The taxonomy table must contain only character" }
  
  if(is.null(Error))
  {
    for(i in 1:ncol(taxo))
    {
      level = levels(taxo[,i])
      nb = length(level)
      if(nb==1 && level=="NA"){ Error = "At least one column contains only NA"}
    }
    if(MGSTable && length(which(toupper(colnames(taxo))%in%"MGS")) != 1){
      Error="The taxonomy table must contain a column named MGS providing the MGS association of each gene"
      }
  }
  
  ## Annotated features without counts
  if(!any(rownames(taxo)%in%rownames(counts)) && is.null(Error)){ Error = "Some annotated features are not in the count table"}
  
  return(list(Error=Error,Warning=Warning))
}



CheckTargetModel <- function(input,target,labeled,CT)
{
  Error = NULL
  HowTo = NULL
  InterVar = input$InterestVar
  labels = rownames(target)
  ind = which(colnames(CT)%in%labels)
#   InterVar%in%
#   uniq_column = (length(which(sapply(target[InterVar], function(x) length(unique(x))) == 1)) > 0)
#   uniq_column_names = names(which(sapply(target[InterVar], function(x) length(unique(x))) == 1))
  
  ## At least one variable selected
  if(is.null(Error) && length(ind)<=1){
    Error = "Less than two samples names fit with the counts table" 
    HowTo = "Check the samples names in the target file. They must be in the first column and must correspond EXACTLY to the names in the count table."
  }
  ## At least one variable selected
  if(is.null(Error) && length(InterVar)==0){
    Error = "At least one variable must be selected for the model" 
    HowTo = "Add at least one variable in the 'Select the variables' widget"
  }
  
  ## Names of samples correct ?
  if(is.null(Error) && labeled==0){
    Error = "The names of the samples in the target file do not fit the counts table" 
    HowTo = "Check the samples names in the target file. They must be in the first column and must correspond EXACTLY to the names in the count table."
    }
  
  ## Number of columns
  if(is.null(Error) && ncol(target)<2){
    Error = "The number of columns of the target file must be at least 2"
    HowTo = "Add at least one additional variable to describe your samples"
  }
  
  if(is.null(Error) && min(sapply(apply(target,2,unique),length)) <=1){
    Error = "One of the variable has the same value for all the samples" 
    HowTo = "Remove the variable from the target file"
  }
  
  
  ## Number of rows
  if(is.null(Error) && nrow(target)<=1){
    Error = "The number of rows if the target file must be at least 2"
    HowTo = "Add information about more than 2 samples"
    }
  
  ## NA values
  if(is.null(Error) && (any(is.na(target)) || any(target ==""))){
    Error = "NA's or missing values are not allowed in the target file" 
    HowTo = "Remove all the samples for which one or more variables are NA or missing"
  }
  
  ## contrasts can be applied only to factors with 2 or more levels
  
#   if(is.null(Error) && (uniq_column)){
#     Error = "Contrasts can be applied only to factors with 2 or more levels."
#     HowTo = paste("Remove all variables with only one factor:", uniq_column_names, sep=" ")
#   }
#   
  
  ## Full rank matrix
  if(is.null(Error) && length(InterVar)>0)
  {
    design = GetDesign(input,target)
    testRank = CheckMatrixRank(design,target)
    if(!testRank){
        Error = "The model matrix is not full rank. One or more variables or interaction terms 
        are linear combinations of the others and must be removed." 
        HowTo = "Remove variable(s) that provide the same information, i.e, if the value of a variable is totaly determine by an other variable remove one of them."
        }
  }
    
  return(list(Error=Error,HowTo=HowTo))
}




CheckContrast <- function(contrastFile,dds)
{
  Error = NULL
  Warning = NULL
  parameterNames = resultsNames(dds)
  if(is.null(contrastFile) && is.null(Error)){Error = "The format of the contrast file is not supported by SHAMAN" }
  
  
  if(ncol(contrastFile)<1 && is.null(Error)){Error = "The contrast file seems to be empty" }
  if(nrow(contrastFile)!=length(parameterNames) && is.null(Error)){Error = "The contrast file does not fit with the model parameters" }

  if(TRUE%in%sapply(contrastFile,is.na) && is.null(Error)){Error = "NA values are considered as 0 is the counts table"; contrastFile[sapply(contrastFile,is.na)]=0}
  
  
  return(list(Error=Error,Warning=Warning,contrastFile=contrastFile))
}


## Check the format of the tree file (for Unifrac distance)
CheckTreeFile <- function(tree)
{
  Error = NULL
  Warning = NULL
  if(!is.phylo(tree) && is.null(Error)){Error = "The loaded file is not a phylogenetic tree"; tree = NULL}
  if(!is.rooted(tree) && is.null(Error) ){
    Warning = "The tree has been rooted using midpoint method";
    roottree = try(midpoint.root(tree), TRUE)
    if (class(roottree) == "try-error"){
      D <- cophenetic(tree)
      dd <- max(D)
      ii <- which(D == dd)[1]
      ii <- c(ceiling(ii/nrow(D)), ii%%nrow(D))
      if (ii[2] == 0) ii[2] <- nrow(D)
      spp <- rownames(D)[ii]
      nn <- which(tree$tip.label == spp[2])
      tree <- reroot(tree, nn, tree$edge.length[which(tree$edge[,2] == nn)])
    } 
    else tree=roottree
  }
  return(list(Error=Error,Warning=Warning,tree=tree))
}



## Check Masque Input
CheckMasque <- function(input, values, check_mail=FALSE)
{
  Error = NULL
  HowTo = NULL
 
  ## Check password
  
  # if(is.null(Error) && input$password == ""){
  #   Error = "<h6><strong>Empty key field </strong></h6>"
  #   HowTo = "<h6><strong>Make sure that you have click the &laquo Get key &raquo button and that you have pasted the key sent by mail </strong></h6>"
  # }
  # 
  if(is.null(Error) && is.null(values$login_email) && check_mail){
       Error = "<h6><strong>Invalid key </strong></h6>";
       HowTo = "<h6><strong>Make sure that you have click the &laquo Get key &raquo button </strong></h6>"
  }

  ## At least one fastq is detected
  # if(is.null(Error) && input$LoadFiles>0 && length(values$fastq_names_only)==0){
  #   Error = "<h6><strong>The selected directory must contain at least one file in the following format : fastq, fastq.gz, or fq.</strong></h6>" 
  #   HowTo = "<h6><strong>Change the working directory and check the format of the files</strong></h6>"
  # }
  
  if(is.null(Error) && input$PairedOrNot=='y' && input$MatchFiles_button>0){
    if(length(values$R2fastQ) !=length(values$R2fastQ)){
        Error = "<h6><strong>The number of fastq files for R1 and R2 must be the same</strong></h6>" 
        HowTo = "<h6><strong>Add/Remove some files or change the suffix to identify the pairs</strong></h6>"
    }
    
    if(length(values$R2fastQ)>0  && length(values$R2fastQ)>0){
      tmpR1 = gsub(input$R1files,x=values$R1fastQ,"") 
      tmpR2 = gsub(input$R2files,x=values$R2fastQ,"")
      
      dup_files = c(values$R1fastQ[duplicated(tmpR1)],values$R2fastQ[duplicated(tmpR2)])
      if(length(dup_files)>0){
        Error = paste("<h6><strong>These fastq files corresponds to the same sample names:</strong></h6>" ,dup_files)
        HowTo = "<h6><strong>Change the suffix to identify the pairs</strong></h6>"
      }
      
      if(!isValidPrimer(input$R1primer)){ Error = "<h6><strong>The primer (forward) must only contain letters from A to Z</strong></h6>" }
      if(!isValidPrimer(input$R2primer)){ Error = "<h6><strong>The primer (reverse) must only contain letters from A to Z</strong></h6>" }
      
    }
    
  }
  if(is.null(Error) && !isValidEmail(input$to)) Error = "<h6><strong>The email address is not valid</strong></h6>"
  

  
  
  if(is.null(Error) && !isValidPrimer(input$primerSingle)){ Error = "<h6><strong>The primer must only contain letters from A to Z</strong></h6>" }
  
  if(is.null(Error)) {
    
    res = SamplesMasque(input,values)
    if(length(res$samples)==0) {
      Error = "<h6><strong>0 sample detected</strong></h6>"
      if(input$PairedOrNot=='y') HowTo = '<h6><strong>Make sure that you click the &laquo Match &raquo button. <br /> Change the working directory and/or verify the pairs matching.</strong></h6>'
      if(input$PairedOrNot=='n') HowTo = '<h6><strong>Make sure that your samples have the correct extension (.fastq, .fq, .fastq.gz or .fq.gz). <br /> Change the working directory.</strong></h6>'
    }
    
  }
  
  if(is.null(Error)) {
    error_file = paste(values$curdir,"www","masque","error",paste('file',values$masque_key,"_error.txt",sep=""),sep= .Platform$file.sep)
    if(file.exists(error_file)){
      Error = "<h6><strong>An error happened during the workflow progress. Please check your email.</strong></h6>" 
    }
  }
  
  return(list(Error=Error,HowTo=HowTo))
}



isValidEmail <- function(x) {
  grepl("\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>", as.character(x), ignore.case=TRUE)
}

isValidPrimer <- function(x) {
  !any(!grepl("[A-Z]",unlist(strsplit(x,""))))
}


SamplesMasque <- function(input,values)
{
  samples_removed = NULL
  samples = NULL
  
  if(input$PairedOrNot=='y')
  {
    tmpR1 = gsub(input$R1files,x=values$R1fastQ,"") 
    tmpR2 = gsub(input$R2files,x=values$R2fastQ,"")
    
    
    R1samples = tmpR1[tmpR1%in%tmpR2]; R1samples_removed = values$R1fastQ[!tmpR1%in%tmpR2]
    R2samples = tmpR2[tmpR2%in%tmpR1]; R2samples_removed = values$R2fastQ[!tmpR2%in%tmpR1]
    
    samples = unique(c(R1samples,R2samples))
    samples_removed = c(R1samples_removed,R2samples_removed)
  } else {samples = unique(values$fastq_names_only)}
  
  return(list(samples=samples,samples_removed=samples_removed))
}



CreateJSON <- function(input,values){

  tmp = tempdir()
  path_fasta = paste(tmp,paste(basename(file_path_sans_ext(values$json_name)),"_contaminant.fasta",sep=""),sep = .Platform$file.sep)
  
  if(input$PairedOrNot=='n')
  {
    path_fastq = paste(tmp,"Masque_files",sep= .Platform$file.sep)
    
    df = data.frame("paired"=FALSE,
                    "path"=path_fastq,
                    "host"=input$HostName,
                    "type"=input$DataTypeMasque,
                    "mail"=values$login_email,
                    "contaminant"= path_fasta,
                    "phredthres" = input$phredthres,
                    "mincorrect" = input$mincorrect,
                    "minreadlength" = input$minreadlength,
                    "dreptype" = input$dreptype,
                    "maxampliconlength" = input$maxampliconlength,
                    "minampliconlength" = input$minampliconlength,
                    "minabundance" = input$minabundance,
                    "clusteringthreshold" = input$clusteringthreshold,
                    "clusteringstrand" = input$clusteringstrand,
                    "annotationstrand" = input$annotationstrand,
                    "aKmin" = input$annotationKingdomthreshold,
                    "aPmin" = input$annotationPhylumthreshold[1],
                    "aPmax" = input$annotationPhylumthreshold[2],
                    "aCmin" = input$annotationClassthreshold[1],
                    "aCmax" = input$annotationClassthreshold[2],
                    "aOmin" = input$annotationOrderthreshold[1],
                    "aOmax" = input$annotationOrderthreshold[2],
                    "aFmin" = input$annotationFamilythreshold[1],
                    "aFmax" = input$annotationFamilythreshold[2],
                    "aGmin" = input$annotationGenusthreshold[1],
                    "aGmax" = input$annotationGenusthreshold[2],
                    "aSmin" = input$annotationSpeciethreshold
                    )

    df %>% jsonlite::toJSON() %>% write_lines(values$json_name)
  }
  if(input$PairedOrNot=='y')
  {
    path_fastq_R1 = paste(tempdir(),"Masque_files_R1",sep= .Platform$file.sep)
    path_fastq_R2 = paste(tempdir(),"Masque_files_R2",sep= .Platform$file.sep)

    df = data.frame("paired"=TRUE,
                    "path_R1"=path_fastq_R1,
                    "path_R2"=path_fastq_R2,
                    "host"=input$HostName,
                    "type"=input$DataTypeMasque,
                    "mail"=values$login_email,
                    "contaminant"= path_fasta,
                    "pattern_R1"= input$R1files,
                    "phredthres" = input$phredthres,
                    "mincorrect" = input$mincorrect,
                    "minoverlap" = input$minoverlap,
                    "minreadlength" = input$minreadlength,
                    "dreptype" = input$dreptype,
                    "maxampliconlength" = input$maxampliconlength,
                    "minampliconlength" = input$minampliconlength,
                    "minabundance" = input$minabundance,
                    "clusteringthreshold" = input$clusteringthreshold,
                    "clusteringstrand" = input$clusteringstrand,
                    "annotationstrand" = input$annotationstrand,
                    "aKmin" = input$annotationKingdomthreshold,
                    "aPmin" = input$annotationPhylumthreshold[1],
                    "aPmax" = input$annotationPhylumthreshold[2],
                    "aCmin" = input$annotationClassthreshold[1],
                    "aCmax" = input$annotationClassthreshold[2],
                    "aOmin" = input$annotationOrderthreshold[1],
                    "aOmax" = input$annotationOrderthreshold[2],
                    "aFmin" = input$annotationFamilythreshold[1],
                    "aFmax" = input$annotationFamilythreshold[2],
                    "aGmin" = input$annotationGenusthreshold[1],
                    "aGmax" = input$annotationGenusthreshold[2],
                    "aSmin" = input$annotationSpeciethreshold
                    )
    df %>% jsonlite::toJSON() %>% write_lines(values$json_name)
  }
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
  taxo = NULL
  counts = NULL
  taxoCreated = FALSE
  ## Counts table
  counts = biom_data(dataBIOM)
  if(!is.null(counts))
  {
    counts = as.matrix(counts)
    ## Change of - to . is risky
    colnames(counts) = gsub("-",".",colnames(counts))
    counts = as.data.frame(counts)
  }
  CheckCounts = CheckCountsTable(counts)
  counts = CheckCounts$counts
  
  ## Taxonomy table
  tmp = observation_metadata(dataBIOM)
  if(!is.null(tmp))
  {
    if(is.data.frame(tmp)) taxo = as.data.frame(tmp)
    if(!is.data.frame(tmp)) taxo = as.data.frame(t(sapply(observation_metadata(dataBIOM),FUN=function(x){x[1:7]})))
    
    OTUnames = rownames(taxo)
    ## Modif taxo table (remove p__,... and change the colnames)
    taxo = as.data.frame(sapply(taxo,gsub,pattern="^.*__",replacement=""))
    colnames(taxo) = c("Kingdom", "Phylum","Class","Order","Family","Genus","Species")
    rownames(taxo) = OTUnames
    ## Remove empty row
    taxo[taxo==""] = NA
    taxo[taxo=="Unassigned"] = NA
    taxo=taxo[rowSums(is.na(taxo))!=dim(taxo)[2], ]
  }
  if(is.null(tmp) && !is.null(counts)) {taxo = data.frame(rownames(counts),row.names = rownames(counts));names(taxo)=NA; taxoCreated = TRUE}
  
  CheckTaxo = CheckTaxoTable(taxo,counts,taxoCreated)
  
  ## Pourcentage of annotation
  tmp = PercentAnnot(counts,taxo)
  
  return(list(counts=counts,taxo=taxo,CheckCounts=CheckCounts,CheckTaxo=CheckTaxo,Percent=tmp$Percent,CheckPercent=tmp$Error))
}

## Check the data
GetDataFromCT <-function(dataC,dataT, MGSTable)
{
  
  ## Counts table
  counts = dataC
  CheckCounts = CheckCountsTable(counts, MGSTable)
  counts = CheckCounts$counts
  
  
  ## Taxonomy table
  taxo = as.data.frame(dataT)
  CheckTaxo = CheckTaxoTable(taxo,counts, MGSTable)
  
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
  CT_Norm = NULL
  Error = NULL
  
  ## Counts and taxo tables
  CT = dataInput$counts
  taxo = dataInput$taxo
  namesTaxo = colnames(taxo)
  # save(CT,target,taxo,file="testMerge.RData")
  
  ## Select cols in the target
  labels = rownames(target)
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
    
    ## Order CT according to the target
    CT = OrderCounts(counts=CT,labels=labels)$CountsOrder
    CT_noNorm = CT
    RowProd = sum(apply(CT_noNorm,1,prod))
    merged_table = merge(CT, taxo, by="row.names")
    CT = as.data.frame(merged_table[,2: (dim(CT)[2]+1)])
    taxo = as.data.frame(merged_table[,(dim(CT)[2]+2):dim(merged_table)[2]])
    
    rownames(CT) = merged_table[,1]
    rownames(taxo) = merged_table[,1]
    colnames(taxo) = namesTaxo
    #ordOTU = order(rownames(taxo))
    counts_annot = CT
    if(0%in%colSums(counts_annot)){Error = "At least one of the column of the counts table is 0" }
    else{
      ## Create the dds object
      dds <- DESeqDataSetFromMatrix(countData=CT, colData=target, design=design,ignoreRank=TRUE)
      #save(dds,file="testdds.RData")
      if(is.null(VarNorm)){
        ## Counts normalisation
        ## Normalisation with or without 0
        if(input$AccountForNA=="NonNull" || RowProd==0) dds = estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)),geoMeans=GeoMeansCT(CT))
        if(input$AccountForNA=="All" && RowProd!=0) dds = estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)))
        if(input$AccountForNA=="Weighted" && input$AccountForNA!="NonNull" ) {dds = estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)),geoMeans=GeoMeansCT(CT)); sizeFactors(dds) = w.sizefactor(CT)}
        if(input$AccountForNA=="Total counts") { sizeFactors(dds) = colSums(CT)/mean(colSums(CT))}
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
            dds_tmp <- DESeqDataSetFromMatrix(countData=CT_tmp, colData=target_tmp, design=~labels,ignoreRank=TRUE)
            if(input$AccountForNA=="NonNull") {dds_tmp = estimateSizeFactors(dds_tmp,locfunc=eval(as.name(input$locfunc)),geoMeans=GeoMeansCT(CT_tmp)); normFactors[indgrp] = sizeFactors(dds_tmp)}
            if(input$AccountForNA=="All") {dds_tmp = estimateSizeFactors(dds_tmp,locfunc=eval(as.name(input$locfunc))); normFactors[indgrp] = sizeFactors(dds_tmp)}
            if(input$AccountForNA=="Weighted" && input$AccountForNA!="NonNull" ) {dds_tmp = estimateSizeFactors(dds_tmp,locfunc=eval(as.name(input$locfunc)),geoMeans=GeoMeansCT(CT_tmp)); normFactors[indgrp] = w.sizefactor(CT_tmp)}
            if(input$AccountForNA=="Total counts") { normFactors[indgrp] = colSums(CT_tmp)/mean(colSums(CT_tmp))}
          }
        } else{
          if(input$AccountForNA=="NonNull" || RowProd==0) dds = estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)),geoMeans=GeoMeansCT(CT))
          if(input$AccountForNA=="All" && RowProd!=0) dds = estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)))
          if(input$AccountForNA=="Weighted" && input$AccountForNA!="NonNull" ) {dds = estimateSizeFactors(dds,locfunc=eval(as.name(input$locfunc)),geoMeans=GeoMeansCT(CT)); sizeFactors(dds) = w.sizefactor(CT)}
          if(input$AccountForNA=="Total counts") { sizeFactors(dds) = colSums(CT)/mean(colSums(CT))}
          normFactors = sizeFactors(dds)
        }
        
        sizeFactors(dds) = normFactors
      }
      
      ## Keep normalized OTU table
      CT_Norm = counts(dds, normalized=TRUE)
      
      # Only interesting OTU
      # merged_table = merge(CT, taxo[order(rownames(CT)),], by="row.names")
      
      #     merged_table = merge(CT, taxo, by="row.names")
      #     CT = as.data.frame(merged_table[,2: (dim(CT)[2]+1)])
      #     taxo = as.data.frame(merged_table[,(dim(CT)[2]+2):dim(merged_table)[2]])
      #     
      #     rownames(CT) = merged_table[,1]
      #     rownames(taxo) = merged_table[,1]
      #     #ordOTU = order(rownames(taxo))
      #     counts_annot = CT
      #       ordOTU = order(rownames(taxo))
      #       indOTU_annot = which(rownames(CT)%in%rownames(taxo))
      #       counts_annot = CT[indOTU_annot[ordOTU],]
      ## Aggregate matrix
      if(taxoSelect=="OTU/Gene") counts = counts_annot
      else{
        if(input$TypeTable == "MGS" && input$FileFormat!="fileBiom"){
          MGS_taxocol = which(toupper(colnames(taxo))%in%"MGS")
          taxoS = taxo[,MGS_taxocol]
          counts = aggregate(counts_annot,by=list(Taxonomy = taxoS),mean)
          rownames(counts)=counts[,1]
          counts=counts[,-1]
          counts_int=t(apply(counts,1,as.integer))
          rownames(counts_int)=rownames(counts)
          colnames(counts_int)=colnames(counts)
          counts=counts_int
        }
        if(taxoSelect != "MGS" || input$FileFormat=="fileBiom"){
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
  }
  return(list(counts=counts,CheckTarget=CheckTarget,normFactors=normFactors, CT_noNorm=CT_noNorm, CT_Norm =CT_Norm, Error = Error))
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

## Get weighted size factors
w.sizefactor <- function(CT)
{
  sf = c()
  CT = as.matrix(CT)
  nbsamp = ncol(CT)
  CT[which(CT<1)]=NA
  gm = apply(CT,1,geometric.mean,na.rm=TRUE)
  weights = nbsamp - apply(CT,1,FUN=function(x){tmp =length(which(is.na(x))) ;return(tmp)})
  weights_tmp = weights
  
  for(i in 1:ncol(CT))
  {
    ind = which(is.na(CT[,i]))
    gm_tmp = gm
    tmp = CT[,i]
    if(length(ind)>0) {tmp = CT[-ind,i]; gm_tmp = gm[-ind]; weights_tmp = weights[-ind]}
    sf[i] = w.median(tmp/gm_tmp,weights_tmp, na.rm = TRUE)
  }
  names(sf) = colnames(CT)
  return(sf)
}

## Calculated the weighted median
w.median <- function (x, w, na.rm = TRUE) 
{
  if (missing(w)) 
    w <- rep.int(1, length(x))
  else {
    if (length(w) != length(x)) 
      stop("'x' and 'w' must have the same length")
    if (any(is.na(w))) 
      stop("NA weights not allowed")
    if (any(w < 0)) 
      stop("Negative weights not allowed")
  }
  if (is.integer(w)) 
    w <- as.numeric(w)
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  if (all(w == 0)) {
    warning("All weights are zero")
    return(NA)
  }
  o <- order(x)
  x <- x[o]
  w <- w[o]
  p <- cumsum(w)/sum(w)
  n <- sum(p < 0.5)
  if (p[n + 1] > 0.5) 
    x[n + 1]
  else (x[n + 1] + x[n + 2])/2
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


## Order the counts 
Filtered_feature <- function(counts,th.samp,th.abund)
{
  ind = NULL
  
  ## Total abundance over samples
  Tot_abundance = log(rowSums(counts)+1)
  ind.ab = which(Tot_abundance<=th.abund)
  
  ## Get the numbre of non zero sample
  counts.bin = as.matrix(counts)
  counts.bin[which(counts>0)] = 1
  nbSampByFeat = rowSums(counts.bin)
  
  ind.samp = which(nbSampByFeat<=th.samp)
  
  ind = unique(c(ind.samp,ind.ab))
  
  return(list(ind=ind,Tot_abundance=Tot_abundance,ind.ab=ind.ab,counts.bin=counts.bin,ind.samp=ind.samp,nbSampByFeat=nbSampByFeat))
}






## Order the counts 
plot_filter <- function(counts,th.samp,th.abund,type="Scatter")
{
  res = NULL
  
  ## Initial plot for plotly
  if(type == 'Abundance' || type == 'Samples'){ 
    dataNull = data.frame(x=c(0,0),y=c(1,2),col=c("white","white"))
    res = ggplot(dataNull,aes(x=x,y=y))+geom_point(aes(colour = col))+theme_bw()+ scale_color_manual(values = "white")
  }
  
  res_filter = Filtered_feature(counts,th.samp,th.abund)
  if(type == 'Abundance' && !is.null(th.samp) && !is.null(th.abund) )
  {
    state = rep("Kept",nrow(counts))
    ## Total abundance over samples
    Tot_abundance = res_filter$Tot_abundance
    
    ind = res_filter$ind.ab
    ord = order(Tot_abundance,decreasing = FALSE)
    
    ## Modify the state
    state[ind] = "Removed"
    
    ## Create the data.frame for ggplot
    df = data.frame(lab = rownames(counts)[ord],y = Tot_abundance[ord],State=state[ord])
    df$lab = factor(df$lab,levels = rownames(counts)[ord])
    df$State = factor(df$State,levels = c("Kept","Removed"))
    
    ## plot the results
    gg = ggplot(df,aes(lab,y,fill=State)) + geom_bar(stat='identity') + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))
    gg = gg + geom_hline( yintercept = th.abund,linetype = "longdash") + xlab("")
    gg = gg + scale_fill_manual(values = c('springgreen3','firebrick'))
    if(!"Kept"%in%df$State ) gg = gg + scale_fill_manual(values = 'firebrick')
    if(!"Removed"%in%df$State ) gg = gg + scale_fill_manual(values = 'springgreen3')
    
    res = gg
  }
  
  if(type == 'Samples' && !is.null(th.samp) && !is.null(th.abund))
  {
    state = rep("Kept",nrow(counts))
    
    ## Get the number of non zero sample
    nbSampByFeat = res_filter$nbSampByFeat
    ind = res_filter$ind.samp
    ord = order(nbSampByFeat,decreasing = FALSE)
    
    state[ind] = "Removed"
    
    df = data.frame(lab = rownames(counts)[ord],y = nbSampByFeat[ord],State=state[ord])
    df$lab = factor(df$lab,levels = rownames(counts)[ord])
    df$State = factor(df$State,levels = c("Kept","Removed"))
    ## plot the results
    
    gg = ggplot(df,aes(lab,y,fill=State)) + geom_bar(stat='identity') + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))
    gg = gg + geom_hline( yintercept = th.samp,linetype = "longdash") + xlab("")
    gg = gg + scale_fill_manual(values = c('springgreen3','firebrick'))  
    if(!"Kept"%in%df$State ) gg = gg + scale_fill_manual(values = 'firebrick')
    if(!"Removed"%in%df$State ) gg = gg + scale_fill_manual(values = 'springgreen3')
    
    res = gg
    
  }
  
  if(type == 'Scatter')
  {
    state = rep("Kept",nrow(counts))
    
    ## Get the number of non zero sample
    nbSampByFeat = res_filter$nbSampByFeat
    Tot_abundance = res_filter$Tot_abundance
    
    ## Get the selected features (under the thresholds)
    ind = res_filter$ind
    
    ## Modify the state
    state[ind] = "Removed"
    
    ## Create the data.frame for ggplot
    df = data.frame(lab =rownames(counts), y = nbSampByFeat,x = Tot_abundance,State=state)
    df$lab = factor(df$lab,levels = rownames(counts))
    df$State = factor(df$State,levels = c("Kept","Removed"))
    
    x_var = df$x
    y_var = df$y
    State = df$State
    PointSize = 2
    colors_scat = list(Kept="#00CD66",Removed="#b22222")

    res = scatterD3(x = x_var,
                     y = y_var,
                     lab = rownames(df),
                     xlab = "Abundance in log",
                     ylab = "Number of samples",
                     col_var = State,
                     colors = colors_scat,
                     size_lab = PointSize,
                     key_var = rownames(df),
                     point_opacity = 0.7,
                     transitions = TRUE)
    
    
#     gg = ggplot(df,aes(x,y,color=State,group = lab)) + geom_point() + theme_bw()
#     gg = gg + geom_vline( xintercept = th.samp,linetype = "longdash")
#     gg = gg + geom_hline( yintercept = th.abund,linetype = "longdash")
#     gg = gg + scale_color_manual(values = c('springgreen3','firebrick'))
#     ggplotly(gg)
#     return(gg)
  }
  
  return(res)
}



######################################################
## NAME: SelectThreshAb 
##
## INPUT:
##    infile : data matrix (counts, rows: taxo)
##    lambda : Tuning parameter (default is 500)
##
## OUTPUT:
##    Cut-off value 
##
######################################################

SelectThreshAb <- function(infile,lambda=500,graph=TRUE){
  
  rs <- rowSums(infile)
  test_Filtre <- sapply(c(min(rs):lambda),FUN=function(x) table(rs>x))
  if(!is.list(test_Filtre))
  {
    x <- c(min(rs):lambda)
    reslm <- lm(test_Filtre[1,]~x)$coefficients
    val <- which(test_Filtre[1,]>reslm[1])[1]
    if (graph==TRUE){
      plot(test_Filtre[1,],ylab="Nb especes avec moins de x comptages sur tous les echantillons",xlab="x")
      abline(v=val,col="red")
    } 
  } else val = min(rs)
  return(max(val, min(rs)+1))
}

