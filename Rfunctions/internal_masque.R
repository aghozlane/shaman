Project_status <- function(masque_key,curdir){

  passOK = FALSE;status = NULL;file = NULL
  
  json_files = list.files(paste(curdir,"www","masque",sep= .Platform$file.sep),pattern = "json",recursive = TRUE)
  allpass = gsub(gsub(json_files,pattern = ".*file",replacement = ""),pattern = ".json",replacement = "")
  
  
  if(length(allpass)>0){
    passOK = any(masque_key==allpass)
    if(passOK){
      ind = which(masque_key==allpass)[1]
      file = paste(curdir,"www","masque",json_files[ind],sep= .Platform$file.sep)
      status = gsub(json_files[ind],pattern = "/.*",replacement = "")
    }
  }
  
  return(list(status=status,file=file,passOK=passOK))
}




Project_box_result <- function(masque_key,curdir){
  res = NULL
  PS = Project_status(masque_key,curdir)
  
  
  if(PS$status=="done")
  {
    res = list()
    json_file = PS$file
    folder_name = paste('file',masque_key,sep="")
    
    ### Paste file name as folder
    annot_process = paste(curdir,"www","masque","done",folder_name,"shaman_annotation_process.tsv",sep= .Platform$file.sep)
    
    
    ## Waiting for file creation (max 3min)
    start = Sys.time(); diff = 0
    while(!file.exists(annot_process) && diff<180){
      annot_process = paste(curdir,"www","masque","done",folder_name,"shaman_annotation_process.tsv",sep= .Platform$file.sep)
      tmp = Sys.time()
      diff = tmp-start
    }
    
    if(file.exists(annot_process))
    {
      ap = read.csv(annot_process,sep="\t")
      build_process = paste(curdir,"www","masque","done",folder_name,"shaman_annotation_process.tsv",sep= .Platform$file.sep)
      
      res[[1]] = fluidRow(
        HTML('<center><h1><strong>Your project is done !</strong></h1> <br/> <em><h4> Hereafter is a summary of the building and annotation processes</h4> </em><br/> </center>'),
        br(),
        column(width=5,
               valueBox(ap$Count[1],tags$strong(tags$h5("Number of amplicons", style = "width: 70%;")), color = "light-blue",width=NULL,icon = uiOutput("amplicon_icon")),
               valueBox(ap$Count[2],tags$strong(tags$h5("Remaining amplicons after dereplication", style = "width: 70%;")), color = "light-blue",width=NULL,icon = uiOutput("dereplication_icon")),
               valueBox(ap$Count[3],tags$strong(tags$h5("Remaining amplicons after removing singletons", style = "width: 70%;")), color = "light-blue",width=NULL,icon = uiOutput("singleton_icon")),
               valueBox(ap$Count[4],tags$strong(tags$h5("Remaining amplicons after removing chimera", style = "width: 70%;")), color = "light-blue",width=NULL,icon = uiOutput("chimera_icon"))
        ),
        column(width=5,offset = 2,
        box(title="Upload the results",width = NULL, background = "light-blue",
            selectInput("masque_database","Select the database",choices=c("Silva" = "silva","Greengenes" = "greengenes")),
            tags$style(type='text/css', "#masque-database { width:100%; margin-top: 5px;}"),
            actionButton("RunResMasque",label = "Upload the results",icon=icon('upload')),
            tags$style(type='text/css', "#RunResMasque { width:100%; margin-top: 15px;}")
        ),
        
        box(id="box-zip",title="Download .zip file",width = NULL, status = "success",
            downloadButton('Download_masque_zip', 'Download the results'),
            tags$style(type='text/css', "#Download_masque_zip { width:100%; margin-top: 15px;}")
        )
        
        
        )
        
        )
      
      
    } else{res =HTML('<center><h1><strong>Your project is done !</strong></h1> <br/> <em><h4> Hereafter is a summary of the building and annotation processes</h4> </em> </center>')}
    
    
    build_process = paste(curdir,"www","masque","done",folder_name,"shaman_build_process.tsv",sep= .Platform$file.sep)
    if(file.exists(build_process))
    {
      ap = read.csv(build_process,sep="\t")
    
      res[[2]] = fluidRow(column(width=12,strong(h3("Building process table")),DT::dataTableOutput("build_process_table")))
    }
    
    
  }
  
  if(PS$status=="doing"){
    res = fluidRow(
      HTML('<center><h1><strong>Your project is currently running !</strong></h1> <br/> <br/> </center>'),
      inlineCSS(gaugeCSS),
      gaugeOutput("gaugeMasque_progress", width = "100%", height = "100%")
    )
  }
  
  if(PS$status=="error"){
    
    json_file = PS$file
    error_file = paste(curdir,"www","masque","error",paste('file',masque_key,"_error",".txt",sep=""),sep= .Platform$file.sep)

    if(file.exists(error_file)){error_message = read_lines(error_file)}
    
    res = fluidRow(
              HTML('<center><h1><strong>Sorry, the workflow failed during progression</strong></h1> <br/> <em><h4> Hereafter is the message error.</h4> </em> <br/> </center>'),
      
              column(width = 12,
                     div(style = "background-color: white; margin: 0 auto;width: 50%; text-align:center;border:1px solid red",
                         h4(strong("Error message")),
                         hr(style = "width: 70%;"),
                         HTML(paste(error_message,collapse = " <br/> ")),
                         br()
                     )
              )
        )
  }
  
  return(list(box=res,PS=PS))
  
}