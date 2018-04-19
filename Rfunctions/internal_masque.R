Project_status <- function(masque_key,curdir){

  passOK = FALSE;status = NULL;file = NULL
  json_files = list.files(paste(curdir,"www","masque",sep= .Platform$file.sep),pattern = paste("file", masque_key, ".json", sep=""),recursive = TRUE)
  allpass = gsub(gsub(json_files,pattern = ".*file",replacement = ""),pattern = ".json",replacement = "")

  if(length(allpass)>0){
    passOK = any(masque_key==allpass)
    if(passOK){
      ind = which(masque_key==allpass)[1]
      file = paste(curdir,"www","masque",json_files[ind],sep= .Platform$file.sep)
      status = gsub(json_files[ind],pattern = "/.*",replacement = "")
    }
  }
  
  # return(list(status=status,file=json_file,passOK=passOK))
  return(list(status=status,file=file,passOK=passOK))
}




Project_box_result <- function(masque_key,curdir){

  res = NULL
  PS = Project_status(masque_key,curdir)
  
  if(PS$passOK){
    if(PS$status=="done")
    {
      res = list()
      json_file = PS$file
      folder_name = paste('file',masque_key,sep="")
      json_file = paste(curdir,"www","masque","done",paste(folder_name, ".json", sep=""),sep= .Platform$file.sep)
      json_data = rjson::fromJSON(file=json_file)
      
      ### Paste file name as folder
      annot_process = paste(curdir,"www","masque","done",folder_name,"shaman_process_annotation.tsv",sep= .Platform$file.sep)

      ## Waiting for file creation (max 3min)
      start = Sys.time(); diff = 0
      while(!file.exists(annot_process) && diff<180){
        annot_process = paste(curdir,"www","masque","done",folder_name,"shaman_process_annotation.tsv",sep= .Platform$file.sep)
        tmp = Sys.time()
        diff = tmp-start
      }

      if(file.exists(annot_process))
      {
        ap = read.csv(annot_process,sep="\t")
        if(json_data$type == "16S") db_choices = c("Silva" = "silva","Greengenes" = "greengenes", "RDP"= "rdp")
        else if(json_data$type == "23S_28S" || json_data$type == "18S") db_choices = c("Silva" = "silva","RDP"= "rdp")
        else db_choices = c("Findley" = "findley", "Underhill"= "underhill", "Unite"= "unite", "RDP"= "rdp")
        #build_process = paste(curdir,"www","masque","done",folder_name,"shaman_process_annotation.tsv",sep= .Platform$file.sep)
        res[[1]] = fluidRow(
          #HTML('<center><h1><strong>Your project is done !</strong></h1> <br/> <em><h4> Hereafter is a summary of the building and annotation processes</h4> </em><br/> </center>'),
          HTML('<center><h1><strong>Your project is done !</strong></h1></center>'),
          br(),
          column(width=4,
                 h3("OTU building process"),
                 shinydashboard::valueBox(ap$Count[1],tags$strong(tags$h5("Number of amplicons", style = "width: 70%;")), color = "light-blue", width = NULL, icon = uiOutput("amplicon_icon")),
                 shinydashboard::valueBox(ap$Count[2],tags$strong(tags$h5("Remaining amplicons after dereplication", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("dereplication_icon")),
                 shinydashboard::valueBox(ap$Count[3],tags$strong(tags$h5("Remaining amplicons after removing singletons", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("singleton_icon")),
                 shinydashboard::valueBox(ap$Count[4],tags$strong(tags$h5("Remaining amplicons after removing chimera", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("chimera_icon")),
                 shinydashboard::valueBox(ap$Count[5],tags$strong(tags$h5("Number of OTU", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("otu_icon"))
          ),
          #column(width=5, div(img(src = "masque.png",width = "50%",height = "20%"))),
          if(json_data$type == "16S"){
            column(width=5,
                   strong(h3("16S annotation process")),
                   shinydashboard::valueBox(ap$Count[6],tags$strong(tags$h5("Number of OTU annotated by SILVA", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("silva_icon")),
                   shinydashboard::valueBox(ap$Count[7],tags$strong(tags$h5("Number of OTU annotated by Greengenes", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("greengenes_icon")),
                   shinydashboard::valueBox(ap$Count[8],tags$strong(tags$h5("Number of OTU annotated by RDP", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("rdp_icon"))
            )
          }
          else if(json_data$type == "18S"){
            column(width=5,
                   strong(h3("18S annotation process")),
                   shinydashboard::valueBox(ap$Count[6],tags$strong(tags$h5("Number of OTU annotated by SILVA", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("silva_icon")),
                   shinydashboard::valueBox(ap$Count[7],tags$strong(tags$h5("Number of OTU annotated by RDP", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("rdp_icon"))
            )
          }
          else if(json_data$type == "23S_28S"){
            column(width=5,
                   strong(h3("23S/28S annotation process")),
                   shinydashboard::valueBox(ap$Count[6],tags$strong(tags$h5("Number of OTU annotated by SILVA", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("silva_icon")),
                   shinydashboard::valueBox(ap$Count[7],tags$strong(tags$h5("Number of OTU annotated by RDP", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("rdp_icon"))
            )
          }
          else{
            column(width=5,
                   strong(h3("ITS annotation process")),
                   shinydashboard::valueBox(ap$Count[6],tags$strong(tags$h5("Number of OTU annotated by Unite", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("unite_icon")),
                   shinydashboard::valueBox(ap$Count[7],tags$strong(tags$h5("Number of OTU annotated by Findley", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("findley_icon")),
                   shinydashboard::valueBox(ap$Count[8],tags$strong(tags$h5("Number of OTU annotated by Underhill", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("underhill_icon")),
                   shinydashboard::valueBox(ap$Count[9],tags$strong(tags$h5("Number of OTU annotated by RDP", style = "width: 70%;")), color = "light-blue", width = NULL,icon = uiOutput("rdp_icon"))
            )
         },
         #offset = 4,
         column(width=3,
                strong(h3("Start statistical analysis")),
                box(title="Load the results",width = NULL, background = "light-blue",
                    selectInput("masque_db","Select the database",choices=db_choices),
                    conditionalPanel(condition="input.masque_db=='rdp'",numericInput("rdp_thres",h6(strong("Threshold:")),0.5,step=0.01,min=0.01,max=1)),
                    actionButton("LoadResMasque", "Upload the results",icon=icon('upload')),
                    tags$style(type='text/css', "#LoadResMasque { width:100%; }"),
                    receiveSweetAlert(messageId = "LoadResMasque")
                ),
                #shinydashboard::box(title="Upload the results",width = NULL, background = "light-blue",
                #selectInput("masque_database","Select the database",choices=c("Silva" = "silva","Greengenes" = "greengenes", "MARDE"="merde")),
                #tags$style(type='text/css', "#masque_database { width:100%; margin-top: 5px;}"),
                #actionButton("RunResMasque",label = "Upload the results",icon=icon('upload')),
                #tags$style(type='text/css', "#RunResMasque { width:100%; margin-top: 15px;}")
                #),
                #receiveSweetAlert(messageId = "WTF2"),
                box(id="box-zip",title="Download .zip file",width = NULL, status = "success",
                    downloadButton('Download_masque_zip', 'Download the results'),
                    tags$style(type='text/css', "#Download_masque_zip { width:100%;}")
                )
         )
        )
      } else{res =HTML('<center><h1><strong>Your project is done !</strong></h1> <br/> <em><h4>Result can not be shown...</h4> </em> </center>')}
      
      build_process = paste(curdir,"www","masque","done",folder_name,"shaman_process_build.tsv",sep= .Platform$file.sep)
      if(file.exists(build_process))
      {
        ap = read.csv(build_process,sep="\t")
      
        res[[2]] = fluidRow(column(width=12,strong(h3("Detailed process table")),DT::dataTableOutput("build_process_table")))
      }
      
      
    }
    
    if(PS$status=="doing"){
      res = fluidRow(
        HTML('<center><h1><strong>Your project is currently running !</strong></h1> <br/> <br/> </center>'),
        inlineCSS(gaugeCSS),
        gaugeOutput("gaugeMasque_progress", width = "100%", height = "100%")
        #jscode <- "shinyjs.refresh = function() { history.go(0); }",
        #useShinyjs(),
        #extendShinyjs(text = jscode),
      )
    }
    
    if(PS$status=="error"){
      error_message = "Failed"
      json_file = PS$file
      error_file = paste(curdir,"www","masque","error",paste('file',masque_key,"_error.txt",sep=""),sep= .Platform$file.sep)
  
      if(file.exists(error_file)){error_message = read_lines(error_file)}
      
      res = fluidRow(
                HTML('<center><h1><strong>The workflow failed during progression</strong></h1> <br/> </center>'),
        
                column(width = 12,
                       div(style = "background-color: white; margin: 0 auto;width: 50%; text-align:center;border:1px solid red",
                           h4(strong("Workflow error message")),
                           hr(style = "width: 70%;text-align:left;"),
                           HTML(paste(error_message,collapse = " <br/> ")),
                           br()
                       )
                )
          )
    }
  }
  
  return(list(box=res,PS=PS))
  
}