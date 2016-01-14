if(!require(shinydashboard)){
  install.packages('shinydashboard')
  library(shinydashboard)
}

if (!require(psych)) {
  install.packages('psych')
  
  library(psych)
  
}

if (!require(ggplot2)) {
  
  install.packages('ggplot2')
  
  library(ggplot2)
  
}

if (!require(vegan)) {
  
  install.packages('vegan')
  
  library(vegan)
  
}

if (!require(dendextend)) {
  
  install.packages('dendextend')
  
  library(dendextend)
  
}

if (!require(circlize)) {
  
  install.packages('circlize')
  
  library(circlize)
  
}

if (!require(biom)) {
  
  install.packages('biom')
  
  library(biom)
  
}

if (!require(DT)) {
  
  install.packages('DT')
  
  library(DT)
  
}

if (!require(RColorBrewer)) {
  
  install.packages('RColorBrewer')
  
  library(RColorBrewer)
  
}

if (!require(gplots)) {
  
  install.packages('gplots')
  
  library(gplots)
  
}

if (!require(DESeq2)) {
  
  source("https://bioconductor.org/biocLite.R")
  
  biocLite("DESeq2")
  
  library(DESeq2)
  
}

if (!require(ade4)) {
  
  install.packages('ade4')
  
  library(ade4)
  
}

sidebar <- dashboardSidebar(
  
  sidebarMenu(
    
    menuItem("Home", tabName = "Home", icon = icon("home")),
    
    menuItem("Tutorial", tabName = "Tutorial", icon = icon("book")),
    
    menuItem("Upload your data", tabName = "Upload", icon = icon("upload")),
    
    menuItemOutput("dymMenu")
    
  )
  
)

body <- dashboardBody(
  
  tabItems(
    
    tabItem(tabName = "Home",
            
            tabBox(title="Welcome to Meta16S", id="tabset1", height = "700", width = 12,
                   
                   tabPanel("About", p("Meta16S is a web interface for the analysis of metagenomic data including the normalization,
                                       
                                       the differential analysis and mutiple visualization.",style = "font-family: 'times'; font-si16pt"),
                            
                            p("Meta16S is based on DESeq2 R package", a("[Anders and Huber 2010]", href="http://www.ncbi.nlm.nih.gov/pubmed/20979621"), "for the analysis of metagenomic data, as suggested in", a("[McMurdie and Holmes 2014]",href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3974642/"),
                              
                              ". Meta16S robustly identifies the differential abundant genera with the Generalized Linear Model implemented in DESeq2", a("[Love 2014]", href="http://www.ncbi.nlm.nih.gov/pubmed/25516281"),".
                              
                              Resulting p-values are adjusted according to the Benjamini and Hochberg procedure [Benjamini and Hochberg 1995].
                              
                              The PCOA is performed with the", a("ade4 R package",href="http://pbil.univ-lyon1.fr/ade4/"), "and plots are generated with", a("ggplot2",href="http://ggplot2.org/"), "or", a("D3.js packages",href="http://d3js.org/"), ".", style = "font-family: 'times'; font-si16pt")),
                   
                   tabPanel("Authors", h3("The main contributors to Meta16s:"),
                            
                            p(a("Stevenn Volant", href="mailto:stevenn.volant@pasteur.fr"), "(Initiator, coding, testing, documentation, evaluation)"),
                            
                            p(a("Amine Ghozlane",href="mailto:amine.ghozlane@pasteur.fr"), "(Coding, testing, documentation, evaluation)"),
                            
                            p(a("Hugo Varet",href="mailto:hugo.varet@pasteur.fr"), "(Coding, testing, feature suggestions)"),
                            
                            p(a("Pierre Lechat",href="mailto:pierre.lechat@pasteur.fr"), "(Coding, testing, feature suggestions)"),
                            
                            p(a("Marie-Agnès Dillies",href="mailto:marie-agnes.dillies@pasteur.fr"), "(Evaluation)"),
                            
                            p(a("Sean Kennedy",href="mailto:sean.kennedy@pasteur.fr"), "(Evaluation)"),
                            
                            h3("Acknowledgements"),
                            
                            p("Thanks to the following people for patches and other suggestions for improvements:"),
                            
                            p(a("Christophe Malabat",href="mailto:christophe.malabat@pasteur.fr"))),
                   
                   tabPanel("Citing Meta16S",p("No papers about Meta16s have been published yet, but a manuscript is in preparation.",style = "font-family: 'times'; font-si16pt"),
                            
                            p("For now, if you have any comments, questions or suggestions, or want to use meta16s please contact", a("Dr. Marie-Agnès Dillies", href="mailto:marie-agnes.dillies@pasteur.fr"),".", style = "font-family: 'times'; font-si16pt; color:red"))
            ),
                   img(src = "logo.jpg", height = 90, width = 600,align="right")
    ),      

    
    tabItem(tabName = "Tutorial",
            
            h2("How to !")
            
    ),
    
    tabItem(tabName = "Upload",
            tags$style(type='text/css', ".well { max-width: 20em; }"),
            fluidRow(
              column(width=3,valueBoxOutput("valueErrorPercent",width=NULL)),
              column(width=3,infoBoxOutput("InfoErrorCounts",width=NULL)),
              column(width=3,infoBoxOutput("InfoErrorTaxo",width=NULL))
            ),
            br(),
             fluidRow(
                box(title="Select your file format",width = 3,status = "success", solidHeader = TRUE,collapsible = FALSE,
                  selectInput("FileFormat","",c("Counts table & taxonomy"="fileCounts","BIOM file"="fileBiom"),selected="fileCounts")
                ),
                conditionalPanel(condition="input.FileFormat=='fileCounts'",
                  box(title="Load the counts table",width = 3,height = "250px", status = "primary", solidHeader = TRUE,collapsible = FALSE,
                    fileInput('fileCounts', h6(strong('Select your file')),width="100%")
                  ),
                  
                  box(title="Load the taxonomy file",width = 3,height = "250px", status = "primary", solidHeader = TRUE,collapsible = FALSE,
                      fluidRow(
                        column(width=6,radioButtons("TypeTaxo",h6(strong("Format:")),c("Table"="Table","RDP"="RDP"))),
                        column(width=6,
                             conditionalPanel(condition="input.TypeTaxo=='RDP'",numericInput("RDP_th",h6(strong("Threshold:")),0.5,step=0.01,min=0.01,max=1))
                        )
                      ),
                      fileInput('fileTaxo', h6(strong('Select your file')),width="100%")
                  )
                  
                ),
                
                conditionalPanel(condition="input.FileFormat=='fileBiom'",
                                 box(title="Load the BIOM file",width = 3,height = "250px", status = "primary", solidHeader = TRUE,collapsible = FALSE,
                                     fileInput('fileBiom', h5(strong('Select your file')),width="100%")
                                 )           
                )
             ),
              column(12,uiOutput("TabBoxData"))

              
          
    ),
    
  #### Statistical analysis

    tabItem(tabName = "RunDiff",
            fluidRow(
              column(width=3,infoBoxOutput("RowTarget",width=NULL)),
              column(width=3,infoBoxOutput("InfoTaxo",width=NULL)),
              column(width=3,infoBoxOutput("InfoDESeq",width=NULL)),
              column(width=3,conditionalPanel(condition="input.RunDESeq>=1",infoBoxOutput("InfoContrast",width=NULL)))
            ),            
            fluidRow(
              column(width=5,
                box(title="Experimental design",width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = FALSE,
                  fluidRow(
                    column(width=6,fileInput('fileTarget', h6(strong('Select your target file')),width="100%")),
                    column(width=6,uiOutput("SelectTaxo"))
                  ),
                  fluidRow( 
                    column(width=6,uiOutput("SelectInterestVar")),
                    column(width=6,uiOutput("SelectInteraction2")),
                    column(width=6,uiOutput("RunButton"))
                  )
                ),
                uiOutput("BoxTarget"),
                uiOutput("BoxCountsMerge")
              ),
       
              column(width=7,
                box(title="Options",width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = TRUE,
                  fluidRow(
                    column(width=3,
                      radioButtons("TransType",h6(strong("Type of transformation")),choices = c("VST"="VST","rlog"="rlog"))
                    ), 
                    column(width=3,
                      radioButtons("IndFiltering",h6(strong("Independent filtering")),choices = c("True"=TRUE,"False"=FALSE))
                    ),
                    column(width=3,
                      radioButtons("AdjMeth",h6(strong("p-value adjustement")),choices = c("BH"="BH","BY"="BY"))
                    ),
                    column(width=3,
                      textInput("AlphaVal",h6(strong("Level of significance")),value=0.05)
                    )
                  ),
                  fluidRow(
                    column(width=3,
                           radioButtons("CooksCutOff",h6(strong("Cooks cut-off")),choices = c("Auto"='Auto',"No cut-off"=Inf,"Value"="val")),
                           conditionalPanel(condition="input.CooksCutOff=='val'",textInput("CutOffVal",h6("Cut-off:"),value=0))
                    ),
                    
                    column(width=3,
                      radioButtons("locfunc",h6(strong("Local function")),choices = c("Median"="median","Shorth"="shorth"))
                    ),  
                    column(width=3,
                      radioButtons("fitType",h6(strong("Relationship")),choices = c("Parametric"="parametric","Local"="local"))
                    ),
                    column(width=3,
                      conditionalPanel(condition="input.FileFormat=='fileCounts' && input.TypeTaxo=='RDP'",
                        sliderInput("ThreshProba",h6(strong("Probability threshold (rdp annotation)")),min=0.01, max=1,value=0.5,step = 0.01)
                      )
                    ),
                    column(width=3,uiOutput("RefSelect"))
                  )
                ),
                fluidRow(
                column(width=8,
                       uiOutput("contrastBox"),
                       uiOutput("contrastBoxAdvanced")
                       ),
                column(width=4,
                       uiOutput("contrastDefined")
                )  
                ) 
              )
            )
            
    ),
    tabItem(tabName = "DiagPlotTab",
            fluidRow(
              column(width=9,
                     
                plotOutput("PlotDiag",height="100%"),
                br(),
                conditionalPanel(condition="input.DiagPlot=='SfactorsVStot'",
                  box(title = "Size factors",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                    dataTableOutput("SizeFactTable")
                  )
                ),
                  
                conditionalPanel(condition="input.DiagPlot=='pcaPlot'",
                                 box(title = "Eigen values",  width = 6, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= FALSE,
                                    plotOutput("PlotEigen",height="100%")
                                 )
                ),
                conditionalPanel(condition="input.DiagPlot=='pcoaPlot'",
                                 box(title = "Eigen values",  width = 6, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= FALSE,
                                     plotOutput("PlotpcoaEigen",height="100%")
                                 )
                )
                
              ),
              column(width=3,
                box(title = "Select your plot",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = FALSE,collapsed= FALSE,
                  selectInput("DiagPlot","",c("Total barplot"="barplotTot","Nul barplot"="barplotNul",
                                              "Maj. taxonomy"="MajTax", "Density"="densityPlot", 
                                              "Size factors VS total"="SfactorsVStot", "PCA"="pcaPlot", "PCoA"="pcoaPlot","Clustering" = "clustPlot",
                                              "SERE" = "SERE"))
                    ),
                box(title = "Options",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= FALSE,
                    
                    conditionalPanel(condition="input.DiagPlot!='Sfactors' && input.DiagPlot!='SfactorsVStot' ",uiOutput("VarIntDiag")),
                    conditionalPanel(condition="input.DiagPlot=='pcoaPlot'",
                                     h5(strong("Select the modalities")),
                                     uiOutput("ModMat")
                    ),
                    conditionalPanel(condition="input.DiagPlot=='pcoaPlot' || input.DiagPlot=='SERE' || input.DiagPlot=='clustPlot' ",
                      selectInput("DistClust","Distance",c("euclidean", "SERE"="sere", "canberra", "bray", "kulczynski", "jaccard", 
                                                         "gower", "altGower", "morisita", "horn","mountford","raup","binomial",
                                                         "chao","cao","mahalanobis"),selected="jaccard")
                    )
                    
#                 conditionalPanel(condition="input.RadioPlotBi=='Nuage'",selectInput("ColorBiplot", "Couleur",choices=c("Bleue" = 'blue',"Rouge"='red',"Vert"='green', "Noir"='black'),width="50%")),
#                 sliderInput("TransAlphaBi", "Transparence",min=1, max=100, value=50, step=1),
#                 conditionalPanel(condition="input.RadioPlotBi!='Nuage'", radioButtons("SensGraphBi","Sens du graph",choices=c("Vertical"="Vert","Horizontal"="Hori"))),
#                 conditionalPanel(condition="input.RadioPlotBi=='box'", checkboxInput("CheckAddPointsBoxBi","Ajouter les données",value=FALSE)) 
               ),
                box(
                  title = "Appearance",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
#                     conditionalPanel(condition="input.DiagPlot=='Sfactors'",
#                                      h6(strong("Layout")),
#                                      numericInput("NbcolSfactors", h6("Columns"),min=1,value = NA)
#                     ),
                  sliderInput("heightDiag", "Height",min=100,max=1500,value = 500,step =10),

                  conditionalPanel(condition="input.DiagPlot=='clustPlot'",
                                   h6(strong("Layout")),
                                   selectInput("typeHculst", h6("Type"),c("Horizontal"="hori","Fan"="fan")),
                                   checkboxInput("colorHC","Add color",value=TRUE)
                  ),
                  conditionalPanel(condition="input.DiagPlot=='pcoaPlot'",  selectInput("labelPCOA","Label type",c("Group", "Sample"),selected="Group"),
                                   #checkboxInput("colorgroup","Same color for the group",value=FALSE),
                                   sliderInput("cexcircle", "Circle size",min=0,max=2,value = 0.9,step =0.1),
                                   sliderInput("cexpoint", "Point size",min=0,max=3,value = 1,step =0.1),
                                   sliderInput("cexstar", "Star height",min=0,max=1,value = 0.95,step =0.1)
                                   
                  ),
                  conditionalPanel(condition="input.DiagPlot=='SfactorsVStot'",
                    checkboxInput("addLabelSFact","Add label",FALSE)
                  ),

                  fluidRow(
                    column(width=12, p(strong("Size"))),
                    column(width=6,sliderInput("cexTitleDiag", h6("Axis"),min=0,max=5,value = 1,step =0.1)),
                    conditionalPanel(condition="input.DiagPlot=='SfactorsVStot' || input.DiagPlot=='pcaPlot'",column(width=6,sliderInput("cexLabelDiag", h6("Points"),min=0,max=5,value = 1,step =0.1)))
                    
                  )

#                   sliderInput("widthDiag", "width",min=100,max=1500,value = 1000,step =10)
                 

              ),
              box(title = "Export",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                  selectInput("Exp_format",h5(strong("Export format")),c("png"="png","pdf"="pdf","eps"="eps","svg"="svg"), multiple = FALSE),
                  fluidRow(
                  column(width=6,numericInput("heightDiagExport", "Height (in px)",min=100,max=NA,value = 500,step =1)),
                  column(width=6,numericInput("widthDiagExport", "Width (in px)",min=100,max=NA,value = 500,step =1))
                  ),
                  downloadButton("exportdiag", "Export")
                  
#                   downloadButton("exportPDFdiag", "Download pdf"),
#                   downloadButton("exportPNGdiag", "Download png"),
#                   downloadButton("exportEPSdiag", "Download eps"),
#                   downloadButton("exportSVGdiag", "Download svg"),
                  
              )
            )
        )
    ),
    tabItem(tabName = "TableDiff",
            fluidRow(
              column(width=9,uiOutput("TabBoxDataDiff")),
              column(width=3,
                     box(title = "Select your contrast",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= FALSE,
                         selectInput("ContrastList_table",h6(strong("Contrast list")),"", multiple = FALSE),
                         htmlOutput("ContrastOverviewTable")
                     ),
                     box(title = "Export",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                         fluidRow(column(width=12,selectInput("WhichExportTable", "Select the table to export",c("Complete"="Complete","Up"="Up","Down"="Down")))),
                         uiOutput("ExportTableButton")   
                     )
            )
            )
    ),
  
  #### Data Viz
  
  tabItem(tabName = "Visu",
          fluidRow(
            column(width=9,
                  uiOutput("plotVisu")
            ),

            column(width=3,
              box(title = "Select your plot",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = FALSE,collapsed= FALSE,
                  selectizeInput("PlotVisuSelect","",c("Barplot"="Barplot","Heatmap"="Heatmap","Boxplot"="Boxplot","Diversity"="Diversity","Rarefaction"="Rarefaction"),selected = "Barplot")
              ),
              
              
              
              ########################################################################
              ###
              ###               Options Visualization
              ###
              ########################################################################
              
              
              box(title = "Options",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= FALSE,
                  
                  conditionalPanel(condition="input.PlotVisuSelect!='Rarefaction'",
                                    uiOutput("VarIntVisu")
                  ),
                  conditionalPanel(condition="input.PlotVisuSelect!='Rarefaction' && input.PlotVisuSelect!='Diversity'",
                                   radioButtons("SelectSpecifTaxo","Select the features",c("Most abundant"="Most","All"="All", "Differential features" = "Diff"))
                  ),
                  conditionalPanel(condition="input.PlotVisuSelect!='Rarefaction' && input.PlotVisuSelect!='Diversity' && input.SelectSpecifTaxo=='Diff'",
                                   selectizeInput("ContrastList_table_Visu","",choices = "", multiple = FALSE)
                  ),
                  conditionalPanel(condition="input.PlotVisuSelect!='Rarefaction' && input.PlotVisuSelect!='Diversity'",
                                   uiOutput("TaxoToPlotVisu")
                  ),

                  
                  
                ##################
                ## BARPLOT
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Barplot'",
                                 selectizeInput(inputId = "CountsOrProp",label = h6(strong("Type of data")),choices = c("Proportions" = "prop", "Counts" = "counts"),selected = "prop")
                ),
                
                
                ##################
                ## HEATMAP
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Heatmap'",
                                 selectizeInput(inputId = "HeatMapType",label = h6(strong("Type of data")),choices = c("Counts" = "Counts", "Log2FC" = "Log2FC"),selected = "Counts")                                 
                ),
                conditionalPanel(condition="input.PlotVisuSelect=='Heatmap' && input.HeatMapType=='Log2FC'",
                                 selectizeInput("ContrastList_table_FC",h6(strong("Contrasts (Min = 2)")),choices = "", multiple = TRUE)
                ),
                conditionalPanel(condition="input.PlotVisuSelect=='Heatmap'",
                                 selectizeInput(inputId = "scaleHeatmap",label = h6(strong("Scale:")),choices = c("None" = "none", "Rows" = "row", "Column" = "col"),selected = "none")
                                 
                ),
                
                ##################
                ## BOXPLOT
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Boxplot'",
                                 selectizeInput("typeDataBox",h6(strong("Type of data")),c("Log2"="Log2","Relative"="Relative"))
                ),
                
                ##################
                ## DIVERSITY
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Diversity'",
                                 selectizeInput("WhichDiv",h6(strong("Diversity")),c('Alpha','Beta','Gamma'),selected  = c('Alpha','Beta','Gamma'),multiple=TRUE),
                                 checkboxInput("AddBoxplotDiv","AddBoxplot",value=FALSE)
                ),
                conditionalPanel(condition="input.PlotVisuSelect=='Diversity' && input.AddBoxplotDiv",
                                 uiOutput("SelectVarBoxDiv")
                )
              ),
              
              
              ########################################################################
              ###
              ###               Appearance Visualization
              ###
              ########################################################################
              
              
              box(title = "Appearance",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                sliderInput("heightVisu", h6(strong("Height")),min=100,max=4000,value = 800),


                
                ##################
                ## BOXPLOT
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Boxplot'",
                                 radioButtons("ScaleBoxplot","Scales",c("Fixed"="fixed","Free"="free"),inline=TRUE),
                                 checkboxInput("CheckAddPointsBox","Add points",value=TRUE)
                ),
                
                ##################
                ## DIVERSITY
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Diversity'",
                                 sliderInput("sizePointGlobal", h6(strong("Points size")),min=0.5,max=10,value =3,step=0.5),
                                 checkboxInput("SplitVisuGlobal","Split diversity",value=FALSE)
                ),
                
                ##################
                ## HEATMAP
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Heatmap' && input.HeatMapType!='Log2FC'",
                                 selectInput("colors", label=h6(strong("Gradient of colors")),choices = c("green-blue", "blue-white-red", "purple-white-orange", "red-yellow-green"))
                ),
                conditionalPanel(condition="input.PlotVisuSelect=='Heatmap'",
                                 fluidRow(
                                   column(width=12,h6(strong("Labels options"))),
                                   column(width=6,sliderInput("LabelSizeHeatmap", h6("Size"),min=0.1,max=2,value = 0.7,step = 0.1)),
                                   column(width=6,sliderInput("LabelOrientHeatmap", h6("Orientation"),min=0,max=90,value = 0,step = 5)),
                                   column(width=6,sliderInput("LabelColOffsetHeatmap", h6("Column offset"),min=0,max=4,value = 0,step = 0.5)),
                                   column(width=6,sliderInput("LabelRowOffsetHeatmap", h6("Row offset"),min=0,max=4,value = 0,step = 0.5)),
                                   
                                   column(width=12,h6(strong("Margins options"))),
                                   column(width=6,sliderInput("rightMargin", h6("Right"),min=0,max=20,value = 6,step = 1)),
                                   column(width=6,sliderInput("lowerMargin", h6("Lower"),min=0,max=20,value = 6,step = 1))
                                 )
                ),              
                
                ##################
                ## ALL
                ##################

                conditionalPanel(condition="input.PlotVisuSelect!='Rarefaction'",
                                 radioButtons(inputId = "SensPlotVisu",label = h6(strong("Orientation")),choices = c("Vertical" = "Vertical", "Horizontal" = "Horizontal"),selected = "Vertical",inline = TRUE)
                )
                
              ),
              box(title = "Export",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                  ##################
                  ## BARPLOT
                  ##################
                  conditionalPanel(condition="input.PlotVisuSelect=='Barplot'",
                                   radioButtons("positionBarPlot","Position",c("Grouped"="fill","Stacked"="dodge"),inline=TRUE)
                  ),
                  
                  selectInput("Exp_format_Visu",h5(strong("Export format")),c("png"="png","pdf"="pdf","eps"="eps","svg"="svg"), multiple = FALSE),
                  fluidRow(
                    column(width=6,numericInput("heightVisuExport", "Height (in px)",min=100,max=NA,value = 500,step =1)),
                    column(width=6,numericInput("widthVisuExport", "Width (in px)",min=100,max=NA,value = 500,step =1))
                  ),
                  downloadButton("exportVisu", "Export")
                  
              )
              
            )
          )
  ),
  



  #### Krona plot
  tabItem(tabName = "Krona",
 
         fluidRow(
            column(width=3,tableOutput("krona") )
 #includeHTML("file:///home/aghozlan/workspace/META16S_App/www/text.krona.html")
          )  
  )       
              
#              column(width=3,infoBoxOutput("NumberColBox",width=NULL)),
#               column(width=3,infoBoxOutput("NumberRowBox",width=NULL)),
#               column(width=3,infoBoxOutput("NumberQuantiBox",width=NULL)),
#               column(width=3,infoBoxOutput("NumberQualiBox",width=NULL))
#             ),
#             fluidRow(
# 
#               column(9,
#                      box(title = "Table de données", width = NULL, status = "success", solidHeader = TRUE,
#                       dataTableOutput("DataBrutes")
#                      )
#               ),
#               column(width=3,
#                      box(title = "Libéllés", width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
#                          h5(strong('Votre fichier contient t-il ...')),
#                          radioButtons('header', h5(strong('... des entêtes de colonnes ?')),choices=list("Oui" = 1,"Non" = 0),selected=1),
#                          radioButtons('rnames', h5(strong('... des libellés de ligne ?')),choices=list("Oui" = 1,"Non" = 0),selected=1)
#                      ),
#                      
#                      box(title = "Choix des variables", width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
#               column(width=3,infoBoxOutput("NumberColBox",width=NULL)),
#               column(width=3,infoBoxOutput("NumberRowBox",width=NULL)),
#               column(width=3,infoBoxOutput("NumberQuantiBox",width=NULL)),
#               column(width=3,infoBoxOutput("NumberQualiBox",width=NULL))
#             ),
#             fluidRow(
# 
#               column(9,
#                      box(title = "Table de données", width = NULL, status = "success", solidHeader = TRUE,
#                       dataTableOutput("DataBrutes")
#                      )
#               ),
#               column(width=3,
#                      box(title = "Libéllés", width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
#                          h5(strong('Votre fichier contient t-il ...')),
#                          radioButtons('header', h5(strong('... des entêtes de colonnes ?')),choices=list("Oui" = 1,"Non" = 0),selected=1),
#                          radioButtons('rnames', h5(strong('... des libellés de ligne ?')),choices=list("Oui" = 1,"Non" = 0),selected=1)
#                      ),
#                      
#                      box(title = "Choix des variables", width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
#                      box(title = "Choix des variables", width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
#               column(width=3,infoBoxOutput("NumberColBox",width=NULL)),
#               column(width=3,infoBoxOutput("NumberRowBox",width=NULL)),
#               column(width=3,infoBoxOutput("NumberQuantiBox",width=NULL)),
#               column(width=3,infoBoxOutput("NumberQualiBox",width=NULL))
#             ),
#             fluidRow(
# 
#               column(9,
#                      box(title = "Table de données", width = NULL, status = "success", solidHeader = TRUE,
#                       dataTableOutput("DataBrutes")
#                      )
#               ),
#               column(width=3,
#                      box(title = "Libéllés", width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
#                          h5(strong('Votre fichier contient t-il ...')),
#                          radioButtons('header', h5(strong('... des entêtes de colonnes ?')),choices=list("Oui" = 1,"Non" = 0),selected=1),
#                          radioButtons('rnames', h5(strong('... des libellés de ligne ?')),choices=list("Oui" = 1,"Non" = 0),selected=1)
#                      ),
#                      
#                      box(title = "Choix des variables", width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
#                          uiOutput("SelectDataVar"),
#                          actionButton("RefreshData",icon=icon("refresh"),strong("Actualiser"))
#                      ),
#                      
#                      box(title = "Transformer des variables", width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
#                          uiOutput("SelectVarQuanti"),
#                          actionButton("GoQuali",icon=icon("arrow-circle-right"),strong("Transformer")),
#                          uiOutput("QuantiToQuali"))  
#               )
#             )
#             
#     ),
#     tabItem(tabName = "AnaStat",
#             
#           p("De nombreuses analyses statistiques peuvent être réalisées à partir de vos données. Decouvrez l'information
#             contenue dans vos données (analyse unvariée) et étudier les relations entre vos différentes variables.")
#     ),
#     
#     tabItem(tabName = "DescGene",
#             fluidRow(
#               
#               column(width=9,
#                      tabBox(side = "left",id = 'tabs', width=NULL, status = "primary",
#                                  tabPanel("Variables quantitatives",value=1,
#                                           h4(strong("Statistiques descriptives (variables quantitatives)")),
#                                           dataTableOutput("TableQuanti")
#                                         
#                                  ),
#                                  tabPanel("Variables qualitatives",value=2,
#                                           h4(strong("Statistiques descriptives (variables qualitatives)")),
#                                           dataTableOutput("TableQuali")
#                                  ),
#                                  tabPanel("Table des corrélations",value=3,
#                                           h4(strong("Table des corrélations")),
#                                           dataTableOutput("CorTable")
#                                  )
# 
#                                  
#                      ),
#                      
#                      conditionalPanel(condition="input.TestCor==true",
#                                       box(title="Résultat du test",width = NULL, status = "success", solidHeader = TRUE,collapsible = TRUE,
#                                                dataTableOutput("CorTableTest")
#                                                
#                                                
#                                       ))
#               ),
#               column(width=3,
#                      box(title = "Choix de l'analyse", width = NULL, status = "primary", solidHeader = TRUE,
#                          conditionalPanel(condition = "input.tabs==1",
#                                           selectInput("IndicQuanti",p(strong("Choisissez les indicateurs"),h6(em("Sélection multiple avec CTRL"))),
#                                                       c("Nb valeurs"=1,
#                                                         "Nb manquants"=2,
#                                                         "Somme"=3,
#                                                         "Min"=4,
#                                                         "1er Quartile"=5,
#                                                         "Mediane"=6,
#                                                         'Moyenne'=7,
#                                                         "3eme Quartile"=8,
#                                                         "Max"=9,
#                                                         "Variance"=10,
#                                                         "Ecart-type"=11,
#                                                         "Coeff Variation"=12),
#                                                       selected=c(1,2,3,4,5,6,7,10),multiple=TRUE,size=2,selectize=FALSE),
#                                           actionButton("RefreshStat",icon=icon("refresh"),strong("Refresh"))
#                          ),
#                          conditionalPanel(condition = "input.tabs==2",
#                                           selectInput("IndicQuali",p(strong("Choisissez les indicateurs"),h6(em("Sélection multiple avec CTRL"))),
#                                                       c("Effectifs"=1,
#                                                         "%"=2,
#                                                         "% cumulés"=3),
#                                                       selected=c(1,2,3),multiple=TRUE,size=2,selectize=FALSE),
#                                           actionButton("RefreshStatQuali",icon=icon("refresh"),strong("Refresh"))
#                          ),
#                          
#                          conditionalPanel(condition = "input.tabs==3",
#                                           radioButtons("CorelMeth","Choix de la corrélation",
#                                                        choices=c("Pearson"="pearson","Spearman"="spearman","Kendall"="kendall")),
#                                           checkboxInput("TestCor","Test de corrélation",value=FALSE)
#                          )
#                      ),
#                      box(title = "Aide", width = NULL, status = "warning", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE
#                      )
#               )
#             )
#     ),
#     
#     tabItem(tabName = "AnalUni",
#             fluidRow(
#               column(width = 9,
#               box(
#                 title = "Représentation",status = "success", solidHeader = TRUE,collapsible = TRUE,width=NULL,
#                 plotOutput("plotuni"),
#                 p(Align='right',
#                   downloadButton("exportPDFuni", "Download pdf"),
#                   downloadButton("exportPNGuni", "Download png")
#                 )
#               ),
#               conditionalPanel(condition="input.TestMoyCible==true", infoBoxOutput("ResTTestBox",width=6))    
#             ),
#             column(width = 3,
#                        box(
#                          title = "Choix de l'analyse", width = NULL, status = "primary", solidHeader = TRUE,
#                          
#                          uiOutput("SelectUniVar"),
#                          uiOutput("RadioUniPlot"),
#                          checkboxInput("TestMoyCible","T-test",value=FALSE),
#                          conditionalPanel(condition="input.TestMoyCible==true", 
#                                           textInput("ValCibleTtest", label ="Valeur cible", value = 0,width="50%"),
#                                           textInput("alphaTtest", label ="Seuil (en %)", value = 5,width="50%"),
#                                           actionButton("ExecuteTtestCible",icon=icon("sign-in"),strong("Exécuter"))
#                          )
#                        ),
#                        box(
#                          title = "Options du graphique",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
#                          
#                          conditionalPanel(condition="input.RadioPlotUni=='BarPlot'",sliderInput("widthBarPlot", "Taille des barres",min=1, max=100, value=50, step=1)),
#                          conditionalPanel(condition="input.RadioPlotUni=='Pie'", sliderInput("PieWidth", "Largeur", min=1, max=100, value=100, step=1)),
#                          conditionalPanel(condition="input.RadioPlotUni!='Pie'", sliderInput("SizeQQplot", "Taille",min=1, max=10, value=2, step=0.5)),
#                          
#                          selectInput("ColorUniplot", "Couleur",choices=c("Bleue" = 'blue',"Rouge"='red',"Vert"='green', "Noir"='black'),width="50%"),
#                          sliderInput("TransAlphaUni", "Transparence",min=1, max=100, value=50, step=1),
#                          conditionalPanel(condition="input.RadioPlotUni=='hist'",sliderInput("binwidth", "Taille de la fenètre glissante",min=1, max=100, value=1, step=1)),
#                          conditionalPanel(condition="input.RadioPlotUni=='hist'", radioButtons("HistDens","Ordonnées",choices=c("Comptages"="counts","Fréquences"="freq"))),
#                          conditionalPanel(condition="input.RadioPlotUni=='hist'", checkboxInput("CheckDens","Ajouter la densité",value=FALSE)),
#                          conditionalPanel(condition="input.RadioPlotUni!='qqplot'", radioButtons("SensGraph","Sens du graph",choices=c("Vertical"="Vert","Horizontal"="Hori"))),
#                          conditionalPanel(condition="input.RadioPlotUni=='box'", checkboxInput("CheckAddPointsBox","Ajouter les données",value=FALSE)),
#                          conditionalPanel(condition="input.RadioPlotUni=='BarPlot'",checkboxInput("BarCircular","Représentation circulaire",FALSE))
#                          
#                        )
#               )
#             )
#             
#             
#             
#             ),
#     
#     tabItem(tabName = "AnalBi",
#     
#             fluidRow(
#               column(width=9, 
#                      box(
#                        title = "Représentation",  width = NULL, status = "success", solidHeader = TRUE,collapsible = TRUE,
#                        
#                      plotOutput("biplot",
#                                 dblclick = dblclickOpts(
#                                   id = "biplot_dblclick"
#                                 ),
#                                 hover = hoverOpts(
#                                   id = "biplot_hover"
#                                 ),
#                                 brush = brushOpts(
#                                   id = "biplot_brush",
#                                   resetOnNew = TRUE
#                                 )    
#                      ),
#                      p(Align='right',
#                      downloadButton("exportPDFbi", "Download pdf"),
#                      downloadButton("exportPNGbi", "Download png")
#                      )
#                      ),
#                      conditionalPanel(condition="input.TestMoy==true", infoBoxOutput("ResTTest2sampBox",width=6)), 
#                      conditionalPanel(condition="input.RadioPlotBi=='Nuage'",infoBoxOutput("biplot_info",width=6))
#               ),
#               column(width=3,
#                      box(
#                      title = "Choix de l'analyse",  width = NULL, status = "primary", solidHeader = TRUE,
#                        
#                      uiOutput("SelectDataVarBi1"),
#                      uiOutput("SelectDataVarBi2"),
#                      uiOutput("RadioBiPlot"),
#                      uiOutput("CheckTestBi"),
#                      conditionalPanel(condition="input.TestBi!=''",textInput("alphaTtestMoy", label ="Seuil (en %)", value = 5,width="50%")),
#                      conditionalPanel(condition="input.RadioPlotBi=='Nuage'", h5(strong("Modèle linéaire")),checkboxInput("CheckLM", "Droite de régression",value = FALSE))
#                      
#                      ),
#                      box(
#                        title = "Options du graphique",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
#                        
#                        sliderInput("SizePoint", "Taille",min=1, max=10, value=2, step=0.5),
#                        conditionalPanel(condition="input.RadioPlotBi=='hist'",sliderInput("binwidthBi", "Taille de la fenètre glissante",min=1, max=100, value=1, step=1)),
#                        conditionalPanel(condition="input.RadioPlotBi=='Nuage'",selectInput("ColorBiplot", "Couleur",choices=c("Bleue" = 'blue',"Rouge"='red',"Vert"='green', "Noir"='black'),width="50%")),
#                        sliderInput("TransAlphaBi", "Transparence",min=1, max=100, value=50, step=1),
#                        conditionalPanel(condition="input.RadioPlotBi!='Nuage'", radioButtons("SensGraphBi","Sens du graph",choices=c("Vertical"="Vert","Horizontal"="Hori"))),
#                        conditionalPanel(condition="input.RadioPlotBi=='box'", checkboxInput("CheckAddPointsBoxBi","Ajouter les données",value=FALSE))
#                        
#                      )
#                    
#               )
#               
#               
#               
#             )
             )
     
    
    
  )


# Put them together into a dashboardPage
dashboardPage(skin="blue",
  dashboardHeader(title = "Meta16S"),
  sidebar,
  body
)


