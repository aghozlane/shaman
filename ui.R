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

if (!require(genefilter)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("genefilter")
  library(genefilter)
}

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Home", tabName = "Home", icon = icon("home")),
    menuItem("Tutorial", tabName = "Tutorial", icon = icon("book")),
    menuItem("Upload your data", tabName = "Upload", icon = icon("upload")),
    menuItemOutput("dymMenu"),
    img(src = "logo.jpg", height = 49, width = 220,style="position:absolute;bottom:0;margin:0 0 15px 10px;")
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "Home",
            div(style="width:100% ; max-width: 1200px",
              tabBox(title="Welcome to SHAMAN", id="tabset1", width=NULL,
                   tabPanel("About", tags$script(type="text/javascript", language="javascript", src="google-analytics.js"),
                            p("SHAMAN is a SHiny application for Metagenomic ANalysis including the normalization,
                                       the differential analysis and mutiple visualization.",style = "font-family: 'times'; font-si16pt"),
                            p("SHAMAN is based on DESeq2 R package", a("[Anders and Huber 2010]", href="http://www.ncbi.nlm.nih.gov/pubmed/20979621"), "for the analysis of metagenomic data, as suggested in", a("[McMurdie and Holmes 2014]",href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3974642/"),
                              ". SHAMAN robustly identifies the differential abundant genera with the Generalized Linear Model implemented in DESeq2", a("[Love 2014]", href="http://www.ncbi.nlm.nih.gov/pubmed/25516281"),".
                              Resulting p-values are adjusted according to the Benjamini and Hochberg procedure [Benjamini and Hochberg 1995].
                              The PCOA is performed with the", a("ade4 R package",href="http://pbil.univ-lyon1.fr/ade4/"), "and plots are generated with", a("ggplot2",href="http://ggplot2.org/"), "or", a("D3.js packages",href="http://d3js.org/"), ".
                              A presentation about SHAMAN is available", a("here",target="_blank",href="shaman_presentation.pdf"),style = "font-family: 'times'; font-si16pt"),
                            p("Hereafter is the global workflow of the SHAMAN application:"),
                            div(img(src = "Workflow.png",width = "100%",style="max-width: 600px"),Align="center")
                            ),
                   tabPanel("Authors", h3("The main contributors to SHAMAN:"),
                            p(a("Stevenn Volant", href="mailto:stevenn.volant@pasteur.fr"), "(Initiator, coding, testing, documentation, evaluation)"),
                            p(a("Amine Ghozlane",href="mailto:amine.ghozlane@pasteur.fr"), "(Coding, testing, documentation, evaluation)"),
                            p(a("Hugo Varet",href="mailto:hugo.varet@pasteur.fr"), "(Coding, testing, feature suggestions)"),
                            p(a("Pierre Lechat",href="mailto:pierre.lechat@pasteur.fr"), "(Coding, testing, feature suggestions)"),
                            p(a("Marie-Agnès Dillies",href="mailto:marie-agnes.dillies@pasteur.fr"), "(Evaluation)"),
                            p(a("Sean Kennedy",href="mailto:sean.kennedy@pasteur.fr"), "(Evaluation)"),
                            h3("Acknowledgements"),
                            p("Thanks to the following people for patches and other suggestions for improvements:"),
                            p(a("Christophe Malabat",href="mailto:christophe.malabat@pasteur.fr"),a(", Julien Tap",href="mailto:julien.tap@danone.com"),a(", Anna Zhukova",href="mailto:anna.zhukova@pasteur.fr"), a(", Rachel Torchet",href="mailto:rachel.torchet@pasteur.fr"))
                          ),
                   tabPanel("Citing SHAMAN",p("No papers about SHAMAN have been published yet, but a manuscript is in preparation.",style = "font-family: 'times'; font-si16pt"),
                            p("For now, if you have any comments, questions or suggestions, or want to use SHAMAN please contact authors", a("here", href="mailto:shaman@pasteur.fr"),".", style = "font-family: 'times'; font-si16pt; color:red")
            )))
    ),

    tabItem(tabName = "Tutorial",
            div(style="width:100% ; max-width: 1200px",
            tabBox(title="How to use SHAMAN", id="tabset1", width =NULL,
            tabPanel("Introduction",
            p(" You can test SHAMAN with the dataset from", a("[Tap et al. 2015]",href="http://www.ncbi.nlm.nih.gov/pubmed/26235304"),
              ", which is available", a("here",target="_blank",href="Alimintest.zip"),"."),
            p("The zip archive contains 4 different files :", br(),
              "- the otu count matrix : Alimintest_otu_table.csv,", br(),
              "- the otu annotation table : Alimintest_otu_table.csv,", br(),
              "- the experimental design : Alimintest_target.csv,", br(),
              "- the contrast table : Alimintest_contrasts.csv."),
            p("Two groups of person follow two  strict diet periods that involve the intake of 40g following 10g of fiber per day, or 10g of fiber after a 40g fiber intake period:"),
            img(src = "tutorial/FigS1.png",width = "100%",style="max-width: 900px"),
            p("The 16S rRNA (V3 - V4 regions) from fece samples was sequenced at time stamp : 2, 3, 4 and 5.", br(),
              "The analysis will consider the different impact of the different fiber intake and the comparison to patient metabolic data.")),
            tabPanel("1-Load 16S data",
            p("The first step consists to load the count matrix and the annotation table as follow :"),
            p("- Select 'Upload your data'", br(),
              "- Load the count table :",br(), img(src = "tutorial/tutorial_upload1.png",width = "100%",style="max-width: 900px"),hr(),
              "- Load the annotation table :", br(), img(src = "tutorial/tutorial_upload2.png",width = "100%",style="max-width: 900px"),hr(),
              "- When successfully loaded, the tables are accessible as bellow :",br(),
              img(src = "tutorial/tutorial_upload3.png",width = "100%",style="max-width: 600px"),img(src = "tutorial/tutorial_upload4.png",width = "100%",style="max-width: 600px"))),
            tabPanel("2-Differential analysis",
            p("The second step consists to load the experimental design and the contrast table as follow :"),
            p("- Select 'Run differential analysis'",br(),
              "- Load the target file :",br(),img(src = "tutorial/tutorial_target.png",width = "100%",style="max-width: 800px"),hr(),
              "- Identify the taxonomy level where the analysis will be performed :",br(),img(src = "tutorial/tutorial_target1.png",width = "100%",style="max-width: 800px"),hr(),
              "- Identify the interactions :",br(),img(src = "tutorial/tutorial_target2.png",width = "100%",style="max-width: 800px"),hr(),
              "- Run the analysis :",br(),img(src = "tutorial/tutorial_target3.png",width = "100%",style="max-width: 800px"),hr(),
              "- When successfully loaded, the tables are accessible as bellow :",br(),img(src = "tutorial/tutorial_target4.png",width = "100%",style="max-width: 800px")),hr(),
            p("- Finally, load the contrast file :",br(),img(src = "tutorial/tutorial_contraste.png",width = "100%",style="max-width: 800px"),hr(),
              "- Contrasts can be visualized as follow :",br(),img(src = "tutorial/tutorial_contraste1.png",width = "100%",style="max-width: 400px"))),
            tabPanel("3-Diagnostic plots",
            p("'Diagnostic plots' section provides several visualizations to control the analysis",br(),
              "- Total mapped read count",br(),img(src="tutorial/tutorial_total_barplot.png",width = "100%",style="max-width: 900px"),hr(),
              "- Nul barplot count",br(),img(src="tutorial/tutorial_nul_barplot.png",width = "100%",style="max-width: 900px"),hr(),
              "- Maj taxonomy count",br(),img(src="tutorial/tutorial_maj_taxonomy.png",width = "100%",style="max-width: 900px"),hr(),
              "- Density of counts",br(),img(src="tutorial/tutorial_density.png",width = "100%",style="max-width: 900px"),hr(),
              "- Size factors vs total number of reads",br(),img(src="tutorial/tutorial_size_factor.png",width = "100%",style="max-width: 900px"),hr(),
              "- PCA",br(),img(src="tutorial/tutorial_pca.png",width = "100%",style="max-width: 900px"),hr(),
              "- PCOA",br(),img(src="tutorial/tutorial_pcoa.png",width = "100%",style="max-width: 900px"),hr(),
              "- Clustering",br(),img(src="tutorial/tutorial_clustering.png",width = "100%",style="max-width: 900px"))),
            tabPanel("4-Differential analysis results",
            img(src = "tutorial/tutorial_table.png",width = "100%",style="max-width: 900px")),
            tabPanel("5-Visualization",
                     p("'Diagnostic plots' section provides several visualization to control the analysis",br(),
                     "- Barplot",br(),img(src="tutorial/tutorial_barplot.png",width = "100%",style="max-width: 900px"),hr(),
                     "- Heatmap",br(),img(src="tutorial/tutorial_heatmap.png",width = "100%",style="max-width: 900px"),hr(),
                     "- Boxplot",br(),img(src="tutorial/tutorial_boxplot.png",width = "100%",style="max-width: 900px"),hr(),
                     "- Diversity",br(),img(src="tutorial/tutorial_diversity.png",width = "100%",style="max-width: 900px"),hr(),
                     "- Rarefaction",br(),img(src="tutorial/tutorial_rarefaction.png",width = "100%",style="max-width: 900px")))))
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
                  selectInput("FileFormat","",c("Count table & taxonomy"="fileCounts","BIOM file"="fileBiom"),selected="fileCounts")
                ),
                conditionalPanel(condition="input.FileFormat=='fileCounts'",
                  box(title="Load the count table",width = 3,height = "250px", status = "primary", solidHeader = TRUE,collapsible = FALSE,
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
              column(width=3,infoBoxOutput("InfoContrast",width=NULL))
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
                    column(width=3,checkboxInput("AccountForNA","Compute geometric mean without 0",value=TRUE))
                    # column(width=3,uiOutput("RefSelect"))
                  )
                ),
                fluidRow(
                column(width=8,
#                        uiOutput("contrastBox"),
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
                                              "Maj. taxonomy"="MajTax", "Density"="densityPlot", "Dispersion" = "DispPlot",
                                              "Size factors VS total"="SfactorsVStot", "PCA"="pcaPlot", "PCoA"="pcoaPlot","Clustering" = "clustPlot"))
                    ),
                box(title = "Options",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= FALSE,
                    conditionalPanel(condition="input.DiagPlot!='clustPlot' && input.DiagPlot!='pcaPlot' && input.DiagPlot!='SfactorsVStot' && input.DiagPlot!='DispPlot'",
                                     radioButtons("CountsType","Counts:",c("Normalized"="norm","Raw"="raw"),inline = TRUE)
                                     ),
                    conditionalPanel(condition="input.DiagPlot!='Sfactors' && input.DiagPlot!='SfactorsVStot' ",uiOutput("VarIntDiag")),
                    conditionalPanel(condition="input.DiagPlot=='pcoaPlot' || input.DiagPlot=='pcaPlot'",
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
                    conditionalPanel(condition="input.DiagPlot=='SfactorsVStot' || input.DiagPlot=='pcaPlot' || input.DiagPlot=='pcoaPlot'",column(width=6,sliderInput("cexLabelDiag", h6("Points"),min=0,max=5,value = 1,step =0.1)))
                    
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
                                    uiOutput("VarIntVisu"),
                                   h5(strong("Select the modalities")),
                                   uiOutput("ModVisu")
                  ),
                  
                  conditionalPanel(condition="input.PlotVisuSelect!='Rarefaction' && input.PlotVisuSelect!='Diversity'",
                                   radioButtons("SelectSpecifTaxo","Select the features",c("Most abundant"="Most","All"="All", "Differential features" = "Diff", "Non differential features" = "NoDiff"))
                  ),
                  conditionalPanel(condition="input.PlotVisuSelect!='Rarefaction' && input.PlotVisuSelect!='Diversity' && (input.SelectSpecifTaxo=='Diff' || input.SelectSpecifTaxo=='NoDiff') ",
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
                                 uiOutput("ColBoxplot"),
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
                                 selectInput("colors", label=h6(strong("Gradient of colors")),choices = c("green-blue", "blue-white-red", "purple-white-orange", "red-yellow-green"),selected = "blue-white-red")
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
            column(width=3,
                   p(strong("Krona plot")),
                   img(src="Krona.png",height = 200, width = 220),
                   a(href = "test_krona.html",target="_blank", "Click Here!")
#                    tableOutput("krona") 

),
            column(width=3,
            p(strong("Tree abundance")),
            img(src="Tree.png",height = 200, width = 220),
            a(href = "http://genopole.pasteur.fr/SynTView/flash/TreeAbundance",target="_blank", "Click Here!")
            )
        )
 #includeHTML("file:///home/aghozlan/workspace/SHAMAN_App/www/text.krona.html")
  )
 )
)
  ## GOOGLE ANALYTIC
 #tags$head(includeScript("google-analytics.js"))
  ## Logo SHAMAN
  dbHeader <- dashboardHeader(title = "SHAMAN")
  dbHeader$children[[2]]$children <-  tags$a(tags$img(src='akuaku.png',height='40',width='50',style="margin:5px 0 5px 0;",align='left'), tags$h3("SHAMAN",style="font-family:Purisa; margin:15px 25px 5px 0;color:white;"))
  
# Put them together into a dashboardPage
  dashboardPage(skin="blue",
    dbHeader,
    sidebar,
    body
  )


