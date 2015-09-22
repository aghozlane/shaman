
library(shinydashboard)
library(DT)
library(biom)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Home", tabName = "Home", icon = icon("home")),
    menuItem("Upload your data", tabName = "Upload", icon = icon("upload")),
    menuItemOutput("dymMenu")
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "Home",
            h2("Bienvenue !")
    ),
    tabItem(tabName = "Upload",
            tags$style(type='text/css', ".well { max-width: 20em; }"),
            # Tags:
            tags$head(
              tags$style(type="text/css", "select[multiple] { width: 100%; height:10em}"),
              tags$style(type="text/css", "select { width: 100%}"),
              tags$style(type="text/css", "input { width: 19em; max-width:100%}")
            ),
            fluidRow(
              column(12,
                  h3(strong("Instructions")),
                  p("Décrire le format des différents fichiers")
              ),
              column(3,
                box(title="Select your file format",width = NULL, status = "success", solidHeader = TRUE,collapsible = FALSE,
                  selectInput("FileFormat","",c("Counts table & taxonomy"="fileCounts","BIOM file"="fileBiom"),selected="fileCounts"),
                  uiOutput("LoadButton")
                )
              ),
              column(9,
                conditionalPanel(condition="input.FileFormat=='fileCounts'",
                  box(title="Load the counts table",width = 4, status = "primary", solidHeader = TRUE,collapsible = FALSE,
                    fileInput('fileCounts', h5(strong('Select your file')),width="100%")
                  ),
                  box(title="Load the taxonomy file",width = 4, status = "primary", solidHeader = TRUE,collapsible = FALSE,
                      fileInput('fileTaxo', h5(strong('Select your file')),width="100%")
                      
                  )
                  
                ),
                
                conditionalPanel(condition="input.FileFormat=='fileBiom'",
                                 box(title="Load the BIOM file",width = 4, status = "primary", solidHeader = TRUE,collapsible = FALSE,
                                     fileInput('fileBiom', h5(strong('Select your file')),width="100%")
                                 )           
                )
                
              )
              
            )
    ),
    tabItem(tabName = "datatable",
            
            fluidRow(
              column(width=3,infoBoxOutput("NumberColBox",width=NULL)),
              column(width=3,infoBoxOutput("NumberRowBox",width=NULL)),
              column(width=3,infoBoxOutput("NumberQuantiBox",width=NULL)),
              column(width=3,infoBoxOutput("NumberQualiBox",width=NULL))
            ),
            fluidRow(

              column(9,
                     box(title = "Table de données", width = NULL, status = "success", solidHeader = TRUE,
                      dataTableOutput("DataBrutes")
                     )
              ),
              column(width=3,
                     box(title = "Libéllés", width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                         h5(strong('Votre fichier contient t-il ...')),
                         radioButtons('header', h5(strong('... des entêtes de colonnes ?')),choices=list("Oui" = 1,"Non" = 0),selected=1),
                         radioButtons('rnames', h5(strong('... des libellés de ligne ?')),choices=list("Oui" = 1,"Non" = 0),selected=1)
                     ),
                     
                     box(title = "Choix des variables", width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                         uiOutput("SelectDataVar"),
                         actionButton("RefreshData",icon=icon("refresh"),strong("Actualiser"))
                     ),
                     
                     box(title = "Transformer des variables", width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                         uiOutput("SelectVarQuanti"),
                         actionButton("GoQuali",icon=icon("arrow-circle-right"),strong("Transformer")),
                         uiOutput("QuantiToQuali"))  
              )
            )
            
    ),
    tabItem(tabName = "AnaStat",
            
          p("De nombreuses analyses statistiques peuvent être réalisées à partir de vos données. Decouvrez l'information
            contenue dans vos données (analyse unvariée) et étudier les relations entre vos différentes variables.")
    ),
    
    tabItem(tabName = "DescGene",
            fluidRow(
              
              column(width=9,
                     tabBox(side = "left",id = 'tabs', width=NULL, status = "primary",
                                 tabPanel("Variables quantitatives",value=1,
                                          h4(strong("Statistiques descriptives (variables quantitatives)")),
                                          dataTableOutput("TableQuanti")
                                        
                                 ),
                                 tabPanel("Variables qualitatives",value=2,
                                          h4(strong("Statistiques descriptives (variables qualitatives)")),
                                          dataTableOutput("TableQuali")
                                 ),
                                 tabPanel("Table des corrélations",value=3,
                                          h4(strong("Table des corrélations")),
                                          dataTableOutput("CorTable")
                                 )

                                 
                     ),
                     
                     conditionalPanel(condition="input.TestCor==true",
                                      box(title="Résultat du test",width = NULL, status = "success", solidHeader = TRUE,collapsible = TRUE,
                                               dataTableOutput("CorTableTest")
                                               
                                               
                                      ))
              ),
              column(width=3,
                     box(title = "Choix de l'analyse", width = NULL, status = "primary", solidHeader = TRUE,
                         conditionalPanel(condition = "input.tabs==1",
                                          selectInput("IndicQuanti",p(strong("Choisissez les indicateurs"),h6(em("Sélection multiple avec CTRL"))),
                                                      c("Nb valeurs"=1,
                                                        "Nb manquants"=2,
                                                        "Somme"=3,
                                                        "Min"=4,
                                                        "1er Quartile"=5,
                                                        "Mediane"=6,
                                                        'Moyenne'=7,
                                                        "3eme Quartile"=8,
                                                        "Max"=9,
                                                        "Variance"=10,
                                                        "Ecart-type"=11,
                                                        "Coeff Variation"=12),
                                                      selected=c(1,2,3,4,5,6,7,10),multiple=TRUE,size=2,selectize=FALSE),
                                          actionButton("RefreshStat",icon=icon("refresh"),strong("Refresh"))
                         ),
                         conditionalPanel(condition = "input.tabs==2",
                                          selectInput("IndicQuali",p(strong("Choisissez les indicateurs"),h6(em("Sélection multiple avec CTRL"))),
                                                      c("Effectifs"=1,
                                                        "%"=2,
                                                        "% cumulés"=3),
                                                      selected=c(1,2,3),multiple=TRUE,size=2,selectize=FALSE),
                                          actionButton("RefreshStatQuali",icon=icon("refresh"),strong("Refresh"))
                         ),
                         
                         conditionalPanel(condition = "input.tabs==3",
                                          radioButtons("CorelMeth","Choix de la corrélation",
                                                       choices=c("Pearson"="pearson","Spearman"="spearman","Kendall"="kendall")),
                                          checkboxInput("TestCor","Test de corrélation",value=FALSE)
                         )
                     ),
                     box(title = "Aide", width = NULL, status = "warning", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE
                     )
              )
            )
    ),
    
    tabItem(tabName = "AnalUni",
            fluidRow(
              column(width = 9,
              box(
                title = "Représentation",status = "success", solidHeader = TRUE,collapsible = TRUE,width=NULL,
                plotOutput("plotuni"),
                p(Align='right',
                  downloadButton("exportPDFuni", "Download pdf"),
                  downloadButton("exportPNGuni", "Download png")
                )
              ),
              conditionalPanel(condition="input.TestMoyCible==true", infoBoxOutput("ResTTestBox",width=6))    
            ),
            column(width = 3,
                       box(
                         title = "Choix de l'analyse", width = NULL, status = "primary", solidHeader = TRUE,
                         
                         uiOutput("SelectUniVar"),
                         uiOutput("RadioUniPlot"),
                         checkboxInput("TestMoyCible","T-test",value=FALSE),
                         conditionalPanel(condition="input.TestMoyCible==true", 
                                          textInput("ValCibleTtest", label ="Valeur cible", value = 0,width="50%"),
                                          textInput("alphaTtest", label ="Seuil (en %)", value = 5,width="50%"),
                                          actionButton("ExecuteTtestCible",icon=icon("sign-in"),strong("Exécuter"))
                         )
                       ),
                       box(
                         title = "Options du graphique",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                         
                         conditionalPanel(condition="input.RadioPlotUni=='BarPlot'",sliderInput("widthBarPlot", "Taille des barres",min=1, max=100, value=50, step=1)),
                         conditionalPanel(condition="input.RadioPlotUni=='Pie'", sliderInput("PieWidth", "Largeur", min=1, max=100, value=100, step=1)),
                         conditionalPanel(condition="input.RadioPlotUni!='Pie'", sliderInput("SizeQQplot", "Taille",min=1, max=10, value=2, step=0.5)),
                         
                         selectInput("ColorUniplot", "Couleur",choices=c("Bleue" = 'blue',"Rouge"='red',"Vert"='green', "Noir"='black'),width="50%"),
                         sliderInput("TransAlphaUni", "Transparence",min=1, max=100, value=50, step=1),
                         conditionalPanel(condition="input.RadioPlotUni=='hist'",sliderInput("binwidth", "Taille de la fenètre glissante",min=1, max=100, value=1, step=1)),
                         conditionalPanel(condition="input.RadioPlotUni=='hist'", radioButtons("HistDens","Ordonnées",choices=c("Comptages"="counts","Fréquences"="freq"))),
                         conditionalPanel(condition="input.RadioPlotUni=='hist'", checkboxInput("CheckDens","Ajouter la densité",value=FALSE)),
                         conditionalPanel(condition="input.RadioPlotUni!='qqplot'", radioButtons("SensGraph","Sens du graph",choices=c("Vertical"="Vert","Horizontal"="Hori"))),
                         conditionalPanel(condition="input.RadioPlotUni=='box'", checkboxInput("CheckAddPointsBox","Ajouter les données",value=FALSE)),
                         conditionalPanel(condition="input.RadioPlotUni=='BarPlot'",checkboxInput("BarCircular","Représentation circulaire",FALSE))
                         
                       )
              )
            )
            
            
            
            ),
    
    tabItem(tabName = "AnalBi",
    
            fluidRow(
              column(width=9, 
                     box(
                       title = "Représentation",  width = NULL, status = "success", solidHeader = TRUE,collapsible = TRUE,
                       
                     plotOutput("biplot",
                                dblclick = dblclickOpts(
                                  id = "biplot_dblclick"
                                ),
                                hover = hoverOpts(
                                  id = "biplot_hover"
                                ),
                                brush = brushOpts(
                                  id = "biplot_brush",
                                  resetOnNew = TRUE
                                )    
                     ),
                     p(Align='right',
                     downloadButton("exportPDFbi", "Download pdf"),
                     downloadButton("exportPNGbi", "Download png")
                     )
                     ),
                     conditionalPanel(condition="input.TestMoy==true", infoBoxOutput("ResTTest2sampBox",width=6)), 
                     conditionalPanel(condition="input.RadioPlotBi=='Nuage'",infoBoxOutput("biplot_info",width=6))
              ),
              column(width=3,
                     box(
                     title = "Choix de l'analyse",  width = NULL, status = "primary", solidHeader = TRUE,
                       
                     uiOutput("SelectDataVarBi1"),
                     uiOutput("SelectDataVarBi2"),
                     uiOutput("RadioBiPlot"),
                     uiOutput("CheckTestBi"),
                     conditionalPanel(condition="input.TestBi!=''",textInput("alphaTtestMoy", label ="Seuil (en %)", value = 5,width="50%")),
                     conditionalPanel(condition="input.RadioPlotBi=='Nuage'", h5(strong("Modèle linéaire")),checkboxInput("CheckLM", "Droite de régression",value = FALSE))
                     
                     ),
                     box(
                       title = "Options du graphique",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                       
                       sliderInput("SizePoint", "Taille",min=1, max=10, value=2, step=0.5),
                       conditionalPanel(condition="input.RadioPlotBi=='hist'",sliderInput("binwidthBi", "Taille de la fenètre glissante",min=1, max=100, value=1, step=1)),
                       conditionalPanel(condition="input.RadioPlotBi=='Nuage'",selectInput("ColorBiplot", "Couleur",choices=c("Bleue" = 'blue',"Rouge"='red',"Vert"='green', "Noir"='black'),width="50%")),
                       sliderInput("TransAlphaBi", "Transparence",min=1, max=100, value=50, step=1),
                       conditionalPanel(condition="input.RadioPlotBi!='Nuage'", radioButtons("SensGraphBi","Sens du graph",choices=c("Vertical"="Vert","Horizontal"="Hori"))),
                       conditionalPanel(condition="input.RadioPlotBi=='box'", checkboxInput("CheckAddPointsBoxBi","Ajouter les données",value=FALSE))
                       
                     )
                   
              )
              
              
              
            )
            )
    
    
    
  )
)

# Put them together into a dashboardPage
dashboardPage(skin="blue",
  dashboardHeader(title = "Meta16S"),
  sidebar,
  body
)


