function(request) {
sidebar <- dashboardSidebar(
  useShinyjs(),
  inlineCSS(appCSS),
#   tags$head(
#     tags$script(src = "custom.js")
#   ),
  div(id = "loading-content-bar",
      p()),
  div(
    id = "app-content-bar",
    sidebarMenu(id = "side",
      menuItem("Home", tabName = "Home", icon = icon("home")),
      menuItem("Tutorial", tabName = "Tutorial", icon = icon("book")),
      menuItem("Download/Install", tabName = "Download", icon = icon("download")),
      menuItem("Raw data", tabName = "RawData", icon = icon("upload")),
      menuItem("Upload your data", tabName = "Upload", icon = icon("upload")),
      #bookmarkButton(),
      menuItemOutput("dymMenu"),
      img(src = "logo.jpg", height = 49, width = 220,style="position:absolute;bottom:0;margin:0 0 15px 10px;")
    )
  )
)

body <- dashboardBody(
  tags$style(type="text/css", Errorcss),
  useToastr(),
  useShinyjs(),
  inlineCSS(appCSS),
  div(
    id = "loading-content",
    br(),
    br(),
    br(),
    h2("Please wait while SHAMAN is loading...")),
  div(
    id = "app-content-bar",
  tabItems(
    tabItem(tabName = "Home",
            fluidRow(
              column(width=9,
            div(style="width:100% ; max-width: 1200px; height: 550px",
                
              tabBox(title="Welcome to SHAMAN", id="tabset1", width=NULL,
                   # tags$script(type="text/javascript", language="javascript", src="google-analytics.js"),
                   tabPanel("About", 
                            p("SHAMAN is a SHiny application for Metagenomic ANalysis including the normalization,
                                       the differential analysis and mutiple visualization.",style = "font-family: 'times'; font-si16pt"),
                            p("SHAMAN is based on DESeq2 R package", a("[Anders and Huber 2010]", href="http://www.ncbi.nlm.nih.gov/pubmed/20979621"), "for the analysis of metagenomic data, as suggested in", a("[McMurdie and Holmes 2014,",href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3974642/"),a("Jonsson2016]",href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4727335/"),
                              ". SHAMAN robustly identifies the differential abundant genera with the Generalized Linear Model implemented in DESeq2", a("[Love 2014]", href="http://www.ncbi.nlm.nih.gov/pubmed/25516281"),".
                              Resulting p-values are adjusted according to the Benjamini and Hochberg procedure [Benjamini and Hochberg 1995].
                              The PCOA is performed with the", a("ade4 R package",href="http://pbil.univ-lyon1.fr/ade4/"), "and plots are generated with", a("ggplot2",href="http://ggplot2.org/"), "or", a("D3.js packages",href="http://d3js.org/"), ".
                              A presentation about SHAMAN is available", a("here",target="_blank",href="shaman_presentation.pdf")," and a poster", a("here.",target="_blank",href="shaman_poster.pdf"), br(),
                              "SHAMAN is compatible with standard formats for metagenomic analysis. We also provide a complete pipeline for OTU picking and annotation named",a("MASQUE", href="https://github.com/aghozlane/masque") ,"used in production at Institut Pasteur.",style = "font-family: 'times'; font-si16pt"),
                            p("Hereafter is the global workflow of the SHAMAN application:"),
                            div(img(src = "Workflow_sh.png",width = "100%",height = "100%",style="max-width: 800px;"),Align="center")
                            ),
                   tabPanel("Authors", h3("The main contributors to SHAMAN:"),
                            p(a("Stevenn Volant", href="mailto:stevenn.volant@pasteur.fr"), "(Initiator, coding, testing, documentation, evaluation)"),
                            p(a("Amine Ghozlane",href="mailto:amine.ghozlane@pasteur.fr"), "(Coding, testing, documentation, evaluation, packaging)"),
                            p("Pierre Lechat", "(Coding, testing, feature suggestions)"),
                            h3("Acknowledgements"),
                            p("Thanks to the following people for patches and other suggestions for improvements:"),
                            p("Carine Rey, ", "Hugo Varet,", "Julien Tap, ","Anna Zhukova.")
                          ),
                   tabPanel("Citing SHAMAN",
                   p("If you use SHAMAN for your project, please cite our first application of SHAMAN in Quereda et al. 2016.",style = "font-family: 'times'; font-si16pt"),
                   p("Publication using SHAMAN :",style = "font-family: 'times'; font-si18pt; font-style: strong"),
                   p(a("Diverse laboratory colonies of Aedes aegypti harbor the same adult midgut bacterial microbiome.", href="https://www.ncbi.nlm.nih.gov/pubmed/29587819"), "Dickson LB, Ghozlane A, Volant S, Bouchier C, Ma L, Vega-Rúa A, Dusfour I, Jiolle D, Paupy C, Mayanja MN, Kohl A, Lutwama JJ, Duong V, Lambrechts L; Parasit Vectors 2018",style = "font-family: 'times'; font-si16pt"),
                   p(a("Characteristics of Fecal Microbiota in Pediatric Crohn’s Disease and Their Dynamic Changes During Infliximab Therapy.", href="https://www.ncbi.nlm.nih.gov/pubmed/29194468"), "Wang Y, Gao X, Ghozlane A, Hu H, Li X, Xiao Y, Li D, Yu G, Zhang T; Journal of Crohn's & colitis  2017",style = "font-family: 'times'; font-si16pt"),
                   p(a("Carryover effects of larval exposure to different environmental bacteria drive adult trait variation in a mosquito vector.", href="https://www.ncbi.nlm.nih.gov/pubmed/28835919"), "Dickson LB, Jiolle D, Minard G, Moltini-Conclois I, Volant S, Ghozlane A, Bouchier C, Ayala D, Paupy C, Moro CV, Lambrechts L; Science Advances 2017",style = "font-family: 'times'; font-si16pt"),
                   p(a("A bacteriocin from epidemic Listeria strains alters the host intestinal microbiota to favor infection.", href="http://www.ncbi.nlm.nih.gov/pubmed/27140611"), "Quereda JJ, Dussurget O, Nahori MA, Ghozlane A, Volant S, Dillies MA, Regnault B, Kennedy S, Mondot S, Villoing B, Cossart P, Pizarro-Cerda J.; PNAS 2016",style = "font-family: 'times'; font-si16pt"),
                   p("Reporting bugs, ask for help",style = "font-family: 'times'; font-si18pt; font-style: strong"),
                   p("If you have any comments, questions or suggestions, or need help to use SHAMAN, please contact us at", a("shaman@pasteur.fr", href="mailto:shaman@pasteur.fr"),"and please provide us with enough information that we can recreate the problem. Useful things to include are:", style = "font-family: 'times'; font-si16pt;"),
                   tags$ul(
                     tags$li("Input data (or examples, a small test case sufficient to recreate the problem)"), 
                     tags$li("Information about which system your are using: web version, docker or R installation")
                   )
            )))),
              column(width=3,
            box(
              title = "What's new in SHAMAN", width = NULL, status = "primary",
              div(style = 'overflow-y: scroll; height: 550px',
                  addNews("April 17th 2018","Bioinformatics","The bioinformatic treatment offers a larger access to parameters. We also worked a lot on the documentation."),
                  addNews("September 4th 2017","Bioinformatics","The bioinformatic treatment of fastq reads is now available in SHAMAN. For now, SHAMAN allows to compute OTU, build an OTU table and annotate them with the last version of the available database. This application is for 16S/18S/23S/28S/ITS sequencing."),
                  addNews("July 18th 2017","Normalization and visualisation","A new method for normalization called total counts was added. More options have been added to the abundance tree."),
                  addNews("May 30th 2017","Bug fixes","Some visualization bug with the abundance tree and phylogenetic tree are now fixed. The export of the relative abundance and normalised abundance are now given in the right level. This update prepares the field for the next major release of shaman for June."),
                  addNews("March 30th 2017","Krona, Phylogeny and bug fixes","Krona and phylogenetic tree plots are now available in visualisation. Several new distance are available in PCOA. The import float count matrices is now ok. We have finaly debugged the export of the relative abundance/normalized matrices."),
                  addNews("Dec 9th 2016","Phylogenetic tree and stress plot","You can now upload a phylogenetic tree to calculate the unifrac distance (only available at the OTU level). 
                          The stress plot has been added to evaluate the goodness of fit of the NMDS."),
                  addNews("Nov 22th 2016","New visualization and bug fix","We have implemented a new visualization called tree abundance. Some bugs have been fixed (thanks Carine Rey from ENS)."),
                  addNews("Oct 12th 2016","Filtering step and bugs fix","You can now apply a filter on the features according to their abundance 
                          and the number of samples. Bugs on confidence intervals for the alpha diversity have been fixed."),
                  addNews("Sep 21th 2016","SHAMAN on docker","The install of SHAMAN is now available with docker.
                           The R install is also updated and passed in release candidate state."),
                  addNews("Sep 14th 2016","Download and install SHAMAN","You can install SHAMAN (beta)."),
                  addNews("Sep 9th 2016","PCA/PCOA","You can select the axes for the PCOA and PCA plots."),
                  addNews("Aug 1st 2016","Biom format","SHAMAN can now support all the Biom format versions."),
                  addNews("Jun 24th 2016","Comparisons plots","The venn diagram and the heatmap of foldchange 
                                                                have been added to compare the results of 2 or more contrasts."), 
                  addNews("Jun 17th 2016","Diversity plots","Enhancement of the visualtisation of the diverties. 
                                                              The shanon and inv. shanon have been added.")
                  )
            )
              )
              )
    ),
    tabItem(tabName = "Tutorial",
            div(style="width:100% ; max-width: 1200px",
            tabBox(title="How to use SHAMAN", id="tabset1", width =NULL,
            tabPanel("Introduction",
            p("First check out our user guide:"),
            tags$ul(
             tags$li(a("User's guide",target="_blank", href="Userguide.pdf")), 
             tags$li(a("Writing a target file for the experimental design",target="_blank", href="experimental_design_guide.pdf")) 
            ),
            p(" You can test SHAMAN with the dataset from", a("[Tap et al. 2015]",href="http://www.ncbi.nlm.nih.gov/pubmed/26235304"),
              ", which is available", a("here",target="_blank",href="Alimintest.zip"),"."),
            p("The zip archive contains 4 different files :", br(),
              "- the otu count matrix : Alimintest_otu_table.csv,", br(),
              "- the otu annotation table : Alimintest_otu_annotation.csv,", br(),
              "- the experimental design : Alimintest_target.csv,", br(),
              "- the contrast table : Alimintest_contrasts.csv."),
            p("Two groups of person follow two  strict diet periods that involve the intake of 40g following 10g of fiber per day, or 10g of fiber after a 40g fiber intake period:"),
            img(src = "tutorial/FigS1.png",width = "100%",style="max-width: 900px"),
            p("The 16S rRNA (V3 - V4 regions) from fece samples was sequenced at time stamp : 2, 3, 4 and 5.", br(),
              "The analysis will consider the different impact of the different fiber intake and the comparison to patient metabolic data.")),
            tabPanel("1-Load 16S data",
            p("The first step consists to load the count table and the annotation table as follow :"),
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
            tabPanel("4-Tables",
                     p("'Tables' section provides the results of the differential analysis.
                       For one given contrast, we have:",br(),
                       "- The id of the given taxonomical level",br(),
                       "- The base mean is the mean normalized count for the given annotation of all samples.",br(),
                       "- The fold-change is a mesure describing how much the abundance varies from one condition to an other. For instance, if the abundance is 100 in condition 1 and 200 for condition 2, the fold-change equals 100/200=0.5.",br(),
                       "- The log2 fold-change is the log2 value of the fold-change.",br(),
                       "- The p-value adjusted (padj) is the pvalue obtained by the Wald test and adjusted by Benjamini & Hochberg procedure (BH) or Benjamini & Yekutieli procedure (see linear model options).",
            img(src = "tutorial/tutorial_table.png",width = "100%",style="max-width: 900px"))),
            tabPanel("5-Visualization",
                     p("'Diagnostic plots' section provides several visualization to control the analysis",br(),
                     "- Barplot",br(),img(src="tutorial/tutorial_barplot.png",width = "100%",style="max-width: 900px"),hr(),
                     "- Heatmap",br(),img(src="tutorial/tutorial_heatmap.png",width = "100%",style="max-width: 900px"),hr(),
                     "- Boxplot",br(),img(src="tutorial/tutorial_boxplot.png",width = "100%",style="max-width: 900px"),hr(),
                     "- Diversity",br(),img(src="tutorial/tutorial_diversity.png",width = "100%",style="max-width: 900px"),hr(),
                     "- Rarefaction",br(),img(src="tutorial/tutorial_rarefaction.png",width = "100%",style="max-width: 900px"))),
            tabPanel("6-Other datasets available",
            p("You can test SHAMAN with two other datasets : "),
            p("- Bacteriocin impact on mice microbiota from ", a("[Quereda et al. 2016]",href="http://www.ncbi.nlm.nih.gov/pubmed/27140611"), ":", a("here",target="_blank",href="listeria.zip"), br(),br(),
              
              "The project is divided into three group :  WT, Delta : Bacteriocin KO, Delta-complemented: Bacteriocin KO + complementation.", br(),
              "Three time points are considered :", br(),
              "- T0 : merged counts from -48H and -24H", br(),
              "- T2 : 6H after listeria infection", br(),
              "- T3 : 24H after listeria infection", br(),
              
              img(src="listeria.png",width = "100%",style="max-width: 500px"),br(),br(),
              
              "Analysis with SHAMAN:", br(),
              "Set variables: condition, time, mice", br(),
              "Set interactions:",
              "1. condition:mice", br(),
              "2. condition:time", br(),hr(),
              "- MOCK communities from ", a("[The NIH HMP Working Group, 2009]",href="http://www.ncbi.nlm.nih.gov/pubmed/19819907"), ":", a("here",target="_blank",href="mock.zip"),br(),
              "Mock communities are composed of 21 species mixed in even or staggered proportions :", br(),br(),
              img(src="mock.png",width = "100%",style="max-width: 600px"),br(),br(),
              "Analysis with SHAMAN: ",br(),
              "Set variables : Community", br(),
              "no interaction"))
              ))
    ),
    
    tabItem(tabName = "Download",
            fluidRow(
              column(width=9,
                     div(style="width:100% ; max-width: 1200px",
                         
                         tabBox(title="Download / Install SHAMAN", id="tabset1", width=NULL,
                                tabPanel("Docker install",
                                p("Docker install is the easiest way to use SHAMAN locally."),
                                p("- Install docker on ubuntu (Linux):",a("here",href="https://docs.docker.com/engine/installation/linux/ubuntulinux/")),
                                p("- Install docker on Windows and Mac:"),
                                p("Download and install docker from",a("https://www.docker.com/", href="https://www.docker.com/"),
                                  "Then, you will need to install the", a("Docker toolbox.", href="https://www.docker.com/products/docker-toolbox"),
                                  "Once installed, run ‘Docker Quickstart Terminal’ application."),
                                p("- Running SHAMAN:"),
                                wellPanel(div(style = 'max-width: 900px',"docker pull aghozlane/shaman && docker run --rm -p 80:80 -p 5438:5438 aghozlane/shaman")),
                                p("Then connect to http://0.0.0.0/ with your favorite web navigator."),
                                p("Failed: port is already allocated ?"),
                                wellPanel(div(style = 'max-width: 900px',"docker run --rm -p 3838:80 -p 5438:5438 aghozlane/shaman")),
                                p("Then connect to http://0.0.0.0:3838/."),
                                p("- Docker update after an update of SHAMAN:"),
                                wellPanel(div(style = 'max-width: 900px',"docker pull aghozlane/shaman"))
                                ),
                                tabPanel("R install (RC)",
                                         p("SHAMAN is available for R>3.1.2. The installation, download and execution can all be performed with a small R script :",style = "font-family: 'times'; font-si16pt"),
                                         wellPanel(style = 'width: 50%; word-wrap: break-word;',"# Load shiny packages",br(),
                                           "if(!require('shiny')){",br(),"  install.packages('shiny')",br(),"  library(shiny)",br(),"}",br(),
                                           "system(\"Rscript -e 'library(\\\"shiny\\\");runGitHub(\\\"pierreLec/KronaRShy\\\",port=5438)'\",wait=FALSE)",
                                              br(),"# Install dependencies,",br(),"# download last version from github,",br(),"# and run SHAMAN in one command:",br(),"runGitHub('aghozlane/shaman')"),
                                         p("This script can also be dowloaded", a("here", target="_blank", href="shamanapp.R"), "and executed as following :"),
                                         wellPanel(style = 'width: 50%; word-wrap: break-word;',"chmod +x ./shamanapp.R && Rscript ./shamanapp.R")),
                                         p("Of note, the R version has an impact on the contrast definition. DESeq2 contrast are harder to define on R>3.2."),
                                         p("Contribution to SHAMAN code are always welcome and can be performed with the", a("github deposit.",href="https://github.com/aghozlane/shaman")
                                ))
                         )
                     )
              )
            ),
    #id="rawdatatab",
    tabItem(tabName = "RawData",
            tags$style(type='text/css', ".well { max-width: 20em; }"),
            # tags$head(tags$style(HTML(InfoBoxCSS))),
            tags$head(tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })")),
            
            fluidRow(
              # column(width=3,valueBoxOutput("valueErrorPercent",width=NULL)),
              # column(width=3,infoBoxOutput("InfoErrorCounts",width=NULL)),
              # column(width=3,infoBoxOutput("InfoErrorTaxo",width=NULL)),
              div(id="masque-infobox",
              column(width=3,
                     withPopup(infoBoxOutput("infoBoxPass",width=NULL),
                                      title="Once you have entered a valid email address, click on this button:",
                                      img_src="icons/GetKey_button.png")
                     ),
              column(width=3,
                    infoBoxOutput("infoBoxFastQ",width=NULL)
                    ),
              column(width=3,      
                     infoBoxOutput("infoBoxFastQ_match",width=NULL)
              ),
              column(width=3,      
                     valueBoxOutput("progressBoxMasque",width=NULL)
              )
            )),
            
                        # fluidRow(
                          # HTML('<center><h1>MASQUE : Metagenomic Analysis with a Quantitative pipeline</h1></center>'),
            #               column(width=8,
            #                 div(style="background-color: white; padding-left:20px;",
            #                 HTML('
            #                <h2>Introduction</h2>
            # 
            #                      <p>The aim of this part is to provide an easy cluster interface to perform <b>targeted metagenomic analysis</b>. This analysis will be done using the MASQUE pipeline which allows :</p>
            # 
            #                      <ul>
            #                      <li>to analyse 16S/18S/23S/28S/ITS data. It builds a count matrix, an annotation table and a phylogeny of the OTU.</li>
            #                      <li>to perform to use a set of parameters already tested on serveral projects for the numerous software used to perform the clustering and the annotation.</li>
            #                      <li>to perform an "uptodate" analysis considering the scientific litterature.</li>
            #                      </ul>
            # 
            #                       <hr></hr>
            #                      <h2>Process</h2>
            # 
            #                      <p>We follow the recommandation described by Robert C. Edgar in <a href="http://www.nature.com/nmeth/journal/v10/n10/full/nmeth.2604.html" target="_blank" >Uparse</a> supplementary paper.<br>
            #                        The clustering process in MASQUE is performed as the following :<br>
            #                        1.  Read quality control<br>
            #                        2.  Dereplication<br>
            #                        3.  Chimera filtering<br>
            #                        4.  Clustering<br>
            #                        5.  Realignment/mapping<br>
            #                        6.  Taxonomical annotation of the OTU<br>
            #                        7.  Quality check of every step  </p>
            # 
            #                        <p>You can find more information in the presentation <a href="/aghozlane/masque/blob/master/tp/Targeted_metagenomics.pdf">here</a>. We try to describe the idea behind each step and a complete TP to do it on your own.</p>'
            #                   )
            # #
            # #
            # #                      <\ui> to analyse 16S/18S/23S/28S/ITS data. It builds a count matrix, an annotation table and a phylogeny of the OTU.
            # #                      to perform to use a set of parameters already tested on serveral projects for the numerous software used to perform the clustering and the annotation.
            # #                      to perform an uptodate analysis considering the scientific litterature."
            #               )),
            # column(width=4)
            # # column(width=4,
            # #         inlineCSS(gaugeCSS),
            # #         gaugeOutput("gaugeMasque", width = "100%", height = "100%")
            # #        )
            # ),
            
            # hr(),
            # HTML('<center><h2 style="color:#053383;"><b>Start your analysis</b></h2></center>'),
            # hr(),
            fluidRow(

              
              # column(width=3,
              #        inlineCSS(gaugeCSS),
              #        gaugeOutput("gaugeMasque", width = "100%", height = "100%")
              # 
              #        ),
              column(width=9,
                     
              div(id="masque-form",
                box(title="About",width = NULL, status = "primary",
  
                  column(width=6,  radioButtons("DataTypeMasque",label = "Type of data",choices = c("16S" = "16S", "18S" = "18S","23S/28S" = "23S_28S","ITS" = "ITS"),inline = TRUE)),
                  column(width=6,  radioButtons("PairedOrNot",label = "Paired-end sequencing ?",choices = c('Yes'="y","No"="n"),selected = "n",inline = TRUE)),
                  column(width=3,  textInput("to", "Email address *", value="yyy@xxx"),
                                    bsTooltip("to", "Enter a valid email address to get your results","bottom",trigger = "hover", options = list(container = "body"))
                         ),
                  column(width=3,  actionButton("checkMail",label = "  Get key",icon=icon('key')),
                                   bsTooltip("checkMail", "Click here to get your key by mail","bottom", options = list(container = "body")),
                                   tags$style(type='text/css', "#checkMail { width:100%; margin-top: 25px;}")
                         ),
                  column(width=6,  selectizeInput("HostName",label = "Select the host",
                                                  choices =   c("None"="", "a.stephensi (mosquito)" = "astephensi"," b.taurus (cow)" = "btaurus", "c.familiaris (dog)" = "cfamiliaris",
                                                                "chiroptera (bat)" = "chiroptera", "c.sabaeus (Apes)" = "csabaeus" , "d.melanogaster (fly)"="dmelanogaster", 
                                                                "e.caballus (Horse)" = "ecaballus", "f.catus (cat)" = "fcatus", "hg18 (human)"="hg18", "hg19 (human)"="hg19", 
                                                                "hg38 (human)" = "hg38", "m.lucifugus (bat)" = "mlucifugus", "mm8 (mouse)"= "mm8", "mm9 (mouse)" = "mm9", 
                                                                "mm10 (mouse)"="mm10", "p.vampyrus (bat)" = "pvampyrus", "s.scrofa (Boar)" = "sscrofa"))),
                  column(width=6,checkboxInput("primer", "Specify the primer"),
                                  bsTooltip("primer", "If no, default values are chosen","bottom", options = list(container = "body"))
                         ),
                  column(width=6,checkboxInput("options", "More workflow options"),
                         bsTooltip("options", "If no, default values are chosen","bottom", options = list(container = "body"))
                  ),
                  conditionalPanel(condition="input.primer && input.PairedOrNot=='y'",
                                   column(width=12,
                                     textInput("R1primer",label = "Forward primer",value = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG"),
                                     textInput("R2primer",label = "Reverse primer",value = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC")
                                   )
                  ),
                  conditionalPanel(condition="input.primer && input.PairedOrNot=='n'",
                                   column(width=12,
                                     textInput("primerSingle",label = "Primer",value = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG")
                                   )
                  ),
                  conditionalPanel(condition="input.options",
                                   column(width=4, h3("Read processing"),
                                          sliderInput("phredthres", "Phred quality score cutoff to trim off low-quality read ends",
                                                      min = 0, max = 40, value = 20),
                                          sliderInput("mincorrect", "Minimum allowed percentage of correctly called nucleotides per reads",
                                                      min = 50, max = 100, value = 80),
                                          numericInput("minreadlength", "Minimum read length", 50,step=1,min=50),
                                          conditionalPanel(condition="input.options && input.PairedOrNot=='y'", numericInput("minoverlap", "Minimum overlap size", 10,step=1,min=1))
                                          
                                   ),
                                   column(width=4, h3("OTU processing"),
                                          selectInput("dreptype", "Dereplication ", choices = list("Prefix" = "--derep_prefix", "Full length" = "--derep_fulllength"), selected = "--derep_prefix"),
                                          numericInput("maxampliconlength", "Maximum OTU length (0 is no limit)", 0,step=1,min=0),
                                          numericInput("minampliconlength", "Minimum OTU length", 50, step=1, min=35),
                                          numericInput("minabundance", "Minimum abundance at dereplication", 4, step=1, min=2),
                                          selectInput("clusteringstrand", "Clustering strand ", choices = list("Both" = "both", "Plus" = "plus"), selected = "both"),
                                          sliderInput("clusteringthreshold", "Clustering threshold", min = 0, max = 1, value = 0.97)),
                                   column(width=4, h3("OTU annotation"),
                                          selectInput("annotationstrand", "Annotation strand ", choices = list("Both" = "both", "Plus" = "plus"), selected = "both"),
                                          sliderInput("annotationKingdomthreshold", "Minimum identity for Kingdom annotation", min = 0, max = 1, value = 0.75, step = 0.005),
                                          sliderInput("annotationPhylumthreshold", "Identity thresholds for Phylum annotation", min = 0, max = 1, value = c(0.75, 0.785), step = 0.005),
                                          sliderInput("annotationClassthreshold", "Identity thresholds for Class annotation", min = 0, max = 1, value = c(0.785, 0.82), step = 0.005),
                                          sliderInput("annotationOrderthreshold", "Identity thresholds for Order annotation", min = 0, max = 1, value = c(0.82, 0.865), step = 0.005),
                                          sliderInput("annotationFamilythreshold", "Identity thresholds for Family annotation", min = 0, max = 1, value = c(0.865, 0.945), step = 0.005),
                                          sliderInput("annotationGenusthreshold", "Identity thresholds for Genus annotation", min = 0, max = 1, value = c(0.945, 0.98), step = 0.005),
                                          sliderInput("annotationSpeciethreshold", "Minimum identity for Specie annotation", min = 0, max = 1, value = 0.98, step = 0.005)
                                  )
                  )
                ),

                box(title="Directory containing the FastQ files ",width = NULL, status = "primary",
                    # column(width=12,verbatimTextOutput("dirSel")),
                    # br(),
                    column(width=12,
                           
                           fileInput("dir",label = 'Select your fastq files',accept = c(".fastq,.fastq.gz,.fgz,.gz"),multiple = TRUE),
                            # shinyDirButton("dir", "Select a directory", "Upload",buttonType = "primary")
                          
                           # tags$input(id = "dir2", webkitdirectory = TRUE, type = "file", onchange="pressed()"),
                            # HTML("&nbsp;"),
                            # actionButton("LoadFiles",'Load',icon=icon("play"))
                           uiOutput("FastQList_out")
                    )

                    
                  ),
                  box(id="box-match",title=" Match the paired files (only for paired-end sequencing)",width = NULL, status = "primary",
                          column(width=6,  textInput("R1files",label = "Suffix R1 (Forward, do not include file extension)",value = "_R1"),
                                selectInput("R1filesList",label = "","",multiple =TRUE,selectize=FALSE)),
                                column(width=6,  
                                        textInput("R2files",label = "Suffix R2 (Reverse, do not include file extension)",value = "_R2"),
                                        selectInput("R2filesList",label = "","",multiple =TRUE,selectize=FALSE)
                                       ),
                                column(width=12,
                                        actionButton("MatchFiles_button",'Match',icon=icon("exchange"),class="btn-primary",style="color: #fff;"),
                                        HTML("&nbsp;"),
                                        actionButton("RemoveFastQbut_R1R2",'Remove file(s)',icon=icon("remove"))
                                      )
                  ),
                div(style = "text-align:center;",
                    actionButton("submit", strong("Check and submit"), icon("chevron-circle-right"),class="btn-primary",style = "color: #fff"),
                    receiveSweetAlert(messageId = "ErrorMasque"),
                    tags$style(type='text/css', "#submit { width:50%; margin-top: 5px;}"),
                    receiveSweetAlert(messageId = "SuccessMasque")
                )
                
              ),
              div(id="project_over",style="text-align:center; display: none;",
                  HTML("<center><h1><strong>Computations are over</strong></h1></center>"),
                  actionButton("Check_project_over", h2("Check results"),icon=icon("check-circle fa-2x"),class="btn-primary",style = "color: #fff"),
                  tags$style(type='text/css', "#Check_project_over {margin-top: 15px;width:50%;}")
              ),
              
              div(id="current-project",style="display: none;",
                  uiOutput("masque_results")
              ),
              div(id="project-over-wait",style="display: none;",
                  HTML('<center><h1><strong> Please wait during the creation of the files...</strong></h1> <br/> <em><h4> This can take about 30 seconds</h4> </em> </center>'),
                  tags$img(src = "gears.gif",id ="loading-spinner")
              )),
              div(id="reload-project",style="display: none;",
                  column(width=12,
                  actionButton("comeback",strong("Close and come back"),icon=icon("backward")),
                  #shinyjs::useShinyjs(),
                  #shinyjs::extendShinyjs(text = "shinyjs.refresh = function() { history.go(0); }"),
                  actionButton("refresh", "Refresh progress"),
                  uiOutput("masque_status_key")
                  
                  )
              ),
              # div(id="MasqueToShaman",style="display: none;",
              #     column(width=3,
              #     #uiOutput("loaddb"),
              #     box(id="box-zip",title="Download .zip file",width = NULL, status = "success",
              #         downloadButton('Download_masque_zip', 'Download the results'),
              #         tags$style(type='text/css', "#Download_masque_zip { width:100%;}")
              #     )
              #     )
              # ),
              div(id="pass",style = "word-wrap: break-word;",
              column(width=3,
                     box(id="boxpass",title = strong("Enter the key"), width = NULL, background = "light-blue",
                         inlineCSS(list(.pwdGREEN = "background-color: #DDF0B3",.pwdRED = "background-color: #F0B2AD")),
                         uiOutput("pass_Arg"),
                         textInput("password","",value = NULL),
                         div(style = "text-align:right;",
                             actionButton("Check_project", "Check project",icon=icon("check-circle")),
                             tags$style(type='text/css', "#Check_project {margin-top: 15px;}")
                         ),
                         bsTooltip("password", 'Fill in the email field and click the "Get key" button',"bottom",trigger = "hover", options = list(container = "body")),
                         tags$style(type='text/css', "#password { width:100%; margin-top: 10px;}")
                     ),
                     tags$style(type='text/css', "#boxpass {margin-bottom: 0px;}"),
                     # valueBoxOutput("progressBoxMasque",width = 12),
                     # infoBox(title = "tete",icon= uiOutput("spinner_icon"),width=12),
                     uiOutput("summary_box_masque")
              )),
                 column(width=3,
                #        inlineCSS(gaugeCSS),
                #        gaugeOutput("gaugeMasque", width = "100%", height = "100%"),
                      uiOutput("InfoMasque"),
                      uiOutput("InfoMasqueHowTo")
                      # uiOutput("MasqueToShaman_button")
                        )
            )
            
    ),

    tabItem(tabName = "Upload",
            tags$style(type='text/css', ".well { max-width: 20em; }"),
            fluidRow(
              column(width=3,valueBoxOutput("valueErrorPercent",width=NULL)),
              column(width=3,infoBoxOutput("InfoErrorCounts",width=NULL)),
              column(width=3,infoBoxOutput("InfoErrorTaxo",width=NULL)),
              column(width=3,infoBoxOutput("InfoTreePhylo_box",width=NULL))
            ),
            br(),
             fluidRow(
                box(title="Select your file format",width = 3,status = "success", solidHeader = TRUE,collapsible = FALSE,
                  # selectInput("FileFormat","",c("Count table & taxonomy (*.csv or *.tsv)"="fileCounts","BIOM file"="fileBiom","Saved project"="fileRData"),selected="fileCounts"),
                  selectInput("FileFormat","",c("Count table & taxonomy (*.csv or *.tsv)"="fileCounts","BIOM file"="fileBiom","Project number"="projnum"),selected="fileCounts"),
                  conditionalPanel(condition="input.FileFormat=='fileCounts'",
                                   checkboxInput("NoTaxoFile","No taxonomy table",value=FALSE)
                  ),
                  conditionalPanel(condition="input.FileFormat=='projnum'",
                                 inlineCSS(list(.pwdGREEN = "background-color: #DDF0B3",.pwdRED = "background-color: #F0B2AD")),
                                 textInput("password_home","Enter the key",value = NULL),
                                 div(style = "text-align:right;",
                                     actionButton("Check_project_home", "Check project number",icon=icon("check-circle")),
                                     tags$style(type='text/css', "#Check_project_home {margin-top: 15px;}")
                                 )
                    
                  )
                         
                ),
                
                conditionalPanel(condition="input.FileFormat=='fileCounts'",
                  box(title="Load the count table",width = 3,height = "260px", status = "primary", solidHeader = TRUE,collapsible = FALSE,
                      fluidRow(
                        column(width=6,radioButtons("TypeTable",h6(strong("Type:")),c("OTU/Gene table"="OTU","MGS table"="MGS"))),
                        column(width=6,selectInput("sepcount", h6(strong("Separator:")),c("Tab" = "\t", "Comma" = ",","Semi-colon" = ";")))
                      ),
                      fileInput('fileCounts', h6(strong('Select your file')),width="100%",accept = c(".csv",".tsv",'.xls','.xlsx','.txt')),
                      tags$script('$( "#fileCounts" ).on( "click", function() { this.value = null; });')
                  ),
                  conditionalPanel(condition="!input.NoTaxoFile",
                    box(title="Load the taxonomy file",width = 3,height = "260px", status = "primary", solidHeader = TRUE,collapsible = FALSE,
                        fluidRow(
                          column(width=6,radioButtons("TypeTaxo",h6(strong("Format:")),c("Table"="Table","RDP"="RDP"))),
                          column(width=6,
                               conditionalPanel(condition="input.TypeTaxo=='RDP'",numericInput("RDP_th",h6(strong("Threshold:")),0.5,step=0.01,min=0.01,max=1))
                          ),
                          column(width=6,
                                 conditionalPanel(condition="input.TypeTaxo=='Table'",selectInput("septaxo", h6(strong("Separator:")),
                                                      c("Tab" = "\t", "Comma" = ",", "Semicolon" = ";")))
                                 )
                        ),
                        fileInput('fileTaxo', h6(strong('Select your file')),width="100%",accept = c(".csv",".tsv",'.xls','.xlsx','.txt')),
                        tags$script('$( "#fileTaxo" ).on( "click", function() { this.value = null; });')
                    )
                  )
                ),
                
                conditionalPanel(condition="input.FileFormat=='fileBiom'",
                                 box(title="Load the BIOM file",width = 3, status = "primary", solidHeader = TRUE,collapsible = FALSE,
                                     fileInput('fileBiom', h5(strong('Select your file')),width="100%",accept = c(".biom")),
                                     tags$script('$( "#fileBiom" ).on( "click", function() { this.value = null; });')
                                 )
                ),
                
                conditionalPanel(condition="input.FileFormat=='projnum'",
                                uiOutput("Project_box_home")
                ),
                
                conditionalPanel(condition="input.FileFormat!='projnum'",
                  box(title="Load phylogenetic tree (optional)",width = 3, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = TRUE,
                      fileInput('fileTree', h6(strong('Select your file (tree)')),width="100%"),
                      tags$script('$( "#fileTree" ).on( "click", function() { this.value = null; });')
                  )
                ),
                
                fluidRow(column(width=3,
                                uiOutput("InfoCountsFile"),
                                uiOutput("InfoTaxoFile"),
                                uiOutput("InfoBIOM")
                          )
                )
  
                
             ),
              column(id="tabboxdata_col",width=12,uiOutput("TabBoxData")),
              receiveSweetAlert(messageId = "ErrorTaxo"),
              receiveSweetAlert(messageId = "ErrorBiom1"),
              receiveSweetAlert(messageId = "ErrorBiom2"),
              receiveSweetAlert(messageId = "ErrorSizeFactor"),
              receiveSweetAlert(messageId = "ErrorCounts"),
              receiveSweetAlert(messageId = "ErrorRDP")
    ),
    
  #### Statistical analysis

    tabItem(tabName = "RunDiff",
            fluidRow(
              column(width=3,valueBoxOutput("RowTarget",width=NULL)),
              column(width=3,infoBoxOutput("InfoTaxo",width=NULL)),
              column(width=3,infoBoxOutput("InfoDESeq",width=NULL)),
              column(width=3,infoBoxOutput("InfoContrast",width=NULL))
            ),            
            fluidRow(
              column(width=5,
                box(title="Experimental design",width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = FALSE,
                  fluidRow(
                    column(width=5,fileInput('fileTarget', h6(strong('Select your target file')),width="100%"),
                           tags$script('$( "#fileTarget" ).on( "click", function() { this.value = null; });')
                           ),
                    column(width=3,selectInput("septarget", h6(strong("Separator:")), c("Tab" = "\t", "Comma" = ",", "Semicolon" = ";"))),
                    column(width=4,uiOutput("SelectTaxo"))
                  ),
                  fluidRow( 
                    column(width=6,uiOutput("SelectInterestVar")),
                    column(width=6,uiOutput("SelectInteraction2"))
                  ),
                  fluidRow( 
                    column(width=6,uiOutput("RunButton"))
                  )
                ),
                fluidRow(uiOutput("InfoModel"),uiOutput("InfoModelHowTo")),
                uiOutput("BoxTarget"),
                uiOutput("BoxCountsMerge")
              ),
              
              column(width=7,
                  box(title="Options",width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed = TRUE,
                    tabBox(title="", id="tabsetOption", width=NULL,
                           tabPanel("Statistical model", 
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
                                             radioButtons("fitType",h6(strong("Relationship")),choices = c("Parametric"="parametric","Local"="local", "Mean"="mean"))
                                      )
                                      # column(width=3,uiOutput("RefSelect"))
                                    )
                           ),
                           tabPanel("Normalization", 
                                    fluidRow(
                                      column(width=3,
                                             selectizeInput("AccountForNA", h6(strong("Normalization method")),choices = c("Usual"="All", "Remove null counts"="NonNull", "Weighted"="Weighted","Total counts"="Total counts"),selected = "NonNull")),
                                      column(width=3,
                                             uiOutput("SelectVarNorm")),
                                      column(width=3,
                                             fileInput('fileSizeFactors', h6(strong('Define your own size factors')),width="100%"),
                                             tags$script('$( "#fileSizeFactors" ).on( "click", function() { this.value = null; });')
                                      ),
                                      column(width=3, selectInput("sepsize", h6(strong("Separator:")), c("Tab" = "\t", "," = "Comma", "Semicolon" = ";"))),
                                      column(width=3,br(),htmlOutput("InfoSizeFactor"))
                                    )
                           ),
                           tabPanel("Filtering",
                                    fluidRow(
                                          column(width=3,
                                                  checkboxInput("AddFilter","Add a filtering step",value = FALSE),
                                                  bsTooltip("AddFilter", "If your count matrix is very sparse, you can add a filter on your data","bottom",trigger = "hover", options = list(container = "body"))
                                          )
                                    ),
                                    fluidRow(
                                      conditionalPanel(condition="input.AddFilter",
                                              column(width=6, 
                                                     plotOutput("Plot_ThSamp"),
                                                     column(width=10,uiOutput("ThSamp"),offset = 1)
                                                     ),
                                              column(width=6,
                                                     plotOutput("Plot_ThAb"),
                                                     column(width=10,uiOutput("ThAb"),offset = 1)
                                                     ),
                                              column(width=12,
                                                     scatterD3Output("Plot_Scatter_Filter")
                                              )
                                      )
                                    )
                           )
                    )
                  ),
                  fluidRow(
                    column(width=8,
                            uiOutput("contrastBox"),
                           uiOutput("contrastBoxAdvanced")
                           ),
                    column(width=4,
                      
                           uiOutput("contrastDefined"),
                           uiOutput("InfoContrast_box")
                      )
                  )
                  
              )
            )
              
            
    ),
    tabItem(tabName = "DiagPlotTab",
            fluidRow(
              column(width=9,
                     tags$head(tags$style(HTML(spincss))),
                     div(id = "plot-container",
                         tags$img(src = "gears.gif",id ="loading-spinner"),
                         plotOutput("PlotDiag",height="100%")
                     ),
                br(),
                conditionalPanel(condition="input.DiagPlot=='SfactorsVStot'",
                  box(title = "Size factors",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                    DT::dataTableOutput("SizeFactTable"),
                    fluidRow( 
                      column(width=3,downloadButton('ExportSizeFactor', 'Export table')),
                      column(width=3,selectInput("sepsizef", h6(strong("Separator:")), c("Tab" = "\t", "Comma" = ",", "Semicolon" = ";")))
                    ),
                    tags$style(type='text/css', "#ExportSizeFactor { margin-top: 37px;}")
                  )
                ),
                  fluidRow(
                conditionalPanel(condition="input.DiagPlot=='pcaPlot'",
                                 box(title = "Eigen values",  width = 6, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= FALSE,
                                    plotOutput("PlotEigen",height="100%")
                                 )
                ),
                
                conditionalPanel(condition="input.DiagPlot=='pcoaPlot'",
                                 box(title = "Eigen values",  width = 6, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= FALSE,
                                     plotOutput("PlotpcoaEigen",height="100%")
                                 )
                ),
                conditionalPanel(condition="input.DiagPlot=='nmdsPlot'",
                                 box(title = "Stress plot",  width = 6, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= FALSE,
                                     plotOutput("PlotnmdsStress",height="100%")
                                 )
                ),
                conditionalPanel(condition="input.DiagPlot=='pcoaPlot' || input.DiagPlot=='nmdsPlot'",
                                 uiOutput("ResPermaTestBox")
                )
                  )
                
              ),
              column(width=3,
                box(title = "Select your plot",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = FALSE,collapsed= FALSE,
                  selectInput("DiagPlot","",c("Total barplot"="barplotTot","Nul barplot"="barplotNul",
                                              "Maj. taxonomy"="MajTax","Boxplots" = "boxplotNorm", "Density"="densityPlot", "Dispersion" = "DispPlot",
                                              "Size factors VS total"="SfactorsVStot", "PCA"="pcaPlot", "PCoA"="pcoaPlot","NMDS"="nmdsPlot","Clustering" = "clustPlot"))
                    ),
                box(title = "Options",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= FALSE,
                    conditionalPanel(condition="input.DiagPlot!='clustPlot' && input.DiagPlot!='pcaPlot' && input.DiagPlot!='SfactorsVStot' && input.DiagPlot!='DispPlot'",
                                     radioButtons("CountsType","Counts:",c("Normalized"="Normalized","Raw"="Raw"),inline = TRUE)
                                     ),
                    conditionalPanel(condition="input.DiagPlot=='boxplotNorm'",
                                    checkboxInput("RemoveNullValue","Remove 0",value = TRUE)
                                    ),
                    conditionalPanel(condition="input.DiagPlot!='Sfactors' && input.DiagPlot!='SfactorsVStot' ",uiOutput("VarIntDiag")),
                    conditionalPanel(condition="input.DiagPlot=='pcoaPlot' || input.DiagPlot=='pcaPlot' || input.DiagPlot=='nmdsPlot'",
                                     h5(strong("Select the modalities")),
                                     uiOutput("ModMat"),
                                     fluidRow(
                                       column(width=5,uiOutput("PC1_sel")),
                                       column(width=2,br(),br(),p("VS")),
                                       column(width=5,uiOutput("PC2_sel"))
                                     )
                                    ),
                    conditionalPanel(condition="input.DiagPlot=='pcoaPlot' || input.DiagPlot=='nmdsPlot' || input.DiagPlot=='SERE' || input.DiagPlot=='clustPlot' ",
                                      uiOutput("DistList")
                                    ),
                    conditionalPanel(condition="input.DistClust=='Unifrac' && (input.DiagPlot=='pcoaPlot' || input.DiagPlot=='nmdsPlot'  || input.DiagPlot=='clustPlot')",
                      column(width=12,
                           selectInput("DistClustUnifrac","Which unifrac distance ?",c("Weighted UniFrac"="WU", "Unweighted UniFrac"="UWU", "Variance adjusted weighted UniFrac"="VAWU"))
                      )
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
                  checkboxInput("modifwidthDiag","Set width",FALSE),
                  conditionalPanel(condition="input.modifwidthDiag",
                  sliderInput("widthDiag", "Width",min=100,max=2500,value = 800,step =10)),

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
                    conditionalPanel(condition="input.DiagPlot=='SfactorsVStot' || input.DiagPlot=='pcaPlot' || input.DiagPlot=='pcoaPlot' || input.DiagPlot=='nmdsPlot' ",column(width=6,sliderInput("cexLabelDiag", h6("Points"),min=0,max=5,value = 1,step =0.1)))
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
                         fluidRow(
                           column(width=8,selectInput("WhichExportTable", "Select the table to export",c("Significant"="Significant","Complete"="Complete","Up"="Up","Down"="Down"))),
                           column(width=4,selectInput("sepexpdiff", "Separator:", c("Tab" = "\t", "Comma" = ",", "Semicolon" = ";")))
                         ),
                         uiOutput("ExportTableButton")
                     )
            )
            ) 
    ),
  
  #### Data Viz
  
  tabItem(tabName = "GlobVisu",
          fluidRow(
            column(width=9,
                   tags$head(tags$style(HTML(spincss))),
                   div(id = "plot-container",
                       conditionalPanel(condition="input.PlotVisuSelect=='Boxplot' || input.PlotVisuSelect=='Diversity' || input.PlotVisuSelect=='Rarefaction'",   
                                          tags$img(src = "gears.gif",id ="loading-spinner")
                                        ),
                       uiOutput("plotVisu"),
                       #conditionalPanel(condition="input.PlotVisuSelect!='Krona'",),
                       #conditionalPanel(condition="input.PlotVisuSelect=='Krona'",tags$script(uiOutput("plotVisu"))),
                       conditionalPanel(condition="input.PlotVisuSelect=='Scatterplot' && !input.AddRegScatter",
                                        p(actionButton("scatterD3-reset-zoom", HTML ("<span class='glyphicon glyphicon-search' aria-hidden='true'></span> Reset Zoom")),Align="right")
                       )
                   ),
                   
                   
                   ### Regression and correlation outputs for the scatter plot
                   conditionalPanel(condition="input.PlotVisuSelect=='Scatterplot' && input.AddRegScatter", 
                                    fluidRow(
                                      column(width=6,
                                      br(),
                                      box(title = "Regression coefficients",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                                          DT::dataTableOutput("lmRegScatter")
                                      )
                                    ),
                                    column(width=6,br(),htmlOutput("lmEquation"))
                                    )
                   ),
                   conditionalPanel(condition="input.PlotVisuSelect=='Scatterplot'",
                                    useShinyjs(),
                                    br(),
                                    box(title = "Correlation table",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                                        DT::dataTableOutput("CorTable")
                                    )
                   ),
                   
                   ## Values of the diversities
                   conditionalPanel(condition="input.PlotVisuSelect=='Diversity'",
                                    br(),
                      box(title = "Diversity values",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                          DT::dataTableOutput("Diversitytable"),
                          fluidRow(
                            column(width=3,downloadButton('ExportDiversitytable', 'Export table')),
                            column(width=3,selectInput("sepdiversity", "Separator:", c("Tab" = "\t", "Comma" = ",", "Semicolon" = ";")))
                          ),
                          tags$style(type='text/css', "#ExportDiversitytable { margin-top: 37px;}")
                      )
                   ) 
                                    
                   ),
            
            column(width=3,
              box(title = "Select your plot",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = FALSE,collapsed= FALSE,
                  selectizeInput("PlotVisuSelect","",c("Barplot"="Barplot","Heatmap"="Heatmap","Boxplot"="Boxplot","Tree"="Tree","Scatterplot"="Scatterplot","Diversity"="Diversity","Rarefaction"="Rarefaction","Krona"="Krona"),selected = "Barplot")
              ),
            

              ########################################################################
              ###
              ###               Options Visualization
              ###
              ########################################################################
              box(title = "Options",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= FALSE,
                  conditionalPanel(condition="input.PlotVisuSelect",
                                   radioButtons("NormOrRaw",label = h5(strong("Type of counts")), c("Normalized" = "norm", "Raw" = "raw"),inline=TRUE)
                  ),
                  conditionalPanel(condition="input.PlotVisuSelect!='Rarefaction' && input.PlotVisuSelect!='Scatterplot'",
                                   uiOutput("VarIntVisu"),
                                   h5(strong("Select the modalities")),
                                   uiOutput("ModVisu")
                  ),
                  # conditionalPanel(condition="input.PlotVisuSelect=='Tree'",
                  #                  uiOutput("VarIntVisuTree")),
                  conditionalPanel(condition="input.PlotVisuSelect=='Scatterplot'",
                                   uiOutput("VarIntVisuScatter"),
                                   radioButtons("TransDataScatter","Data transformation",c("Log2 +1" = "log2","None" = "none"),inline=TRUE),
                                   hr(),
                                   radioButtons("CorMeth","Correlation method",c("Pearson" = "pearson", "Spearman" = "spearman"),inline=TRUE),
                                   checkboxInput("AddRegScatter","Add regression line",FALSE)
                  ),                 
                  conditionalPanel(condition="input.PlotVisuSelect!='Rarefaction' && input.PlotVisuSelect!='Diversity' && input.PlotVisuSelect!='Scatterplot' && input.PlotVisuSelect!='Krona'",
                                   radioButtons("SelectSpecifTaxo","Select the features",c("Most abundant"="Most","All"="All", "Differential features" = "Diff", "Non differential features" = "NoDiff"))
                  ),
                  conditionalPanel(condition="input.PlotVisuSelect!='Rarefaction' && input.PlotVisuSelect!='Diversity' && input.PlotVisuSelect!='Scatterplot' && (input.SelectSpecifTaxo=='Diff' || input.SelectSpecifTaxo=='NoDiff') && input.PlotVisuSelect!='Krona' ",
                                   selectizeInput("ContrastList_table_Visu","",choices = "", multiple = TRUE),
                                   radioButtons("UnionInterContrasts","Union or intersection ?",c("Union"="Union","Intersection"="Inter"),inline = TRUE)
                  ),
                  conditionalPanel(condition="input.PlotVisuSelect!='Rarefaction' && input.PlotVisuSelect!='Diversity' && input.PlotVisuSelect!='Scatterplot' && input.PlotVisuSelect!='Krona'",
                                   uiOutput("TaxoToPlotVisu")
                  ),

                ##################
                ## BARPLOT
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Barplot'",
                                 hr(),
                                 selectizeInput(inputId = "CountsOrProp",label = h6(strong("Type of data")),choices = c("Proportions" = "prop", "Counts" = "counts"),selected = "prop")
                ),
                
                
                ##################
                ## HEATMAP
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Heatmap'",
                                 selectizeInput(inputId = "scaleHeatmap",label = h5(strong("Scale:")),choices = c("None" = "none", "Rows" = "row", "Column" = "col"),selected = "none"),
                                 radioButtons("SortHeatRow","Row Clustering:", c("Yes" ="Yes","No" = "No"),inline=TRUE),
                                 radioButtons("SortHeatColumn","Column Clustering:", c("Yes" ="Yes","No" = "No"), inline=TRUE)
                                 
                ),
                
                ##################
                ## BOXPLOT
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Boxplot'",
                                 hr(),
                                 selectizeInput("typeDataBox",h6(strong("Type of data")),c("Log2"="Log2","Relative"="Relative"))
                ),
                
                ##################
                ## DIVERSITY
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Diversity'",
                                 selectizeInput("WhichDiv",h6(strong("Diversity")),c('Alpha','Beta','Gamma','Shannon','Simpson','Inv.Simpson'),selected  = c('Alpha','Shannon','Simpson','Inv.Simpson'),multiple=TRUE)
                ),
                conditionalPanel(condition="input.PlotVisuSelect=='Diversity'",
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
                checkboxInput("modifwidthVisu","Set width",value=FALSE),
                conditionalPanel(condition="input.modifwidthVisu",
                sliderInput("widthVisu", h6(strong("Width")),min=100,max=4000,value = 800)),
                ##################
                ## BARPLOT
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Barplot'",
                                 sliderInput("rotateXLabel", h6(strong("Rotate X labels (Only vertical orientation)")),min=-90, max=90,value = 0, step = 5)  
                ),
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
                                 radioButtons("DivScale","Scales",c("Fixed"="fixed","Free"="free"),selected = "free",inline=TRUE),
                                radioButtons("DivAddError","Add Error bars",c("Add"="Add","Remove"="Remove"),selected = "Add",inline=TRUE)
                ),
                ##################
                ## HEATMAP
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Heatmap'",
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
                                   column(width=6,sliderInput("rightMargin", h6("Right"),min=0,max=30,value = 6,step = 1)),
                                   column(width=6,sliderInput("lowerMargin", h6("Lower"),min=0,max=30,value = 6,step = 1))
                                 )
                ),
                
                ##################
                ## Scatterplot
                ##################
                conditionalPanel(condition="input.PlotVisuSelect=='Scatterplot'",
                                 fluidRow(
                                   column(width=12,sliderInput("SizeLabelScatter", h6("Label size"),min=0,max=50,value = 10,step = 1))
                                  )
                ),
                
                ##################
                ## ALL
                ##################
                conditionalPanel(condition="input.PlotVisuSelect!='Rarefaction' && input.PlotVisuSelect!='Scatterplot' && input.PlotVisuSelect!='Krona' && input.PlotVisuSelect!='Phylogeny'",
                                 radioButtons(inputId = "SensPlotVisu",label = h6(strong("Orientation")),choices = c("Vertical" = "Vertical", "Horizontal" = "Horizontal"),selected = "Vertical",inline = TRUE)
                )
              ),
              conditionalPanel(condition="input.PlotVisuSelect!='Krona' && input.PlotVisuSelect!='Phylogeny'",
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
          )
  ),

  tabItem(tabName = "CompPlot",
          fluidRow(
            column(width=9,
                   uiOutput("plotVisuComp"),
                   conditionalPanel(condition="input.PlotVisuSelectComp=='Venn'",
                                    DT::dataTableOutput("DataVenn")
                   )
            ),
            column(width=3,
                   box(title = "Select your plot",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = FALSE,collapsed= FALSE,
                        selectizeInput("PlotVisuSelectComp","",c("Venn diagram"="Venn","Heatmap"="Heatmap_comp"),selected = "Heatmap_comp")
                   ),
                   box(title = "Options",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= FALSE,
                          selectizeInput("ContrastList_table_FC",h6(strong("Contrasts (Min = 2)")),choices = "", multiple = TRUE),
                          conditionalPanel(condition="input.PlotVisuSelectComp=='Heatmap_comp'",
                                          radioButtons("SelectSpecifTaxoComp","Select the features",c("Most abundant"="Most","All"="All", "Differential features" = "Diff", "Non differential features" = "NoDiff"))
                                          ),
                          conditionalPanel(condition="input.PlotVisuSelectComp=='Heatmap_comp' && (input.SelectSpecifTaxoComp=='Diff' || input.SelectSpecifTaxoComp=='NoDiff')",
                                          selectizeInput("ContrastList_table_VisuComp","",choices = "", multiple = TRUE),
                                          radioButtons("UnionInterContrastsComp","Union or intersection ?",c("Union"="Union","Intersection"="Inter"),inline = TRUE)
                                          ),
                          conditionalPanel(condition="input.PlotVisuSelectComp=='Heatmap_comp'",
                                          uiOutput("TaxoToPlotVisuComp"),
                                          selectizeInput(inputId = "scaleHeatmapComp",label = h6(strong("Scale:")),choices = c("None" = "none", "Rows" = "row", "Column" = "col"),selected = "none"),
                                          selectInput("SortHeatComp","Sort by:", c("Selection" ="Selection","Values" = "Values","Names"="Names","Auto"="Auto"))
                                          )
                    ),
                   
                   box(title = "Appearance",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                       sliderInput("heightVisuComp", h6(strong("Height")),min=100,max=4000,value = 800),
                       checkboxInput("modifwidthComp","Set width",FALSE),
                       conditionalPanel(condition="input.modifwidthComp",
                                        sliderInput("widthComp", "Width",min=100,max=2500,value = 800,step =10)),
                       ##################
                       ## HEATMAP
                       ##################

                       conditionalPanel(condition="input.PlotVisuSelectComp=='Heatmap_comp'",
                                        radioButtons(inputId = "SensPlotVisuComp",label = h6(strong("Orientation")),choices = c("Vertical" = "Vertical", "Horizontal" = "Horizontal"),selected = "Vertical",inline = TRUE),
                                        fluidRow(
                                          column(width=12,h6(strong("Labels options"))),
                                          column(width=6,sliderInput("LabelSizeHeatmapComp", h6("Size"),min=0.1,max=2,value = 0.7,step = 0.1)),
                                          column(width=6,sliderInput("LabelOrientHeatmapComp", h6("Orientation"),min=0,max=90,value = 0,step = 5)),
                                          column(width=6,sliderInput("LabelColOffsetHeatmapComp", h6("Column offset"),min=0,max=4,value = 0,step = 0.5)),
                                          column(width=6,sliderInput("LabelRowOffsetHeatmapComp", h6("Row offset"),min=0,max=4,value = 0,step = 0.5)),
                                          column(width=12,h6(strong("Margins options"))),
                                          column(width=6,sliderInput("rightMarginComp", h6("Right"),min=0,max=20,value = 6,step = 1)),
                                          column(width=6,sliderInput("lowerMarginComp", h6("Lower"),min=0,max=20,value = 6,step = 1))
                                        )
                       )
                       
                    ),
                   conditionalPanel(condition="input.PlotVisuSelectComp!='Venn'",
                     box(title = "Export",  width = NULL, status = "primary", solidHeader = TRUE,collapsible = TRUE,collapsed= TRUE,
                         ##################
                         ## BARPLOT
                         ##################
                         selectInput("Exp_format_VisuComp",h5(strong("Export format")),c("png"="png","pdf"="pdf","eps"="eps","svg"="svg"), multiple = FALSE),
                         fluidRow(
                           column(width=6,numericInput("heightVisuExportComp", "Height (in px)",min=100,max=NA,value = 500,step =1)),
                           column(width=6,numericInput("widthVisuExportComp", "Width (in px)",min=100,max=NA,value = 500,step =1))
                         ),
                         downloadButton("exportVisuComp", "Export")
                     )
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

          ))
 #includeHTML("file:///home/aghozlan/workspace/SHAMAN_App/www/text.krona.html")
  )
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
}



