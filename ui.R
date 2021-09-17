# ui.R
shinyUI(dashboardPage(
  ######################### Titre ##########################
  dashboardHeader(title = "Flinovo"),
  
  ########################## Menu ##########################
  dashboardSidebar(
    sidebarMenu(
      br(),
      selectInput("patient", h4("Sélectionnez le patient : ", style="color:white"), choices = c("FL08G0293", "FL12C1888", "FL140304", "FL09C1164","FL02G095","FL05G0330","all"), selected = "FL08G0293", selectize = T),
      actionButton(inputId = "actBtnPatient1", label = "Submit",icon = icon("play"), width ='87%'),
      #submitButton("Submit", width = 230),
      br(),
      
      column(1,h4("Information :", style="color:white")),br(),br(), 
      column(1,textOutput("nb_tot_cell")),br(),
      column(1,textOutput("nb_b_cell")),br(),
      column(1,textOutput("nb_other_cell")),br(),br(), 
      column(1,textOutput("nb_Controle_cell")),br(),
      column(1,textOutput("nb_Excipient_cell")),br(),
      column(1,textOutput("nb_RCHOP_cell")),br(),br(),
      verbatimTextOutput("nb_clonotype_cell"),

      menuItem("Analyses", startExpanded = F,
        menuSubItem("Métadonnées",                  tabName = "metadata", icon = icon("poll")),
        menuSubItem("Réduction de dimension",       tabName = "readData", icon = icon("sort-by-attributes-alt", lib = "glyphicon")),
        #menuSubItem("3D Visualisation",            tabName = "3D_RD", icon = icon("indent-right", lib = "glyphicon")),
        #menuSubItem("Analyses PCA",                tabName = "PCA_facteurs", icon = icon("lock", lib = "glyphicon")),
        menuSubItem("Expression et enrichissement", tabName = "hallmark", icon = icon("calendar")),
        #menuSubItem("Trajectoires évolutives",     tabName = "monocle", icon = icon("sort-by-attributes-alt", lib = "glyphicon")),
        menuSubItem("Controle Qualité",             tabName = "mitochondrie", icon = icon("leaf", lib = "glyphicon")),
        menuSubItem("Analyses VDJ",                 tabName = "VDJ", icon = icon("readme")),
        menuSubItem("Subsamples",                   tabName = "Chargement", icon = icon("cog"))),
      menuItem("Méta-Analyses", startExpanded = F,
        menuSubItem("VDJ",                          tabName = "metadata_VDJ", icon = icon("poll")),
        menuSubItem("Cellules",                     tabName = "metadata_cell", icon = icon("poll")),
        menuSubItem("scRepertoir",                  tabName = "metadata_scRepertoir", icon = icon("poll")))
    )
  ),
  
  ########################## Body ##########################
  dashboardBody(
    ## -- Thème -- ##
    dashboardthemes::shinyDashboardThemes(theme = "blue_gradient"),
    tags$head(tags$style(HTML('.main-header .logo {font-family: "system-ui", Times, "Times New Roman", serif;font-weight: italic;font-size: 36px;'))),
    
    ## -- Page -- ##
    tabItems(# Chargement
      tabItem(tabName = "Chargement",
              h1("Subsample :"),

              #HTML("Additional text to disply only when menuItem tab One is selected \n"),
              #fileInput("dataFile",label = NULL, buttonLabel = "Parcourir...", placeholder = "Aucun patient sélectionné", accept = ".RData"),
              #selectInput("patient", "Sélectionnez le patient : ", choices = c("", "FL12C1888", "FL140304", "FL08G0293", "FL09C1164", "FULL"), selected = NULL, selectize = F),
              
              fluidRow(),
              
              fluidRow(
                column(width = 6,
                  box("Metadata filters : ", width = 12, solidHeader = T, collapsible = T,br(),br(),
                    column(6,wellPanel(radioButtons(inputId = "Subgroup", label = NULL, choices = metadata, selected = F))),
                    column(6,wellPanel(uiOutput("Dynamic_Sub_Spe")))
                )),
              
                column(width = 3,
                  box("QC filters : ", width = 12, solidHeader = T, collapsible = T,
                    textInput("maximum", "Nombre maximal de gènes exprimé :", value = NULL, width = NULL, placeholder = NULL),
                    textInput("percent_mt", "Pourcentage seuil d'expression mitochondriale :", value = NULL, width = NULL, placeholder = NULL)
                  ),
                  
                  fluidRow(),
                  box("Gene expression filter : ", width = 12, solidHeader = T, collapsible = T,
                    selectInput('Svariables', "Gène : ", NULL, selectize=TRUE, selected = NULL),
                    textInput("Seuil_variables", "Expression supérieur à :", value = NULL, width = NULL, placeholder = NULL),
                  ),
                  
                  fluidRow(),
                  radioButtons(inputId = "combo", label = "Combinaison to analyse :", choiceNames = c("Excipient/RCHOP", "Pré-greffe/RCHOP", "Excipient/Pré-greffe"), choiceValues = c("EX_RCHOP","PG_RCHOP","PG_EX"),selected = F),
                  
                  fluidRow(),
                  div(actionButton(inputId = "actBtnPatient", label = "Subset",icon = icon("play") ), align = "left", style = "margin-bottom: 10px;", style = "margin-top: -10px;"),
                  div(actionButton(inputId = "resetPatient", label = "Reset",icon = icon("play") ),align = "left", style = "margin-bottom: 10px;", style = "margin-top: -10px;"),
              ))
      ),

      # Métadonnées
      tabItem(tabName = "metadata",
              h1("Présentation du patient et métadonnées : "),
              splitLayout(cellWidths=c("35%","65%"),verticalLayout(fluid = T,
                            tags$head(tags$style(HTML('.shiny-split-layout>div {overflow: hidden;}')),),
                            verbatimTextOutput("nb_clonotype_cell2"),br(),br(),
                            HTML('<center><img src="bcr.png" width="400"></center>')),
                plotlyOutput('dataTable', width = "100%",  height = "800px")
              ),

              downloadButton("PDF_report", "PDF report"),
              downloadButton("HTML_report", "HTML report"),
              downloadButton("Word_report", "Word report")
      ),
      
      # Read data
      tabItem(tabName = "readData",
              
              #a(
              #  href = "https://www.dublinbus.ie/RTPI/Sources-of-Real-Time-Information/",
              #  "Bus stop numbers can be found here.",
              #  target = "_blank"
              #),
              #br(),
              #h3("Custom URL"),
              #p("A custom URL can be used to pre select choices when loading the app. Use the button below to create a URL for the choices currently selected."),
              
              # Paramètre
              sidebarLayout(
                sidebarPanel(h4("Paramètres : "),br(),
                             
                  # Group by choice
                  materialSwitch(inputId = "SwitchGroup", label = "Group by", value = F, status = "primary", inline = T),
                  conditionalPanel(condition = "input.SwitchGroup == true ",
                      column(6,wellPanel(checkboxGroupInput(inputId = "Groupes", label = NULL, choices = metadata, selected = F))),
                      column(6,wellPanel(uiOutput("Dynamic_Group_Spe"))),
                      ),
  
                  # Split by choice
                  fluidRow(),
                  materialSwitch(inputId = "SwitchSplit", label = "Split by", value = F, status = "primary", inline = T),
                  
                  conditionalPanel(
                    condition = "input.SwitchSplit == true ",
                    column(6,wellPanel(checkboxGroupInput(inputId = "Splites", label = NULL, choices = metadata, selected = F),
                              checkboxGroupInput("sclonotype", NULL, choices = list("Clonotype" = "Clonotype"), selected = 1),
                              conditionalPanel(
                                condition = "input.sclonotype == 'Clonotype' ",
                                fluidRow(),
                                sliderInput("NBS_Clonotype", "Clonotype number:", min = 0, max = 5, value = 0),
                              ),
                              fluidRow(),
                    )),
                    column(6,wellPanel(uiOutput("Dynamic_Split_Spe"))),
                  ),
                  fluidRow(),
              ),
              mainPanel(
                  # Visualisation
                  navbarPage("Reduction de dimension",
                             tabPanel("PCA",plotOutput("PCA", width = "100%",  height = "650px"),uiOutput("feature_pca")),
                             tabPanel("UMAP",plotOutput("UMAP", width = "100%",  height = "650px"),uiOutput("feature_umap")),
                             tabPanel("TSNE",plotOutput("TSNE", width = "100%",  height = "650px"),uiOutput("feature_tsne")),
                             tabPanel("Other",plotlyOutput("feature_other", width = "100%",  height = "650px"))
                  ),
                  
                  fluidRow(),
                  column(1,align="right",checkboxGroupInput("Features", NULL, choices = list("Features" = "Features"), selected = 0)),
                  
                  fluidRow(),
                  conditionalPanel(
                    condition = "input.Features == 'Features' ",
                    fluidRow(column(6,selectInput('variables', 'Choices', NULL, multiple=TRUE, selectize=TRUE))),
                  ),
                  #Other
                  #column(11,h3("File preview"),dataTableOutput(outputId = "preview")),
                  #tags$br(),
                  #div(actionButton(inputId = "actBtnVisualisation", label = "Visualisation",icon = icon("play") ), align = "center")
                  
               )),
              
      ),
      
      # 3D Reduction dimention 
      tabItem(tabName = "3D_RD",
              radioButtons("metadata", "Metadata:",inline=T,singlet@tools$meta_variable),
              
              navbarPage("Reduction dimention",
                         tabPanel("PCA",plotlyOutput("D_PCA", width = "100%",  height = "800px"),uiOutput("Dfeature_pca")),
                         tabPanel("UMAP",plotlyOutput("D_UMAP", width = "100%",  height = "800px"),uiOutput("Dfeature_umap")),
                         tabPanel("TSNE",plotlyOutput("D_TSNE", width = "100%",  height = "800px"),uiOutput("Dfeature_tsne"))
              ),
              
              fluidRow(),
              column(1,align="right",checkboxGroupInput("DFeatures", NULL, choices = list("Features" = "Features"), selected = 0)),
              
              fluidRow(),
              conditionalPanel(
                condition = "input.DFeatures == 'Features'",
                fluidRow(column(6,selectInput('Dvariables', 'Choices', NULL, selectize=TRUE, selected = NULL))),
              ),
      ),
      
      # 10 PCA dimensions
      tabItem(tabName = "PCA_facteurs",
              navbarPage("PCA paramètres : ",
                         tabPanel("Facteurs", verbatimTextOutput("PCA10")),
                         tabPanel("ElbowPlot", plotOutput("ElbowPlot", width = "100%",  height = "1000px")),
                         tabPanel("DimHeatmap", plotOutput("DimHeatmap", width = "100%",  height = "1000px")),
                         tabPanel("VizDimLoadings", plotOutput("VizDimLoadings", width = "100%",  height = "1000px")),
                         tabPanel("JackStrawPlot", plotOutput("JackStrawPlot", width = "100%",  height = "1000px"))
              ),
      ),
      
      # Mitochondrie
      tabItem(tabName = "mitochondrie",
              navbarPage("MT paramètres : ",
                         tabPanel("VlnPlot", plotOutput("MT_VlnPlot", width = "100%",  height = "1000px")),
                         tabPanel("FeatureScatter", plotOutput("MT_FeatureScatter", width = "100%",  height = "1000px")),
                         tabPanel("FeatureScatter2", plotOutput("MT_FeatureScatter2", width = "100%",  height = "1000px"))
              ),
      ),
      
      # Monocle
      tabItem(tabName = "monocle",
              navbarPage("Analyses: ",
                         tabPanel("Trajectory", 
                                  radioGroupButtons(inputId = "color_trajectory", label = "Color by :", choices = singlet@tools$meta_variable, justified = TRUE, checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                                  fluidRow(),
                                  plotOutput("cell_trajectory", width = "100%",  height = "600px")
                         ),
                         tabPanel("Ordering", 
                                  radioGroupButtons(inputId = "color_order", label = "Color by :", choices = singlet@tools$meta_variable, justified = TRUE, checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                                  fluidRow(),
                                  selectInput('DDvariables', 'Choices', NULL, selectize=TRUE, multiple=TRUE, selected = NULL),
                                  uiOutput("DDorder_trajectory")
                         )
                         
              ),
      ),
      
      # Hallmark
      tabItem(tabName = "hallmark",
              navbarPage("Analyses : ",
                         tabPanel("Expression Différentielle",
                                  tabsetPanel(
                                    tabPanel("Heatmap",
                                             splitLayout(cellWidths=c("50%","50%"),
                                                         h3("RCHOP vs Excipient :"),
                                                         h3("Pré-greffe vs Excipient :")),
                                             splitLayout(cellWidths=c("50%","50%"),
                                                         plotOutput("DE_Heatmap_RE", width = "95%",  height = "1200px"),
                                                         plotOutput("DE_Heatmap_PE", width = "95%",  height = "1200px"))),
                                    tabPanel("RidgePlot",
                                             splitLayout(cellWidths=c("50%","50%"),
                                                         h3("RCHOP vs Excipient :"),
                                                         h3("Pré-greffe vs Excipient :")),
                                             splitLayout(cellWidths=c("50%","50%"),
                                                         plotOutput("DE_RidgePlot_RE", width = "95%",  height = "1200px"),
                                                         plotOutput("DE_RidgePlot_PE", width = "95%",  height = "1200px"))),
                                    tabPanel("VlnPlot",
                                             splitLayout(cellWidths=c("50%","50%"),
                                                         h3("RCHOP vs Excipient :"),
                                                         h3("Pré-greffe vs Excipient :")),
                                             splitLayout(cellWidths=c("50%","50%"),
                                                         plotOutput("DE_VlnPlot_RE", width = "95%",  height = "1200px"),
                                                         plotOutput("DE_VlnPlot_PE", width = "95%",  height = "1200px"))),
                                    tabPanel("DotPlot",
                                             splitLayout(cellWidths=c("50%","50%"),
                                                         h3("RCHOP vs Excipient :"),
                                                         h3("Pré-greffe vs Excipient :")),
                                             splitLayout(cellWidths=c("50%","50%"),
                                                         plotOutput("DE_DotPlot_RE", width = "95%",  height = "600px"),
                                                         plotOutput("DE_DotPlot_PE", width = "95%",  height = "600px"))),
                                    tabPanel("Linear",
                                             splitLayout(cellWidths=c("50%","50%"),
                                                         h3("RCHOP vs Excipient :"),
                                                         h3("Pré-greffe vs Excipient :")),
                                             splitLayout(cellWidths=c("50%","50%"),
                                                         plotOutput("Linear_RE", width = "95%",  height = "600px"),
                                                         plotOutput("Linear_PE", width = "95%",  height = "600px"))),
                                    tabPanel("Info",
                                             splitLayout(cellWidths=c("50%","50%"),
                                                         h3("RCHOP vs Excipient :"),
                                                         h3("Pré-greffe vs Excipient :")),
                                             splitLayout(cellWidths=c("50%","50%"),
                                                         verbatimTextOutput("DE_info_RE"),
                                                         verbatimTextOutput("DE_info_PE"))),

                                    tabPanel("Seurat Standard", 
                                             plotOutput("Heatmap", width = "100%",  height = "1200px"),
                                             verbatimTextOutput("Heatmap_feature"))
                                    )
                                  
                         ),
                         tabPanel("Enrichissement Fonctionnel",
                                  tabsetPanel(
                                    tabPanel("KEGG", plotOutput("KEGG", width = "100%",  height = "1200px")),
                                    tabPanel("GO : Biological Process", plotOutput("GO_Biological", width = "100%",  height = "1200px")),
                                    tabPanel("GO : Cellular Component", plotOutput("GO_Cellular", width = "100%",  height = "1200px")),
                                    tabPanel("GO : Molecular Function", plotOutput("GO_Molecular", width = "100%",  height = "1200px"))
                                  )
                         ),
                         tabPanel("Hallmark enrichissement",
                                  tabsetPanel(
                                  tabPanel("Heatmap",
                                            plotOutput("hallmark_Heatmap", width = "100%",  height = "1200px"),
                                            fluidRow(),
                                            column(1,align="right",checkboxGroupInput("Subsets", NULL, choices = list("Subsets" = "Subsets"), selected = 0)),
                                            
                                            fluidRow(),
                                            conditionalPanel(
                                              condition = "input.Subsets == 'Subsets' ",
                                              wellPanel(h4("Group by :"),
                                                        radioGroupButtons(inputId = "hallmark_order", label = "To order :", choices = singlet@tools$meta_variable, justified = TRUE, checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                                                        fluidRow(),
                                                        
                                                        pickerInput(inputId = "numSelector", label = "To subset :", choices = singlet@tools$hallmarks, multiple = TRUE, options = list(`actions-box` = TRUE)),
                                                        fluidRow(),
                                              ),
          
                                              tags$br(),
                                              div(actionButton(inputId = "actBtnVisualisation", label = "Apply",icon = icon("play") ), align = "center")
                                            )
                                            ),
                                  tabPanel("VlnPlot",
                                            plotOutput("hallmark_VlnPlot", width = "100%",  height = "800px"),
                                            fluidRow(),
                                            wellPanel(h4("Group by :"),
                                                      pickerInput(inputId = "hallmark_order_vln",label = "Hallmarks : ", choices = singlet@tools$hallmarks),
                                                      fluidRow(),
                                                      pickerInput(inputId = "metadata_order_vln",label = "Métadonnées : ", choices = singlet@tools$meta_variable)
                                            )
                                            ),
                                  tabPanel("Hex Density", 
                                            plotOutput("hallmark_HD", width = "100%",  height = "1000px"),
                                            wellPanel(h4("Group by :"),
                                                      pickerInput(inputId = "hallmark_order_X",label = "X : ", choices = singlet@tools$hallmarks),
                                                      fluidRow(),
                                                      pickerInput(inputId = "hallmark_order_Y",label = "Y : ", choices = singlet@tools$hallmarks),
                                                      fluidRow(),
                                                      pickerInput(inputId = "metadata_order_density", label = "Metadata : ", choices = singlet@tools$meta_variable)   
                                            )
                                            ),
                                  tabPanel("Ridge Plot", 
                                            plotOutput("hallmark_RP", width = "100%",  height = "1000px"),
                                            wellPanel(h4("Group by :"),
                                                      pickerInput(inputId = "hallmark_order_RP", label = "Hallmark : ", choices = singlet@tools$hallmarks),
                                                      fluidRow(),
                                                      pickerInput(inputId = "metadata_group_RP", label = "Group by : ", choices = singlet@tools$meta_variable),
                                                      fluidRow(),
                                                      pickerInput(inputId = "metadata_facet_RP", label = "Facet by : ", choices = singlet@tools$meta_variable)  
                                            )))
                                 #tabPanel("PCA", plotOutput("hallmark_PCA", width = "100%",  height = "1000px")),
                         ),
                         
                         tabPanel("Feature Variable",
                                  tabsetPanel(
                                  tabPanel("ScatterPlot", plotOutput("top50", width = "100%",  height = "1000px")),
                                  tabPanel("Top", verbatimTextOutput("Variable_feature")))
                         )
              )
      ),
      
      # VDJ
      tabItem(tabName = "VDJ",
               box("Clonotype", width = 12, plotlyOutput("VDJ_Clonotype", width = "100%",  height = "400px")),
               box("Heavy", width = 7, plotlyOutput("VDJ_Heavy", width = "100%",  height = "400px"),
                        column(1,align="right",checkboxGroupInput("heavy_details", NULL, choices = list("Details" = "Details"), selected = 0)),
                        conditionalPanel(
                          condition = "input.heavy_details == 'Details' ", br(),
                          plotlyOutput("VDJ_DHeavy", width = "100%",  height = "400px"))
                        ),
               box("Lights", width = 5, plotlyOutput("VDJ_Light", width = "100%",  height = "400px"),
                        column(1,align="right",checkboxGroupInput("light_details", NULL, choices = list("Details" = "Details"), selected = 0)),
                        conditionalPanel(
                          condition = "input.light_details == 'Details' ", br(),
                          plotlyOutput("VDJ_DLight", width = "100%",  height = "400px"))
                        ),
               fluidRow(
                 box("V", width = 4, plotlyOutput("V", width = "100%",  height = "400px")),
                 box("D", width = 4, plotlyOutput("D", width = "100%",  height = "400px")),
                 box("J", width = 4, plotlyOutput("J", width = "100%",  height = "400px")))
      ),
      
      ############################################################################################
      # metadata VDJ
      tabItem(tabName = "metadata_VDJ",
              fluidRow(
                box("", width = 6, plotlyOutput("FL08_VDJ", width = "100%",  height = "600px")),
                box("", width = 6, plotlyOutput("FL09_VDJ", width = "100%",  height = "600px"))),
              fluidRow(
                box("", width = 6, plotlyOutput("FL12_VDJ", width = "100%",  height = "600px")),
                box("", width = 6, plotlyOutput("FL14_VDJ", width = "100%",  height = "600px"))),
              fluidRow(
                box("", width = 6, plotlyOutput("FL02_VDJ", width = "100%",  height = "600px")),
                box("", width = 6, plotlyOutput("FL05_VDJ", width = "100%",  height = "600px"))),

      ),
      # metadata cell
      tabItem(tabName = "metadata_cell",
              fluidRow(
                box("", width = 6, plotlyOutput("FL08_Meta", width = "100%",  height = "600px")),
                box("", width = 6, plotlyOutput("FL09_Meta", width = "100%",  height = "600px"))),
              fluidRow(
                box("", width = 6, plotlyOutput("FL12_Meta", width = "100%",  height = "600px")),
                box("", width = 6, plotlyOutput("FL14_Meta", width = "100%",  height = "600px"))),
              fluidRow(
                box("", width = 6, plotlyOutput("FL02_Meta", width = "100%",  height = "600px")),
                box("", width = 6, plotlyOutput("FL05_Meta", width = "100%",  height = "600px"))),
      ),
      # metadata scRepertoir
      tabItem(tabName = "metadata_scRepertoir",
              navbarPage("Figure : ",
                         tabPanel("Homeostasis",HTML('<left><img src="scRepertoir/clonalHomeostasis.png" width="1000"></left>')),
                         tabPanel("Proportion",HTML('<left><img src="scRepertoir/clonalProportion.png" width="900"></left>')),
                         tabPanel("Contig length",HTML('<left><img src="scRepertoir/lengthContig.png" width="900"></left>')),
                         tabPanel("Contig quant",HTML('<left><img src="scRepertoir/quantContig_output.png" width="700"></left>')))
      )
    ))
))
