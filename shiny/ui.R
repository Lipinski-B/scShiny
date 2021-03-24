# ui.R
shinyUI(dashboardPage(
  ######################### Titre ##########################
  dashboardHeader(title = "Flinovo"),
  
  ########################## Menu ##########################
  dashboardSidebar(
    sidebarMenu(
      menuItem(text = "Lecture", startExpanded = T,
        menuSubItem("Chargement des données",       tabName = "Chargement", icon = icon("readme")),
        menuSubItem("Métadonnées",                  tabName = "metadata", icon = icon("poll"))),
      
      menuItem("Analyses", startExpanded = F,
        menuSubItem("Réduction de dimension",       tabName = "readData", icon = icon("poll")),
        menuSubItem("3D Visualisation",             tabName = "3D_RD", icon = icon("indent-right", lib = "glyphicon")),
        menuSubItem("Analyses PCA",                 tabName = "PCA_facteurs", icon = icon("lock", lib = "glyphicon")),
        menuSubItem("Expression et enrichissement", tabName = "hallmark", icon = icon("calendar")),
        menuSubItem("Trajectoires évolutives",      tabName = "monocle", icon = icon("sort-by-attributes-alt", lib = "glyphicon")),
        menuSubItem("Analyses mitochondriales",     tabName = "mitochondrie", icon = icon("leaf", lib = "glyphicon")))
    )
  ),
  
  ########################## Body ##########################
  dashboardBody(
    ## -- Thème -- ##
    shinyDashboardThemes(theme = "blue_gradient"),
    tags$head(tags$style(HTML('
      .main-header .logo {
        font-family: "system-ui", Times, "Times New Roman", serif;
        font-weight: italic;
        font-size: 36px;
      }
    '))),
    
    ## -- Page -- ##
    tabItems(
      # Chargement
      tabItem(tabName = "Chargement",
              h1("Lecture des données"),
              HTML("Additional text to disply only when menuItem tab One is selected \n"),
              fileInput("dataFile",label = NULL,buttonLabel = "Browse...", placeholder = "No file selected"),
              selectInput("patient", "Sélectionnez le patient : ", choices = c("hFL_180008B","hFL_130337", "FL12C1888", "FL140304"), selected = F),
              div(actionButton(inputId = "actBtnPatient", label = "Apply",icon = icon("play") ), align = "center")
      ),

      # Métadonnées
      tabItem(tabName = "metadata",
              h2("Visualisation des métadonnées"),
              imageOutput('dataTable')
      ),
      
      # Read data
      tabItem(tabName = "readData",

              a(
                href = "https://www.dublinbus.ie/RTPI/Sources-of-Real-Time-Information/",
                "Bus stop numbers can be found here.",
                target = "_blank"
              ),
              br(),
              h3("Custom URL"),
              p("A custom URL can be used to pre select choices when loading the app. Use the button below to create a URL for the choices currently selected."), 
              
              
               
              # Paramètre
              sidebarLayout(
                sidebarPanel(
                  h4("Parameters : "),
                  br(),
                    # Group by choice
                  materialSwitch(inputId = "SwitchGroup", label = "Group by", value = F, status = "primary", inline = T),
                  conditionalPanel(
                      condition = "input.SwitchGroup == true ",
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
                  navbarPage("Reduction dimention",
                             tabPanel("PCA",
                                      plotOutput("PCA", width = "100%",  height = "650px"),
                                      uiOutput("feature_pca")
                             ),
                             tabPanel("UMAP",
                                      plotOutput("UMAP", width = "100%",  height = "650px"),
                                      uiOutput("feature_umap")
                             ),
                             tabPanel("TSNE",
                                      plotOutput("TSNE", width = "100%",  height = "650px"),
                                      uiOutput("feature_tsne")
                             ),
                             tabPanel("Other",
                                      plotlyOutput("feature_other", width = "100%",  height = "650px")
                             )
                  ),
                  
                  fluidRow(),
                  column(1,align="right",checkboxGroupInput("Features", NULL, choices = list("Features" = "Features"), selected = 0)),
                  
                  fluidRow(),
                  conditionalPanel(
                    condition = "input.Features == 'Features' ",
                    fluidRow(column(6,uiOutput('variables'))),
                  ),
                  #Other
                  #column(11,h3("File preview"),dataTableOutput(outputId = "preview")),
                  #tags$br(),
                  #div(actionButton(inputId = "actBtnVisualisation", label = "Visualisation",icon = icon("play") ), align = "center")
                  
               )),
              
      ),
      
      # 3D Reduction dimention 
      tabItem(tabName = "3D_RD",
              radioButtons("metadata", "Metadata:",inline=T,c("seurat_clusters", "SingleR.calls", "Greffe", "clonotype_id","chain", "v_gene", "d_gene", "j_gene","c_gene", "cdr3", "Phase", "HTO_maxID")),
              
              navbarPage("Reduction dimention",
                         tabPanel("PCA",
                                  plotlyOutput("D_PCA", width = "100%",  height = "800px"),
                                  uiOutput("Dfeature_pca")
                                ),
                         tabPanel("UMAP",
                                  plotlyOutput("D_UMAP", width = "100%",  height = "800px"),
                                  uiOutput("Dfeature_umap")
                                  ),
                         tabPanel("TSNE",
                                  plotlyOutput("D_TSNE", width = "100%",  height = "800px"),
                                  uiOutput("Dfeature_tsne")
                                  )
              ),
              
              fluidRow(),
              column(1,align="right",checkboxGroupInput("DFeatures", NULL, choices = list("Features" = "Features"), selected = 0)),
              
              fluidRow(),
              conditionalPanel(
                condition = "input.DFeatures == 'Features'",
                fluidRow(column(6,uiOutput('Dvariables'))),
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
                                  radioGroupButtons(inputId = "color_trajectory", label = "Color by :", choices = meta_variable, justified = TRUE, checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                                  fluidRow(),
                                  plotOutput("cell_trajectory", width = "100%",  height = "600px")
                         ),
                         tabPanel("Ordering", 
                                  radioGroupButtons(inputId = "color_order", label = "Color by :", choices = meta_variable, justified = TRUE, checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                                  fluidRow(),
                                  uiOutput("DDvariables"),
                                  uiOutput("DDorder_trajectory")
                         )
                         
              ),
      ),
      
      # Hallmark
      tabItem(tabName = "hallmark",
              navbarPage("Analyses : ",
                         tabPanel("Expression",
                                  tabsetPanel(
                                  tabPanel("Heatmap", plotOutput("Heatmap", width = "100%",  height = "1400px")),
                                  tabPanel("Info ", verbatimTextOutput("Heatmap_feature")))
                         ),
                         
                         tabPanel("Enrichissement Fonctionnel",
                                  tabsetPanel(
                                  tabPanel("Heatmap",
                                            plotOutput("hallmark_Heatmap", width = "100%",  height = "1000px"),
                                            fluidRow(),
                                            column(1,align="right",checkboxGroupInput("Subsets", NULL, choices = list("Subsets" = "Subsets"), selected = 0)),
                                            
                                            fluidRow(),
                                            conditionalPanel(
                                              condition = "input.Subsets == 'Subsets' ",
                                              wellPanel(h4("Group by :"),
                                                        radioGroupButtons(inputId = "hallmark_order", label = "To order :", choices = meta_variable, justified = TRUE, checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                                                        fluidRow(),
                                                        
                                                        pickerInput(inputId = "numSelector", label = "To subset :", choices = hallmark, multiple = TRUE, options = list(`actions-box` = TRUE)),
                                                        fluidRow(),
                                              ),
          
                                              tags$br(),
                                              div(actionButton(inputId = "actBtnVisualisation", label = "Apply",icon = icon("play") ), align = "center")
                                            )
                                            ),
                                  tabPanel("VlnPlot",
                                            plotOutput("hallmark_VlnPlot", width = "100%",  height = "1000px"),
                                            fluidRow(),
                                            wellPanel(h4("Group by :"),
                                                      pickerInput(inputId = "hallmark_order_vln",label = "Choices : ", choices = hallmark),
                                                      fluidRow(),
                                                      pickerInput(inputId = "metadata_order_vln",label = "Choices : ", choices = meta_variable)
                                            )
                                            ),
                                  tabPanel("Hex Density", 
                                            plotOutput("hallmark_HD", width = "100%",  height = "1000px"),
                                            
                                            wellPanel(h4("Group by :"),
                                                      pickerInput(inputId = "hallmark_order_X",label = "X : ", choices = hallmark),
                                                      fluidRow(),
                                                      pickerInput(inputId = "hallmark_order_Y",label = "Y : ", choices = hallmark),
                                                      fluidRow(),
                                                      pickerInput(inputId = "metadata_order_density", label = "Metadata : ", choices = meta_variable)   
                                            )
                                            ),
                                  tabPanel("Ridge Plot", 
                                            plotOutput("hallmark_RP", width = "100%",  height = "1000px"),
                                            
                                            wellPanel(h4("Group by :"),
                                                      pickerInput(inputId = "hallmark_order_RP", label = "Hallmark : ", choices = hallmark),
                                                      fluidRow(),
                                                      pickerInput(inputId = "metadata_group_RP", label = "Group by : ", choices = meta_variable),
                                                      fluidRow(),
                                                      pickerInput(inputId = "metadata_facet_RP", label = "Facet by : ", choices = meta_variable)  
                                            )))
                                 #tabPanel("PCA", plotOutput("hallmark_PCA", width = "100%",  height = "1000px")),
                         ),
                         
                         tabPanel("Feature Variable",
                                  tabsetPanel(
                                  tabPanel("ScatterPlot", plotOutput("top50", width = "100%",  height = "1000px")),
                                  tabPanel("Top", verbatimTextOutput("Variable_feature")))
                         )
              )
      )  
    ))
))
