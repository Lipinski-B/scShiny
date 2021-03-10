# ui.R
shinyUI(dashboardPage(
  # Titre de l'application
  dashboardHeader(title = "Seurat object"),
  
  # Menu
  dashboardSidebar(
    sidebarMenu(
      menuItem("Réduction de dimension", tabName = "readData", icon = icon("indent-right", lib = "glyphicon")),
      menuItem("3D Visualisation", tabName = "3D_RD", icon = icon("indent-right", lib = "glyphicon")),
      menuItem("Analyses PCA", tabName = "PCA_facteurs", icon = icon("lock", lib = "glyphicon")),
      menuItem("Expression et enrichissement", tabName = "hallmark", icon = icon("calendar")),
      menuItem("Trajectoires évolutives", tabName = "monocle", icon = icon("calendar")),
      #menuItem("Top50 gènes plus variables", tabName = "top50", icon = icon("sort-by-attributes-alt", lib = "glyphicon")),
      menuItem("Analyses mitochondriales", tabName = "mitochondrie", icon = icon("leaf", lib = "glyphicon")),
      menuItem("Lecture des données", tabName = "visualization", icon = icon("readme"))
    )
  ),
  
  # Body
  dashboardBody(
    tabItems(
      tabItem(tabName = "visualization",
              # Input: Checkbox if file has header
              selectInput("patient", "Sélectionnez le patient : ",choices = c("hFL_180008B","hFL_130337")),
              
              h1("Lecture des données"),
              fileInput("dataFile",label = NULL,buttonLabel = "Browse...", placeholder = "No file selected"),
              dataTableOutput('dataTable')
      ),
      # Read data
      tabItem(tabName = "readData",
              h3("Parameters :"),
              #################################################################################################
              wellPanel(h4("Group by :"),
                        useShinyjs(),
                        uiOutput("Dynamic_Group"),
                        fluidRow(),
              ),
              wellPanel(h4("Specifications :"),
                        uiOutput("Dynamic_Group_Spe"),
                        fluidRow(),
              ),

              
              #################################################################################################
              fluidRow(),
              wellPanel(h4("Split by :"),
                        useShinyjs(),
                        uiOutput("Dynamic_Split"),
                        column(1,align="right",checkboxGroupInput("sclonotype", NULL, choices = list("Clonotype" = "Clonotype"), selected = 1)),
                        fluidRow(),
              ),
              wellPanel(h4("Specification :"),
                        uiOutput("Dynamic_Split_Spe"),
                        conditionalPanel(
                          condition = "input.sclonotype == 'Clonotype' ",
                          fluidRow(),
                          sliderInput("NBS_Clonotype", "Clonotype number:", min = 0, max = 5, value = 0),
                        ),
                        
                        fluidRow(),
              ),
              
              
              #################################################################################################
              #Visualisation
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
              
      ),
      # 3D Reduction dimention 
      tabItem(tabName = "3D_RD",
              radioButtons("metadata", "Metadata:",inline=T,
                           c("seurat_clusters", "SingleR.calls", "Greffe", "clonotype_id","chain", "v_gene", "d_gene", "j_gene","c_gene", "cdr3", "Phase", "HTO_maxID")
                           ),
              
              
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
      
      # mitochondrie
      tabItem(tabName = "mitochondrie",
              navbarPage("MT paramètres : ",
                         tabPanel("VlnPlot", plotOutput("MT_VlnPlot", width = "100%",  height = "1000px")),
                         tabPanel("FeatureScatter", plotOutput("MT_FeatureScatter", width = "100%",  height = "1000px")),
                         tabPanel("FeatureScatter2", plotOutput("MT_FeatureScatter2", width = "100%",  height = "1000px"))
              ),
      ),
      
      # monocle
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
      
      
      # hallmark
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
