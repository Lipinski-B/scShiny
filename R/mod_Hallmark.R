mod_Hallmark_ui <- function(id) {
  ns <- NS(id)
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
  )}

mod_Hallmark_server <- function(input, output, session) {
  ns <- session$ns
}


## To be copied in the UI
# mod_Hallmark_ui("Hallmark_ui_1")

## To be copied in the server
# callModule(mod_Hallmark_server, "Hallmark_ui_1")