mod_Read_ui <- function(id) {
  ns <- NS(id)

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
          
  )
  
  }

mod_Read_server <- function(input, output, session) {
  ns <- session$ns
  
  
  output$PCA <- renderPlot({
    Seurat::DimPlot(object = tokeep(), group.by = group(), split.by = split(), label.size = 0.0, pt.size = 2, reduction = 'pca') & theme(title = element_text(size=20),legend.position = "top",legend.title = element_text(size=10),legend.text = element_text(size=10)
    ) & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 6))) & xlab(label = paste0("PCA 1 : ", round(Seurat::Stdev(singlet[["pca"]])[1],2), " %")) & ylab(label = paste0("PCA 2 : ", round(Seurat::Stdev(singlet[["pca"]])[2],2), " %"))
  })
  
  output$UMAP <- renderPlot({Seurat::DimPlot(object = tokeep(), group.by = group(), split.by = split(), label.size = 0.0, pt.size = 2, reduction = 'umap') & theme(legend.position = "top") & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4)))})
  
  output$TSNE <- renderPlot({Seurat::DimPlot(object = tokeep(), group.by = group(), split.by = split(), label.size = 0.0, pt.size = 2, reduction = 'tsne') & theme(legend.position = "top") & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4)))})
  
  
}


## To be copied in the UI
# mod_Read_ui("Read_ui_1")

## To be copied in the server
# callModule(mod_Read_server, "Read_ui_1")