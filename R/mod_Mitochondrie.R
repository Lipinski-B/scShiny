mod_Mitochondrie_ui <- function(id) {
  ns <- NS(id)
  
  # Mitochondrie
  tabItem(tabName = "mitochondrie",
          navbarPage("MT paramètres : ",
                     tabPanel("VlnPlot", plotOutput("MT_VlnPlot", width = "100%",  height = "1000px")),
                     tabPanel("FeatureScatter", plotOutput("MT_FeatureScatter", width = "100%",  height = "1000px")),
                     tabPanel("FeatureScatter2", plotOutput("MT_FeatureScatter2", width = "100%",  height = "1000px"))
          )
  )
}

mod_Mitochondrie_server <- function(input, output, session, r) {
  ns <- session$ns

  data <- reactive(r$test$data)
  ## -- Mitochondrie Figure -- ##
  output$MT_VlnPlot <- renderPlot({data()@tools$mitochondrie_all + Seurat::VlnPlot(data(), features = "nFeature_RNA", group.by = "Condition") +  Seurat::VlnPlot(data(), features = "nCount_RNA", group.by = "Condition") + Seurat::VlnPlot(data(), features = "percent.mt", group.by = "Condition")})
  output$MT_FeatureScatter <- renderPlot({Seurat::FeatureScatter(data(), feature1 = "nCount_RNA", feature2 = "percent.mt") + Seurat::FeatureScatter(data(), feature1 = "nCount_RNA", feature2 = "nFeature_RNA")})
  output$MT_FeatureScatter2 <- renderPlot({Seurat::FeatureScatter(data(), feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "HTO_classification") + Seurat::FeatureScatter(data(), feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "HTO_classification")})
  
}


## To be copied in the UI
# mod_Mitochondrie_ui("Mitochondrie_ui_1")

## To be copied in the server
# callModule(mod_Mitochondrie_server, "Mitochondrie_ui_1")