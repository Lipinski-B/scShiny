mod_Patient_QC_ui <- function(id) {
  ns <- NS(id)
  
  # QC
  tabItem(tabName = "Patient_QC",
          h5("Group by : "),
          pickerInput(ns("QC_Groupes"), inline = FALSE, choices = c("Sample","seurat_clusters", "Phénotype", "Phénotype.fine","Phase","Condition", "Clonotype","Chaine","V","D","J","C","CDR3","BCL2_L23L", "BCL2_K22K", "CD79B_Y196H", "EZH2_A682G", "EZH2_A692V","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S"), multiple = FALSE, options = list(`actions-box` = TRUE)),fluidRow(),
          
          
          h4("Controle Quality Summary :"),
          splitLayout(cellWidths = c("50%", "50%"),
            plotOutput(ns("QC_summary"), width = "100%",  height = "500px"),
            plotOutput(ns("QC_summary2"), width = "100%",  height = "500px")
          ), fluidRow(),
          
          h4("Count RNA vs Percent's metrics :"),
          splitLayout(cellWidths = c("50%", "50%"),
                      splitLayout(cellWidths = c("50%", "50%"),
                                  plotOutput(ns("nRNA"), width = "100%",  height = "500px"),
                                  plotOutput(ns("nRNA_mt"), width = "100%",  height = "500px")),
                      splitLayout(cellWidths = c("50%", "50%"),
                                  plotOutput(ns("nRNA_rb"), width = "100%",  height = "500px"),
                                  plotOutput(ns("nRNA_ig"), width = "100%",  height = "500px"))
          ), fluidRow(),
          
          h4("Count SCT vs Percent's metrics :"),
          splitLayout(cellWidths = c("50%", "50%"),
                      splitLayout(cellWidths = c("50%", "50%"),
                                  plotOutput(ns("nSCT"), width = "100%",  height = "500px"),
                                  plotOutput(ns("nSCT_mt"), width = "100%",  height = "500px")),
                      splitLayout(cellWidths = c("50%", "50%"),
                                  plotOutput(ns("nSCT_rb"), width = "100%",  height = "500px"),
                                  plotOutput(ns("nSCT_ig"), width = "100%",  height = "500px"))
          ), fluidRow(),
          
          h4("Percent's metrics correlations :"),
          splitLayout(cellWidths = c("66%", "33%"),
                      splitLayout(cellWidths = c("50%", "50%"),
                                  plotOutput(ns("rb_mt"), width = "100%",  height = "500px"),
                                  plotOutput(ns("ig_mt"), width = "100%",  height = "500px")),
                                  plotOutput(ns("ig_rb"), width = "100%",  height = "500px"),
          ),
  )
}

mod_Patient_QC_server <- function(input, output, session, r) {
  ns <- session$ns

  # output$MT_VlnPlot <- renderPlot({r$dataset@tools$mitochondrie_all + Seurat::VlnPlot(r$dataset, features = "nFeature_RNA", group.by = "Condition") +  Seurat::VlnPlot(r$dataset, features = "nCount_RNA", group.by = "Condition") + Seurat::VlnPlot(r$dataset, features = "percent.mt", group.by = "Condition")})
  # output$MT_FeatureScatter <- renderPlot({Seurat::FeatureScatter(r$dataset, feature1 = "nCount_RNA", feature2 = "percent.mt") + Seurat::FeatureScatter(r$dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")})
  # output$MT_FeatureScatter2 <- renderPlot({Seurat::FeatureScatter(r$dataset, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "HTO_classification") + Seurat::FeatureScatter(r$dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "HTO_classification")})

  reactive({input$QC_Groupes})
  
  output$QC_summary <- renderPlot({Seurat::VlnPlot(r$dataset, group.by = input$QC_Groupes, features = c("percent.mt","percent.rb","percent.ig"),ncol = 3,pt.size = 0.1) & theme(plot.title = element_text(size=10)) })
  output$QC_summary2 <- renderPlot({Seurat::VlnPlot(r$dataset, group.by = input$QC_Groupes, features = c("nFeature_RNA","nCount_RNA","nFeature_SCT","nCount_SCT"),ncol = 4,pt.size = 0.1) & theme(plot.title = element_text(size=10))})
  
  output$nRNA <- renderPlot({Seurat::FeatureScatter(r$dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = input$QC_Groupes)})
  output$nRNA_mt <- renderPlot({Seurat::FeatureScatter(r$dataset, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = input$QC_Groupes)})
  output$nRNA_rb <- renderPlot({Seurat::FeatureScatter(r$dataset, feature1 = "nCount_RNA", feature2 = "percent.rb", group.by = input$QC_Groupes)})
  output$nRNA_ig <- renderPlot({Seurat::FeatureScatter(r$dataset, feature1 = "nCount_RNA", feature2 = "percent.ig", group.by = input$QC_Groupes)})
  
  output$nSCT <- renderPlot({Seurat::FeatureScatter(r$dataset, feature1 = "nCount_SCT", feature2 = "nFeature_SCT", group.by = input$QC_Groupes)})
  output$nSCT_mt <- renderPlot({Seurat::FeatureScatter(r$dataset, feature1 = "nCount_SCT", feature2 = "percent.mt", group.by = input$QC_Groupes)})
  output$nSCT_rb <- renderPlot({Seurat::FeatureScatter(r$dataset, feature1 = "nCount_SCT", feature2 = "percent.rb", group.by = input$QC_Groupes)})
  output$nSCT_ig <- renderPlot({Seurat::FeatureScatter(r$dataset, feature1 = "nCount_SCT", feature2 = "percent.ig", group.by = input$QC_Groupes)})
  
  output$rb_mt <- renderPlot({Seurat::FeatureScatter(r$dataset, feature1 = "percent.rb", feature2 = "percent.mt", group.by = input$QC_Groupes)})
  output$ig_mt <- renderPlot({Seurat::FeatureScatter(r$dataset, feature1 = "percent.ig", feature2 = "percent.mt", group.by = input$QC_Groupes)})
  output$ig_rb <- renderPlot({Seurat::FeatureScatter(r$dataset, feature1 = "percent.ig", feature2 = "percent.rb", group.by = input$QC_Groupes)})
  
}


## To be copied in the UI
# mod_Patient_QC_ui("QC_ui_1")

## To be copied in the server
# callModule(mod_Patient_QC_server, "QC_ui_1")