mod_All_Reduction_Dimension_ui <- function(id) {
  ns <- NS(id)
  
  tabItem(tabName = "All_Reduction_Dimension",
          splitLayout(cellWidths = c("80%", "20%"),
                      tagList(plotOutput(ns("PCA"), width = "100%",  height = "650px"),uiOutput(ns("feature_pca"))),
                      tagList(
                        column(10,shinyBS::bsCollapse(id = "collapse", open = "Paramètres : ",
                          shinyBS::bsCollapsePanel("Paramètres : ",
                          column(12,h5("Reduction :")),column(12,pickerInput(inputId = ns("Reduction"),label = NULL , choices = c("umap","pca","tsne"), multiple = F, options = list(`actions-box` = TRUE))),fluidRow(),
                          column(12,h5("Group by : ")),column(12,pickerInput(ns("test_Groupes"), inline = FALSE, choices = c("Sample","Phénotype_scType", "Heavy", "Light","Phénotype_TRUST","V_Heavy_TRUST","D_Heavy_TRUST","J_Heavy_TRUST","Heavy_TRUST","CDR3_DNA_Heavy_TRUST","CDR3_AA_Heavy_TRUST","V_Light_TRUST","D_Light_TRUST","J_Light_TRUST","Light_TRUST","CDR3_DNA_Light_TRUST","CDR3_AA_Light_TRUST","State","Phénotype", "Phénotype.fine","Phase","Condition","Greffe","Reponse", "Clonotype","Chaine_Light","V_Light","D_Light","J_Light","C_Light","CDR3_Light","CDR3_nt_Light","Chaine_Heavy", "V_Heavy","D_Heavy","J_Heavy","C_Heavy","CDR3_Heavy","CDR3_nt_Heavy","BCL2_L23L", "BCL2_K22K", "CD79B_Y196H", "EZH2_A682G", "EZH2_A692V","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S"), multiple = TRUE, options = list(`actions-box` = TRUE))),fluidRow(),
                          column(12,h5("Split by : ")),column(12,pickerInput(ns("test_Splites"), inline = FALSE, choices = c("Sample","Phénotype_scType", "Heavy", "Light","Phénotype_TRUST","V_Heavy_TRUST","D_Heavy_TRUST","J_Heavy_TRUST","Heavy_TRUST","CDR3_DNA_Heavy_TRUST","CDR3_AA_Heavy_TRUST","V_Light_TRUST","D_Light_TRUST","J_Light_TRUST","Light_TRUST","CDR3_DNA_Light_TRUST","CDR3_AA_Light_TRUST","State","Phénotype", "Phénotype.fine","Phase","Condition","Greffe","Reponse", "Clonotype","Chaine_Light","V_Light","D_Light","J_Light","C_Light","CDR3_Light","CDR3_nt_Light","Chaine_Heavy", "V_Heavy","D_Heavy","J_Heavy","C_Heavy","CDR3_Heavy","CDR3_nt_Heavy","BCL2_L23L", "BCL2_K22K", "CD79B_Y196H", "EZH2_A682G", "EZH2_A692V","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S"), multiple = TRUE, options = list(`actions-box` = TRUE))),fluidRow(),
                          column(12,sliderTextInput(ns("Point_Size"),"Taille des points :",choices = seq(from = 0,to = 1,by = 0.1),grid = TRUE,selected = 0.3)),fluidRow(),
                          column(12,awesomeRadio(inputId = ns('Assays'), label= 'Assay : ', choices = c('RNA','SCT','integrated'),inline = T)),
                          column(12,awesomeCheckbox(ns("Label_element"),"Show Label",TRUE)),fluidRow(),
                          column(12,selectizeInput(inputId = ns('variables'), label= 'Expression : ', choices = NULL, multiple=TRUE)),br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br()
                          ))))
          ))
  
  }

mod_All_Reduction_Dimension_server <- function(input, output, session, r) {
  ns <- session$ns
  
  reactive({input$variables})
  reactive({input$test_Groupes})
  reactive({input$test_Splites})
  reactive({input$Assays})
  
  observe({updateSelectizeInput(session, 'variables', choices = c( "RCHOP_AUC_Score", "percent.mt", "percent.rb",  "percent.ig", "percent.mt.meta",  "percent.rb.meta",  "percent.ig.meta", "nCount_RNA","nFeature_RNA","nCount_SCT","nFeature_SCT","nCount_HTO","nFeature_HTO", "RCHOP_SCT_Score1", "RCHOP_RNA_Score1",  "RCHOP_i_Score1", r$feature), server = TRUE)})
  
  ## -- Réduction de dimensions -- ## 
  output$PCA <- renderPlot({Seurat::DimPlot(object = r$dataset, group.by = input$test_Groupes, split.by = input$test_Splites, label.size = 5, pt.size = input$Point_Size, reduction = input$Reduction, label = input$Label_element) & Seurat::NoLegend() & xlab(label = paste0(input$Reduction, " / PCA 1 : ", round(Seurat::Stdev(r$dataset[["pca"]])[1],2), " %")) & ylab(label = paste0(input$Reduction, " / PCA 2 : ", round(Seurat::Stdev(r$dataset[["pca"]])[2],2), " %"))})
  output$feature_pca = renderUI({if(length(input$variables)>0){plotOutput(ns("FeaturePlot_PCA"), width = "100%",  height = "500px")}})
  output$FeaturePlot_PCA <- renderPlot({
    Seurat::DefaultAssay(r$dataset) <- input$Assays
    Seurat::FeaturePlot(r$dataset, features = input$variables, reduction = input$Reduction, split.by = input$test_Splites, pt.size = input$Point_Size, combine = T , ncol = 3) & Seurat::NoAxes()
    })
}

## To be copied in the UI
# mod_All_Reduction_Dimension_ui("All_Reduction_Dimension_ui_1")

## To be copied in the server
# callModule(mod_All_Reduction_Dimension_server, "All_Reduction_Dimension_ui_1")