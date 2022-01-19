colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
mod_Patient_Expression_ui <- function(id) {
  ns <- NS(id)
  # Expression
  tabItem(tabName = "Patient_Expression",
          navbarPage("Analyses : ",
                     #################################################################################################
                     tabPanel("Expression Différentielle",
                              tabsetPanel(
                                tabPanel("Heatmap",
                                         splitLayout(cellWidths=c("50%","50%"),h3("RCHOP vs Excipient :"),h3("Pré-greffe vs Excipient :")),
                                         splitLayout(cellWidths=c("50%","50%"),plotOutput(ns("DE_Heatmap_RE"), width = "95%",  height = "1200px"),plotOutput(ns("DE_Heatmap_PE"), width = "95%",  height = "1200px"))),
                                tabPanel("RidgePlot",
                                         splitLayout(cellWidths=c("50%","50%"),h3("RCHOP vs Excipient :"),h3("Pré-greffe vs Excipient :")),
                                         splitLayout(cellWidths=c("50%","50%"),plotOutput(ns("DE_RidgePlot_RE"), width = "95%",  height = "1200px"),plotOutput(ns("DE_RidgePlot_PE"), width = "95%",  height = "1200px"))),
                                tabPanel("VlnPlot",
                                         splitLayout(cellWidths=c("50%","50%"),h3("RCHOP vs Excipient :"),h3("Pré-greffe vs Excipient :")),
                                         splitLayout(cellWidths=c("50%","50%"),plotOutput(ns("DE_VlnPlot_RE"), width = "95%",  height = "1200px"),plotOutput(ns("DE_VlnPlot_PE"), width = "95%",  height = "1200px"))),
                                tabPanel("DotPlot",
                                         splitLayout(cellWidths=c("50%","50%"),h3("RCHOP vs Excipient :"),h3("Pré-greffe vs Excipient :")),
                                         splitLayout(cellWidths=c("50%","50%"),plotOutput(ns("DE_DotPlot_RE"), width = "95%",  height = "600px"),plotOutput(ns("DE_DotPlot_PE"), width = "95%",  height = "600px"))),
                                # tabPanel("Linear",
                                #          splitLayout(cellWidths=c("50%","50%"),h3("RCHOP vs Excipient :"),h3("Pré-greffe vs Excipient :")),
                                #          splitLayout(cellWidths=c("50%","50%"),plotOutput(ns("Linear_RE"), width = "95%",  height = "600px"),plotOutput(ns("Linear_PE"), width = "95%",  height = "600px"))),
                                tabPanel("Info",
                                         splitLayout(cellWidths=c("50%","50%"),h3("RCHOP vs Excipient :"),h3("Pré-greffe vs Excipient :")),
                                         splitLayout(cellWidths=c("50%","50%"),verbatimTextOutput(ns("DE_info_RE")),verbatimTextOutput(ns("DE_info_PE")))),
                                tabPanel("Seurat Standard", 
                                         plotOutput(ns("Heatmap"), width = "100%",  height = "1200px"),
                                         verbatimTextOutput(ns("Heatmap_feature")))
                              )
                     ),
                     #################################################################################################
                     tabPanel("Hallmark enrichissement",
                              tabsetPanel(
                                tabPanel("Heatmap",
                                         plotOutput(ns("hallmark_Heatmap"), width = "95%",  height = "800px"),fluidRow(),
                                         column(1,align="right",checkboxGroupInput("Subsets", NULL, choices = list("Subsets" = "Subsets"), selected = 0)),fluidRow(),
                                         
                                         conditionalPanel(
                                           condition = "input.Subsets == 'Subsets' ",
                                           wellPanel(h4("Group by :"),
                                                     #radioGroupButtons(inputId = ns("hallmark_order"), label = "To order :", choices = singlet@tools$meta_variable, justified = TRUE, checkIcon = list(yes = icon("ok",lib = "glyphicon"))),fluidRow(),                                                     
                                                     pickerInput(inputId = ns("numSelector"), label = "To subset :", choices = singlet@tools$hallmarks, selected=singlet@tools$hallmarks, multiple = TRUE, options = list(`actions-box` = TRUE)),fluidRow(),
                                           ),
                                           
                                           #tags$br(), div(actionButton(inputId = ns("actBtnVisualisation"), label = "Apply",icon = icon("play") ), align = "center")
                                         )
                                ),
                                tabPanel("VlnPlot",
                                         plotOutput(ns("hallmark_VlnPlot"), width = "100%",  height = "800px"),fluidRow(),
                                         wellPanel(h4("Group by :"),
                                                   pickerInput(inputId = ns("hallmark_order_vln"),label = "Hallmarks : ", choices = singlet@tools$hallmarks),fluidRow(),
                                                   pickerInput(inputId = ns("metadata_order_vln"),label = "Métadonnées : ", choices = singlet@tools$meta_variable)
                                         )
                                ),
                                tabPanel("Hex Density", 
                                         plotOutput(ns("hallmark_HD"), width = "100%",  height = "1000px"),
                                         wellPanel(h4("Group by :"),
                                                   pickerInput(inputId = ns("hallmark_order_X"),label = "X : ", choices = singlet@tools$hallmarks),fluidRow(),
                                                   pickerInput(inputId = ns("hallmark_order_Y"),label = "Y : ", choices = singlet@tools$hallmarks),fluidRow(),
                                                   pickerInput(inputId = ns("metadata_order_density"), label = "Metadata : ", choices = singlet@tools$meta_variable)   
                                         )
                                ),
                                tabPanel("Ridge Plot", 
                                         plotOutput(ns("hallmark_RP"), width = "100%",  height = "1000px"),
                                         wellPanel(h4("Group by :"),
                                                   pickerInput(inputId = ns("hallmark_order_RP"), label = "Hallmark : ", choices = singlet@tools$hallmarks),fluidRow(),
                                                   pickerInput(inputId = ns("metadata_group_RP"), label = "Group by : ", choices = singlet@tools$meta_variable),fluidRow(),
                                                   pickerInput(inputId = ns("metadata_facet_RP"), label = "Facet by : ", choices = singlet@tools$meta_variable)  
                                         )),
                              tabPanel("PCA", plotOutput(ns("hallmark_PCA"), width = "100%",  height = "1000px")))
                     ),
                     #################################################################################################
                     tabPanel("Feature Variable",
                              tabsetPanel(
                                tabPanel("ScatterPlot", plotOutput(ns("top50"), width = "100%",  height = "1000px")),
                                tabPanel("Top", verbatimTextOutput(ns("Variable_feature"))))
                     )
          )
  )}

mod_Patient_Expression_server <- function(input, output, session, r) {
  ns <- session$ns

  #################################################################################################
  ## -- Expression différentielle -- ##
  observe(Seurat::Idents(r$dataset)<-"Condition")
  
  # RCHOP / Excipient
  output$DE_Heatmap_RE <- renderPlot({Seurat::DoHeatmap(r$dataset, cells = rownames(r$dataset@meta.data)[which(r$dataset@meta.data$Condition==c("Excipient","RCHOP"))], features = rownames(r$dataset@tools$DE_RE)[1:50], size = 3, assay = 'RNA', slot = "data")})
  output$DE_RidgePlot_RE <- renderPlot({Seurat::RidgePlot(r$dataset, idents = c("Excipient","RCHOP"), features = rownames(r$dataset@tools$DE_RE)[1:12], ncol = 3, assay = 'RNA', slot = "data") & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")})
  output$DE_VlnPlot_RE <- renderPlot({Seurat::VlnPlot(r$dataset, idents = c("Excipient","RCHOP"), features = rownames(r$dataset@tools$DE_RE)[1:12], sort = T, ncol = 3, assay = 'RNA', slot = "data") & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom")})
  output$DE_DotPlot_RE <- renderPlot({Seurat::DotPlot(r$dataset, idents = c("Excipient","RCHOP"), features = rownames(r$dataset@tools$DE_RE)[1:12]) + Seurat::RotatedAxis()})
  output$DE_info_RE <- renderPrint({print(r$dataset@tools$DE_RE)})
  # output$Linear_RE <- renderPlot({
  #   result <- rownames(r$dataset@tools$avg.b.cells_RE[which(r$dataset@tools$avg.b.cells_RE$CTRL/r$dataset@tools$avg.b.cells_RE$STIM < 0.75),])
  #   p1 <- ggplot(r$dataset@tools$avg.b.cells_RE, aes(CTRL, STIM)) + geom_point() + ggtitle("B Cells") 
  #   Seurat::LabelPoints(plot = p1, points = result, repel = TRUE, xnudge = 0.2, ynudge = 0.5)
  # })
  
  # Pré-greffe / Excipient
  output$DE_Heatmap_PE <- renderPlot({Seurat::DoHeatmap(r$dataset, cells = rownames(r$dataset@meta.data)[which(r$dataset@meta.data$Condition==c("Excipient","Pré-greffe"))], features = rownames(r$dataset@tools$DE_PE)[1:50], size = 3, assay = 'RNA', slot = "data")})
  output$DE_RidgePlot_PE <- renderPlot({Seurat::RidgePlot(r$dataset, idents = c("Excipient","Pré-greffe"), features = rownames(r$dataset@tools$DE_PE)[1:12], ncol = 3, assay = 'RNA', slot = "data") & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")})
  output$DE_VlnPlot_PE <- renderPlot({Seurat::VlnPlot(r$dataset, idents = c("Excipient","Pré-greffe"), features = rownames(r$dataset@tools$DE_PE)[1:12], sort = T, ncol = 3, assay = 'RNA', slot = "data") & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom")})
  output$DE_DotPlot_PE <- renderPlot({Seurat::DotPlot(r$dataset, idents = c("Excipient","Pré-greffe"), features = rownames(r$dataset@tools$DE_PE)[1:12]) + Seurat::RotatedAxis()})
  output$DE_info_PE <- renderPrint({print(r$dataset@tools$DE_PE)})
  # output$Linear_PE <- renderPlot({
  #   result <- rownames(r$dataset@tools$avg.b.cells_PE[which(r$dataset@tools$avg.b.cells_PE$CTRL/r$dataset@tools$avg.b.cells_PE$STIM > 2),])
  #   p1 <- ggplot(r$dataset@tools$avg.b.cells_PE, aes(CTRL, STIM)) + geom_point() + ggtitle("B Cells") 
  #   Seurat::LabelPoints(plot = p1, points = result, repel = TRUE, xnudge = 0.2, ynudge = 0.2)
  # })
  
  ## -- Heatmap -- ##
  output$Heatmap <- renderPlot({r$dataset ; Seurat::DoHeatmap(r$dataset, features = heatmap()$gene, group.by = "Condition") + NoLegend()})
  output$Heatmap_feature <- renderPrint({print(r$dataset@commands[["FindAllMarkers"]])})
  

  #################################################################################################
  ## -- Enrichissement de gène -- ##
  reactive(input$numSelector)
  reactive(input$hallmark_order)
  reactive(input$hallmark_order_vln)
  reactive(input$metadata_order_vln)
  reactive(input$hallmark_order_X)
  reactive(input$hallmark_order_Y)
  reactive(input$metadata_order_density)
  reactive(input$hallmark_order_RP)
  reactive(input$metadata_group_RP)
  reactive(input$metadata_facet_RP)
  

  
  output$hallmark_Heatmap <- renderPlot({
    gc();gc();gc();gc()
    r$dataset@meta.data$active.idents <- r$dataset@active.ident
    #r$dataset@tools$meta_variable <- c("seurat_clusters",  "Phénotype", "Phase", "orig.ident") # "clonotype_id","Condition", "Greffe"
    dittoSeq::dittoHeatmap(r$dataset, genes = NULL, metas = input$numSelector, heatmap.colors = rev(colorblind_vector(50)),annot.by = r$dataset@tools$meta_variable, cluster_cols = F, fontsize = 10, order.by = input$hallmark_order)
    rm(list = ls())
    
  })
  output$hallmark_VlnPlot <- renderPlot({dittoSeq::dittoPlot(r$dataset , input$hallmark_order_vln, group.by = input$metadata_order_vln, legend.show = FALSE) + theme(title = element_text(size=20), axis.text = element_text(size=15)) +  ylab(label = "Score") })
  output$hallmark_HD <- renderPlot({dittoSeq::dittoScatterHex(r$dataset,x.var = input$hallmark_order_X, y.var = input$hallmark_order_Y, do.contour = TRUE, split.by =  input$metadata_order_density) + theme_classic() + scale_fill_gradientn(colors = rev(colorblind_vector(11))) + geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2)})
  output$hallmark_RP <- renderPlot({
    ES2 <- data.frame(r$dataset[[]], Seurat::Idents(r$dataset))
    colnames(ES2)[ncol(ES2)] <- "cluster"
    escape::ridgeEnrichment(ES2, gene.set = input$hallmark_order_RP, group = input$metadata_group_RP, facet = input$metadata_facet_RP, add.rug = TRUE)
  })
  output$hallmark_PCA <- renderPlot({
    ES2 <- data.frame(r$dataset[[]], Seurat::Idents(r$dataset))
    ES2 <- ES2[which(colnames(ES2) %in% r$dataset@tools$hallmarks)]
    #PCA <- performPCA(enriched = ES2, groups = colnames(singlet@meta.data))
    #pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE)
    #pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = FALSE, facet = "cluster") 
    masterPCAPlot(ES2, PCx = "PC1", PCy = "PC2", top.contribution = 10)
  })  
  
  #################################################################################################
  heatmap <- reactive({
    top10 <- r$dataset@commands[["FindAllMarkers"]] %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    return(top10)
  })
  FeaturesVariable <- reactive({
    annotations <- read.csv("document/annotation_FindAllMarkers.csv")
    fv <- annotations[which(annotations$gene_name %in% heatmap()$gene),] 
    fv <- fv[order(fv$gene_name),] 
    return(fv)
  })
  ## -- The 50 most highly variable genes -- ##
  output$top50 <- renderPlot({
    Seurat::VariableFeaturePlot(r$dataset)
    Seurat::LabelPoints(plot = Seurat::VariableFeaturePlot(r$dataset), points = head(Seurat::VariableFeatures(r$dataset), 50), repel = TRUE)
  })
  output$Variable_feature <- renderPrint({print(FeaturesVariable())})
  
}


## To be copied in the UI
# mod_Patient_Expression_ui("Expression_ui_1")

## To be copied in the server
# callModule(mod_Patient_Expression_server, "Expression_ui_1")