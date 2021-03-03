# server.R

shinyServer(function(input, output, session) {

  load(file = "/home/boris/Documents/analyse/singlet_hFL_180008B.RData")
  
  #################################################################################################
  heatmap <- reactive({
    top10 <- singlet@commands[["FindAllMarkers"]] %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
    return(top10)
  })
  FeaturesVariable <- reactive({
    annotations <- read.csv("/home/boris/Documents/analyse/annotation_FindAllMarkers.csv")
    fv <- annotations[which(annotations$gene_name %in% heatmap()$gene),] 
    fv <- fv[order(fv$gene_name),] 
    return(fv)
  })
  
  group <- reactive({
    v <- c()
    for(i in 1:length(colnames(singlet@meta.data))){if(length(input[[colnames(singlet@meta.data)[i]]])== 1){v = c(v, as.character(input[[colnames(singlet@meta.data)[i]]]))}}
    return(v)
  })
  split <- reactive({
    v <- c()
    for(i in 1:length(colnames(singlet@meta.data))){if(length(input[[paste0("s",colnames(singlet@meta.data)[i])]])== 1){v = c(v, as.character(input[[paste0("s",colnames(singlet@meta.data)[i])]]))}}
    if(length(input$sclonotype) == 1){v = c(v, 'nbs_clonotype')}
    return(v)
  })
  clonotype <- reactive({
    clone=as.data.frame(singlet$clonotype_id[which(singlet$clonotype_id<=input$NBS_Clonotype)])
    rownames(clone)=rownames(as.matrix(singlet$clonotype_id[which(singlet$clonotype_id<=input$NBS_Clonotype)]))
    colnames(clone)="nbs_clonotype"
    singlet <- AddMetaData(object=singlet, metadata = clone)
    return(singlet)
  })
  tokeep <- reactive({
    singlet2 <- singlet
    for (l in split()) {
      
      Idents(singlet2) <- l
      tokeep <- levels(Idents(singlet2))
      t <- c()
      for (i in tokeep){t <- c(t, input[[paste0("s",i)]])}      
      tokeep <- tokeep[tokeep %in% t]
      singlet2 <- subset(singlet2, idents = tokeep)
      Idents(singlet2)<-"seurat_clusters"
      if(length(input$sclonotype) == 1){
        singlet2 <- clonotype()
      }
    }
    
    for (k in group()) {
      Idents(singlet2)<- k
      tokeep <- levels(Idents(singlet2))
      t <- c()
      for (i in tokeep){t <- c(t, input[[i]])}      
      tokeep <- tokeep[tokeep %in% t]
      singlet2 <- subset(singlet2, idents = tokeep)
      Idents(singlet2)<-"seurat_clusters"
    }
    
    return(singlet2)
  })
  D3plot <- reactive({
    singlet2 <- singlet
    singlet2 <- RunUMAP(singlet2, reduction = "pca", dims = 1:40, n.components = 3L)
    singlet2 <- RunTSNE(singlet2, reduction = "pca", dims = 1:40, dim.embed = 3)
    meta_variable <- c("seurat_clusters", "HTO_maxID", "SingleR.calls", "clonotype_id","chain", "v_gene", "d_gene", "j_gene","c_gene", "cdr3", "Phase")
    plot.data <- FetchData(object = singlet2, vars = c(meta_variable, "PC_1", "PC_2", "PC_3", "tSNE_1", "tSNE_2", "tSNE_3", "UMAP_1", "UMAP_2", "UMAP_3"))
    plot.data$label <- paste(rownames(plot.data))
    return(plot.data)
  })
  
  output$Dynamic_Group <- renderUI({
    List <- list()
    for(i in 1:length(colnames(singlet@meta.data))){
      if (length(levels(as.factor(singlet@meta.data[[i]]))) > 1 && length(levels(as.factor(singlet@meta.data[[i]]))) < 25 && is.numeric(levels(as.factor(singlet@meta.data[[1]])))==F ){
        List[[colnames(singlet@meta.data)[i]]] <- list(column(3,align="left",checkboxGroupInput(inputId = colnames(singlet@meta.data)[i], label = NULL, choices = colnames(singlet@meta.data)[i])))
      }
    }
    return(List)                     
  })
  output$Dynamic_Group_Spe <- renderUI({
    List <- list()   
    for(i in 1:length(colnames(singlet@meta.data))){
      for(j in 1:length(levels(as.factor(singlet@meta.data[[i]])))){
        if(length(input[[colnames(singlet@meta.data)[i]]]) == 1){
          List[[levels(as.factor(singlet@meta.data[[i]]))[j]]] <- list(column(3,align="left",checkboxGroupInput(inputId = levels(as.factor(singlet@meta.data[[i]]))[j], label = NULL, choices = levels(as.factor(singlet@meta.data[[i]]))[j])))  
        }
      }
      fluidRow()
    }
    return(List)
  })
  output$Dynamic_Split <- renderUI({
    List <- list()
    for(i in 1:length(colnames(singlet@meta.data))){
      if (length(levels(as.factor(singlet@meta.data[[i]]))) > 1 && length(levels(as.factor(singlet@meta.data[[i]]))) < 25 && is.numeric(levels(as.factor(singlet@meta.data[[1]])))==F ){
        List[[colnames(singlet@meta.data)[i]]] <- list(column(3,align="left",checkboxGroupInput(inputId = paste0("s",colnames(singlet@meta.data)[i]), label = NULL, choices = colnames(singlet@meta.data)[i])))
      }
    }
    return(List)                     
  })
  output$Dynamic_Split_Spe <- renderUI({
    List <- list()   
    for(i in 1:length(colnames(singlet@meta.data))){
      for(j in 1:length(levels(as.factor(singlet@meta.data[[i]])))){
        if(length(input[[paste0("s",colnames(singlet@meta.data)[i])]]) == 1){
          List[[levels(as.factor(singlet@meta.data[[i]]))[j]]] <- list(column(3,align="left",checkboxGroupInput(inputId =  paste0("s",levels(as.factor(singlet@meta.data[[i]]))[j]), label = NULL, choices = levels(as.factor(singlet@meta.data[[i]]))[j])))  
        }
      }
      fluidRow()
    }
    return(List)
  })
  
  feature <- reactive({row.names(as.matrix(singlet[["RNA"]]@counts))})
  output$variables = renderUI({selectInput('in6', 'Choices', feature(), multiple=TRUE, selectize=TRUE)})
  output$feature_pca = renderUI({if(length(input$in6)>0){plotOutput("FeaturePlot_PCA", width = "100%",  height = "650px")}})
  output$feature_umap = renderUI({if(length(input$in6)>0){plotOutput("FeaturePlot_UMAP", width = "100%",  height = "650px")}})
  output$feature_tsne = renderUI({if(length(input$in6)>0){plotOutput("FeaturePlot_TSNE", width = "100%",  height = "650px")}})
  
  #################################################################################################
  ## -- FeaturePlot -- ## 
  output$FeaturePlot_PCA <- renderPlot({FeaturePlot(tokeep(), features = input$in6, reduction = "pca", split.by = split(), pt.size = 2, combine = T) & NoAxes() & NoLegend()})
  output$FeaturePlot_UMAP <- renderPlot({FeaturePlot(tokeep(), features = input$in6, reduction = "umap", split.by = split(), pt.size = 2, combine = T) & NoAxes() & NoLegend()})
  output$FeaturePlot_TSNE <- renderPlot({FeaturePlot(tokeep(), features = input$in6, reduction = "tsne", split.by = split(), pt.size = 2, combine = T) & NoAxes() & NoLegend()})
  output$feature_other <- renderPlotly({FeaturePlot(singlet, features = "CD19", interactive = T)})
  
  ## -- PCA -- ##
  output$PCA <- renderPlot({
    plots <- PCAPlot(object = tokeep(), group.by = group(), split.by = split(), label.size = 0.0, pt.size = 2)
    plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4))) & xlab(label = paste0("PCA 1 : ", round(Stdev(singlet[["pca"]])[1],2), " %")) & ylab(label = paste0("PCA 2 : ", round(Stdev(singlet[["pca"]])[2],2), " %"))
  })
  
  output$D_PCA <- renderPlotly({
    plot_ly(data = D3plot(), x = ~PC_1, y = ~PC_2, z = ~PC_3, 
            color = singlet@meta.data[[input$metadata]], 
            type = "scatter3d", mode = "markers", 
            marker = list(size = 3, width=2),
            text=~singlet@meta.data[[input$metadata]], hoverinfo="text") %>% layout(title=input$metadata,
                                                      scene = list(
                                                        xaxis = list(title = paste0("PCA 1 : ", round(Stdev(singlet[["pca"]])[1],2), " %")),
                                                        yaxis = list(title = paste0("PCA 2 : ", round(Stdev(singlet[["pca"]])[2],2), " %")),
                                                        zaxis = list(title = paste0("PCA 3 : ", round(Stdev(singlet[["pca"]])[3],2), " %"))
                                                      ))
  })
  
  ## -- UMAP -- ##
  output$UMAP <- renderPlot({
    plots <- UMAPPlot(object = tokeep(), group.by = group(), split.by = split(), label.size = 0.0, pt.size = 2)
    plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4))) 
  })
  
  output$D_UMAP <- renderPlotly({
    plot_ly(data = D3plot(), x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
            color = singlet@meta.data[[input$metadata]],
            type = "scatter3d", mode = "markers", 
            marker = list(size = 3, width=2),
            text=~singlet@meta.data[[input$metadata]], hoverinfo="text") %>% layout(title=input$metadata) 
  })
  
  ## -- TSNE -- ##
  output$TSNE <- renderPlot({
    plots <- TSNEPlot(object = tokeep(), group.by = group(), split.by = split(), label.size = 0.0, pt.size = 2)
    plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4)))
  })
  
  output$D_TSNE <- renderPlotly({
    data2 <- D3plot()
    plot_ly(data = data2, x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
            color = singlet@meta.data[[input$metadata]], 
            type = "scatter3d", mode = "markers", 
            marker = list(size = 3, width=2),
            text=~singlet@meta.data[[input$metadata]], hoverinfo="text") %>% layout(title=input$metadata)
  })
  
  ## -- Heatmap -- ##
  output$Heatmap <- renderPlot({DoHeatmap(singlet, features = heatmap()$gene, group.by = "HTO_maxID") + NoLegend()})
  output$Heatmap_feature <- renderPrint({print(singlet@commands[["FindAllMarkers"]])})
  
  ## -- Mitochondrie Figure -- ##
  output$MT_VlnPlot <- renderPlot({VlnPlot(singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)})
  output$MT_FeatureScatter <- renderPlot({FeatureScatter(singlet, feature1 = "nCount_RNA", feature2 = "percent.mt") + FeatureScatter(singlet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")})
  output$MT_FeatureScatter2 <- renderPlot({FeatureScatter(singlet, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "HTO_classification") + FeatureScatter(singlet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "HTO_classification")})
  
  ## -- The 50 most highly variable genes -- ##
  output$top50 <- renderPlot({
    VariableFeaturePlot(singlet)
    LabelPoints(plot = VariableFeaturePlot(singlet), points = head(VariableFeatures(singlet), 50), repel = TRUE)
  })
  output$Variable_feature <- renderPrint({print(FeaturesVariable())})
  
  # -- PCA Figure -- ##
  output$PCA10 <- renderPrint({print(singlet[["pca"]], dims = 1:10, nfeatures = 10)})
  output$ElbowPlot <- renderPlot({ElbowPlot(singlet, ndims = 50, reduction = "pca")})
  output$DimHeatmap <- renderPlot({DimHeatmap(singlet, dims = 1:10, cells = 100, balanced = TRUE)})
  output$VizDimLoadings <- renderPlot({VizDimLoadings(singlet, dims = 1:5, reduction = "pca")})
  output$JackStrawPlot <- renderPlot({JackStrawPlot(singlet, dims = 1:15)})
  
  output$dataTable = DT::renderDataTable(as.matrix(singlet[["RNA"]]@counts))
})