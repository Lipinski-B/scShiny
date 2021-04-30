# server.R
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
shinyServer(function(input, output, session) {
  #################################################################################################
  observeEvent(input$actBtnPatient,{
    singlet <<- seurat_subset(singlet, input$Subgroup, tosub())
  })
  observeEvent(input$resetPatient,{
    singlet <<- all
  })

  heatmap <- reactive({
    top10 <- singlet@commands[["FindAllMarkers"]] %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    return(top10)
  })
  FeaturesVariable <- reactive({
    annotations <- read.csv("/home/boris/Bureau/scShiny/annotation_FindAllMarkers.csv")
    fv <- annotations[which(annotations$gene_name %in% heatmap()$gene),] 
    fv <- fv[order(fv$gene_name),] 
    return(fv)
  })
  
  group <- reactive({
    v <- c()
    for(i in 1:length(metadata)){if(metadata[i] %in% input$Groupes){v = c(v, as.character(metadata[i]))}}
    return(v)
  })
  split <- reactive({
    v <- c()
    for(i in 1:length(metadata)){if(metadata[i] %in% input$Splites){v = c(v, as.character(metadata[i]))}}
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
    #singlet <- data()
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
  
  sub <- reactive({
    v <- c()
    for(i in 1:length(metadata)){if(metadata[i] %in% input$Subgroup){v = c(v, as.character(metadata[i]))}}
    return(v)
  })
  tosub <- reactive({
    for (l in sub()) {
      Idents(singlet) <- l
      tokeep <- levels(Idents(singlet))
      t <- c()
      for (i in tokeep){t <- c(t, as.character(input[[paste0("p",i)]]))}
      tokeep <- tokeep[tokeep %in% t]
    }
    return(tokeep)
  })
  
  singlet2 <- reactive({
    singlet2 <- singlet
    singlet2 <- RunUMAP(singlet2, reduction = "pca", dims = 1:40, n.components = 3L)
    singlet2 <- RunTSNE(singlet2, reduction = "pca", dims = 1:40, dim.embed = 3)
    return(singlet2)
  })
  
  D3plot <- reactive({
    plot.data <- FetchData(object = singlet2(), vars = c(meta_variable, "PC_1", "PC_2", "PC_3", "tSNE_1", "tSNE_2", "tSNE_3", "UMAP_1", "UMAP_2", "UMAP_3"))
    plot.data$label <- paste(rownames(plot.data))
    return(plot.data)
  })
  rv <- reactiveValues(
    hallmark = all@tools$hallmarks,
    order = NULL
  )

  output$Dynamic_Sub_Spe <- renderUI({
    choice <- list()
    for(i in 1:length(metadata)){
      for(j in 1:length(List[[metadata[i]]])){
        if(metadata[i] %in% input$Subgroup){
          choice[[List[[metadata[i]]][j]]] <- list(checkboxGroupInput(inputId = paste0("p",List[[metadata[i]]][j]), label = NULL, choices = List[[metadata[i]]][j]))  
        }
      }
      fluidRow()
    }
    return(choice)
  })
  output$Dynamic_Group_Spe <- renderUI({
    choice <- list()
    for(i in 1:length(metadata)){
      for(j in 1:length(List[[metadata[i]]])){
        if(metadata[i] %in% input$Groupes){
          choice[[List[[metadata[i]]][j]]] <- list(checkboxGroupInput(inputId = List[[metadata[i]]][j], label = NULL, choices = List[[metadata[i]]][j]))  
        }
      }
      fluidRow()
    }
    return(choice)
  })
  output$Dynamic_Split_Spe <- renderUI({
    choice <- list()   
    for(i in 1:length(metadata)){
      for(j in 1:length(List[[metadata[i]]])){
        if(metadata[i] %in% input$Splites){
          choice[[List[[metadata[i]]][j]]] <- list(checkboxGroupInput(inputId = paste0("s",List[[metadata[i]]][j]), label = NULL, choices = List[[metadata[i]]][j]))  
        }
      }
      fluidRow()
    }
    return(choice)
  })
  
  D3feature <- reactive({
    gene_feature <- input$in7
    singlet2 <- singlet2()
    all_feature <- rownames(as.matrix(singlet2[["RNA"]]@counts))
    plot.data <- FetchData(object = singlet2, vars = c(all_feature, "PC_1", "PC_2", "PC_3", "tSNE_1", "tSNE_2", "tSNE_3", "UMAP_1", "UMAP_2", "UMAP_3"), slot = 'data')
    plot.data$changed <- ifelse(test = plot.data[[gene_feature]] <1, yes = plot.data[[gene_feature]], no = 1)
    plot.data$label <- paste(rownames(plot.data)," - ", plot.data[[gene_feature]], sep="")
    return(plot.data)
  })
  output$Dfeature_pca = renderUI({if(length(input$in7)>0){plotlyOutput("DFeaturePlot_PCA", width = "100%",  height = "800px")}})
  output$Dfeature_umap = renderUI({if(length(input$in7)>0){plotlyOutput("DFeaturePlot_UMAP", width = "100%",  height = "800px")}})
  output$Dfeature_tsne = renderUI({if(length(input$in7)>0){plotlyOutput("DFeaturePlot_TSNE", width = "100%",  height = "800px")}})
  
  output$DFeaturePlot_PCA <- renderPlotly({
    plot_ly(data = D3feature(), x = ~PC_1, y = ~PC_2, z = ~PC_3, 
            color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
            colors = c('darkgreen', 'red'), opacity = .5,
            type = "scatter3d", mode = "markers",
            marker = list(size = 3, width=2), 
            text=~label, hoverinfo="text"
    )# %>% layout(title=gene_feature)
  })
  output$DFeaturePlot_UMAP<- renderPlotly({
    plot_ly(data = D3feature(), x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
            color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
            colors = c('darkgreen', 'red'), opacity = .5,
            type = "scatter3d", mode = "markers",
            marker = list(size = 3, width=2), 
            text=~label, hoverinfo="text"
    )# %>% layout(title=gene_feature)
  })
  output$DFeaturePlot_TSNE<- renderPlotly({
    plot_ly(data = D3feature(), x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
            color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
            colors = c('darkgreen', 'red'), opacity = .5,
            type = "scatter3d", mode = "markers",
            marker = list(size = 3, width=2), 
            text=~label, hoverinfo="text"
    ) #%>% layout(title=gene_feature)
  })
  
  
  feature <- reactive({row.names(as.matrix(singlet[["RNA"]]@counts))})
  output$variables = renderUI({selectInput('in6', 'Choices', feature(), multiple=TRUE, selectize=TRUE)})
  output$Dvariables = renderUI({selectInput('in7', 'Choices', feature(), selectize=TRUE, selected = NULL)})
  output$feature_pca = renderUI({if(length(input$in6)>0){plotOutput("FeaturePlot_PCA", width = "100%",  height = "650px")}})
  output$feature_umap = renderUI({if(length(input$in6)>0){plotOutput("FeaturePlot_UMAP", width = "100%",  height = "650px")}})
  output$feature_tsne = renderUI({if(length(input$in6)>0){plotOutput("FeaturePlot_TSNE", width = "100%",  height = "650px")}})
  
  output$DDvariables = renderUI({selectInput('in8', 'Choices', feature(), selectize=TRUE, multiple=TRUE, selected = NULL)})
  output$DDorder_trajectory = renderUI({if(length(input$in8)>0){plotOutput("order_trajectory", width = "100%",  height = "1000px")}})
  
  #################################################################################################
  ## -- FeaturePlot -- ## 
  output$FeaturePlot_PCA <- renderPlot({FeaturePlot(tokeep(), features = input$in6, reduction = "pca", split.by = split(), pt.size = 2, combine = T) & NoAxes() & NoLegend() & theme(title = element_text(size=20))})
  output$FeaturePlot_UMAP <- renderPlot({FeaturePlot(tokeep(), features = input$in6, reduction = "umap", split.by = split(), pt.size = 2, combine = T) & NoAxes() & NoLegend()})
  output$FeaturePlot_TSNE <- renderPlot({FeaturePlot(tokeep(), features = input$in6, reduction = "tsne", split.by = split(), pt.size = 2, combine = T) & NoAxes() & NoLegend()})
  output$feature_other <- renderPlotly({FeaturePlot(singlet, features = "CD19", interactive = T)})
  
  
  ## -- PCA -- ##
  output$PCA <- renderPlot({
    plots <- PCAPlot(object = tokeep(), group.by = group(), split.by = split(), label.size = 0.0, pt.size = 2)
    plots & theme(title = element_text(size=20),
                  legend.position = "top",
                  legend.title = element_text(size=15),
                  legend.text = element_text(size=15)
                  ) & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 6))) & xlab(label = paste0("PCA 1 : ", round(Stdev(singlet[["pca"]])[1],2), " %")) & ylab(label = paste0("PCA 2 : ", round(Stdev(singlet[["pca"]])[2],2), " %"))
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
  output$Heatmap <- renderPlot({DoHeatmap(singlet, features = heatmap()$gene, group.by = "Condition") + NoLegend()})
  output$Heatmap_feature <- renderPrint({print(singlet@commands[["FindAllMarkers"]])})

  
  ## -- Mitochondrie Figure -- ##
  output$MT_VlnPlot <- renderPlot({
    singlet@tools$mitochondrie_all + VlnPlot(singlet, features = "nFeature_RNA", group.by = "Condition") +  VlnPlot(singlet, features = "nCount_RNA", group.by = "Condition") + VlnPlot(singlet, features = "percent.mt", group.by = "Condition")
  })
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
  
  
  # -- Enrichissement de gène -- ##
  observeEvent(input$actBtnVisualisation,{
    if(length(input$numSelector)==0){rv$hallmark <- all@tools$hallmarks}
    else{rv$hallmark <- input$numSelector}
    rv$order <- input$hallmark_order
  })
  observeEvent(input$Subsets,{
    rv$hallmark <- all@tools$hallmarks
    rv$order <-NULL}, ignoreNULL = FALSE)
  output$hallmark_Heatmap <- renderPlot({
    singlet@meta.data$active.idents <- singlet@active.ident
    
    if(is.null(rv$order)){
      dittoHeatmap(singlet, genes = NULL, metas = all@tools$hallmarks, heatmap.colors = rev(colorblind_vector(50)),
                   annot.by = singlet@tools$meta_variable, cluster_cols = T, fontsize = 12)
    } else{
      dittoHeatmap(singlet, genes = NULL, metas = rv$hallmark, heatmap.colors = rev(colorblind_vector(50)),
                   annot.by = singlet@tools$meta_variable, cluster_cols = F, fontsize = 12, order.by = rv$order)
    }  
  })
  output$hallmark_VlnPlot <- renderPlot({
    dittoPlot(singlet, input$hallmark_order_vln, group.by = input$metadata_order_vln) + scale_fill_manual(values = colorblind_vector(20))})
  output$hallmark_HD <- renderPlot({
    dittoScatterHex(singlet,x.var = input$hallmark_order_X, y.var = input$hallmark_order_Y, do.contour = TRUE, split.by =  input$metadata_order_density) + 
      theme_classic() + scale_fill_gradientn(colors = rev(colorblind_vector(11))) + geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2)  
  })
  output$hallmark_RP <- renderPlot({
    ES2 <- data.frame(singlet[[]], Idents(singlet))
    colnames(ES2)[ncol(ES2)] <- "cluster"
    ridgeEnrichment(ES2, gene.set = input$hallmark_order_RP, group = input$metadata_group_RP, facet = input$metadata_facet_RP, add.rug = TRUE)
  })
  #output$hallmark_PCA <- renderPlot({
  #ES2 <- data.frame(singlet[[]], Idents(singlet))
  #PCA <- performPCA(enriched = ES2, groups = c("cluster", "SingleR.calls"))
  #pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE)
  #pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = FALSE, facet = "cluster") 
  #masterPCAPlot(ES2, PCx = "PC1", PCy = "PC2", top.contribution = 10)
  #})
  
  
  ## -- Monocle -- ##
  output$cell_trajectory <- renderPlot({
    plot_cell_trajectory(data, color_by = input$color_trajectory)
  })
  output$order_trajectory <- renderPlot({
    data_expressed_genes <-  row.names(subset(fData(data), num_cells_expressed >= 10))
    data_filtered <- data[data_expressed_genes,]
    my_genes <- row.names(subset(fData(data_filtered), gene_short_name %in% input$in8))
    cds_subset <- data_filtered[my_genes,]
    p <- plot_genes_in_pseudotime(cds_subset, color_by = input$color_order)
    p + theme(
      text = element_text(size=20),
      legend.title = element_text(color = "black", size = 14),
      legend.text = element_text(color = "black", size = 14),
      axis.text=element_text(size=14),
      axis.title=element_text(size=14)
    ) 
  })
  
  #output$dataTable = DT::renderDataTable({singlet@meta.data}, options = list(
  #  scrollY = '700px', paging = FALSE,scrollX = TRUE
  #))

  output$dataTable2 = renderImage({
    list(src = '/home/boris/Documents/analyse/sunburst2.png', contentType = 'image/png',width = 1100, height = 800,
         alt = "Alternate text")
  }, deleteFile = FALSE)
  
  output$dataTable = renderPlotly({
    plot_ly(singlet@tools$sunburst, ids = ~ids ,labels = ~labels, parents = ~parents, values = ~values, marker = list(colors = c( "#BEBADA", "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462")),
            type = 'sunburst', branchvalues = 'total', hoverinfo = "text", hovertext = paste(singlet@tools$sunburst$labels, ":", round((singlet@tools$sunburst$values/singlet@tools$sunburst$values[1])*100,2),"%"))
  })

})
