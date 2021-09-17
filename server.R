# server.R
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
jscode <- "shinyjs.refresh = function() { location.reload(); }"
shinyServer(function(input, output, session) {
  #################################################################################################
  sortie <- eventReactive(c(input$actBtnPatient1,input$actBtnPatient, input$resetPatient) ,{singlet})
  
  observeEvent(input$actBtnPatient1,{    
    shinybusy::show_modal_spinner(spin = "semipolar",color = "deepskyblue",text = "Please wait...")
    if(nchar(as.character(input$patient))>0){
      load(file = paste0("www/",input$patient,"/",input$patient,".RData"))
      singlet <<- singlet
    }
    shinyWidgets::sendSweetAlert(session = session,title = "Done !",type = "success")
    shinyjs::runjs("window.scrollTo(0, 50)")
    shinybusy::remove_modal_spinner()
  })
  observeEvent(input$actBtnPatient,{    
    shinybusy::show_modal_spinner(spin = "semipolar",color = "deepskyblue",text = "Please wait...")
    
    if(nchar(as.character(input$patient))>0){
      load(file = paste0("www/",input$patient,"/",input$patient,".RData"))
      singlet <<- singlet
    }
    sortie()
    
    if(length(input$Subgroup)>0){singlet <<- seurat_subset(singlet, input$Subgroup, c(tosub()))}
    if(length(input$subFeatures)>0){singlet <<- gene_subset(singlet, input$Svariables, as.integer(input$Seuil_variables))}
    if(nchar(as.character(input$maximum))>0 || nchar(as.character(input$percent_mt))>0){singlet <<- QC_subset(singlet, maximum_sub=input$maximum, percent_mt_sub=input$percent_mt)}
    
    if(length(input$combo)>0){
      load(file = paste0("www/singlet_",input$patient,"_",input$combo,".RData"))
      singlet <<- singlet
    }
    
    shinyWidgets::sendSweetAlert(session = session,title = "Done !",text = "Le subsample a bien été créé !",type = "success")
    shinybusy::remove_modal_spinner()
  })
  
  observeEvent(input$resetPatient,{
    singlet <<- singlet
    shinyWidgets::sendSweetAlert(session = session,title = "Done !",text = "Tous les subsamples ont été supprimé !",type = "success")
  })
  
  heatmap <- reactive({
    top10 <- sortie()@commands[["FindAllMarkers"]] %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    return(top10)
  })
  FeaturesVariable <- reactive({
    annotations <- read.csv("document/annotation_FindAllMarkers.csv")
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
    sortie()
    singlet2 <- singlet
    for (l in split()) {
      Seurat::Idents(singlet2) <- l
      tokeep <- levels(Seurat::Idents(singlet2))
      t <- c()
      for (i in tokeep){t <- c(t, input[[paste0("s",i)]])}      
      tokeep <- tokeep[tokeep %in% t]
      singlet2 <- subset(singlet2, idents = tokeep)
      Seurat::Idents(singlet2)<-"seurat_clusters"
      if(length(input$sclonotype) == 1){
        singlet2 <- clonotype()
      }
    }
    for (k in group()) {
      Seurat::Idents(singlet2)<- k
      tokeep <- levels(Seurat::Idents(singlet2))
      t <- c()
      for (i in tokeep){t <- c(t, input[[i]])}      
      tokeep <- tokeep[tokeep %in% t]
      singlet2 <- subset(singlet2, idents = tokeep)
      Seurat::Idents(singlet2)<-"seurat_clusters"
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
      Seurat::Idents(singlet) <- l
      tokeep <- levels(Seurat::Idents(singlet))
      t <- c()
      for (i in tokeep){t <- c(t, as.character(input[[paste0("p",i)]]))}
      tokeep <- tokeep[tokeep %in% t]
    }
    return(tokeep)
  })
  
  singlet2 <- reactive({
    sortie()
    singlet2 <- singlet
    singlet2 <- RunUMAP(singlet2, reduction = "pca", dims = 1:40, n.components = 3L)
    singlet2 <- RunTSNE(singlet2, reduction = "pca", dims = 1:40, dim.embed = 3)
    return(singlet2)
  })
  
  D3plot <- reactive({
    sortie()
    plot.data <- FetchData(object = singlet2(), vars = c(singlet@tools$meta_variable, "PC_1", "PC_2", "PC_3", "tSNE_1", "tSNE_2", "tSNE_3", "UMAP_1", "UMAP_2", "UMAP_3"))
    plot.data$label <- paste(rownames(plot.data))
    return(plot.data)
  })
  rv <- reactiveValues(
    hallmark = singlet@tools$hallmarks,
    order = NULL
  )
  
  output$nb_tot_cell <- renderText({return(paste("Number of cells : \t", as.character(as.numeric(sum(table(sortie()@meta.data$Phénotype))))))})
  output$nb_b_cell <- renderText({return(paste("B cells : \t", as.character(as.numeric(table(sortie()@meta.data$Phénotype)["B-cells"]))))})
  output$nb_other_cell <- renderText({return(paste("\nOther cells : \t", as.character(as.numeric(sum(table(sortie()@meta.data$Phénotype))) - as.numeric(table(sortie()@meta.data$Phénotype)["B-cells"]))))})
  output$nb_Controle_cell <- renderText({return(paste("\nControle : \t", as.character(as.numeric(table(sortie()@meta.data$Condition)["Pré-greffe"]))))})
  output$nb_Excipient_cell <- renderText({return(paste("\nExcipient : \t", as.character(as.numeric(table(sortie()@meta.data$Condition)["Excipient"]))))})
  output$nb_RCHOP_cell <- renderText({return(paste("\nRCHOP : \t", as.character(as.numeric(table(sortie()@meta.data$Condition)["RCHOP"]))))})
  output$nb_clonotype_cell2 <- renderText({
    sortie()
    return(paste0("Numbre total de cellule détectée : \t", as.character(as.numeric(sum(table(singlet@meta.data$Phénotype)))),
                  "\n\nB cells : \t", as.character(as.numeric(table(singlet@meta.data$Phénotype)["B-cells"])),
                  "\nOther cells : \t", as.character(as.numeric(sum(table(singlet@meta.data$Phénotype))) - as.numeric(table(singlet@meta.data$Phénotype)["B-cells"])),
                  
                  "\n\nControle : \t", as.character(as.numeric(table(singlet@meta.data$Condition)["Pré-greffe"])),
                  "\nExcipient : \t", as.character(as.numeric(table(singlet@meta.data$Condition)["Excipient"])),
                  "\nRCHOP : \t", as.character(as.numeric(table(singlet@meta.data$Condition)["RCHOP"])),
                  
                  "\n\nClonotype majoritaire :", 
                  "\n Proportion : \t", round(singlet@tools$vloupe$proportion[1],3)*100,"% \n Type : \t", singlet@tools$vloupe$type[1],"\n Isotype : \t", singlet@tools$vloupe$igh_c_genes[1], 
                  "\n Heavy : \t", singlet@tools$vloupe$V_lourde[1], "/", singlet@tools$vloupe$D_lourde[1], "/", singlet@tools$vloupe$J_lourde[1], 
                  "\n Light : \t", singlet@tools$vloupe$V_legere[1], "/", singlet@tools$vloupe$J_legere[1], "\n"))
  })
  output$nb_clonotype_cell <- renderText({
    sortie()
    return(paste0("Clonotype majoritaire :", 
                  "\n Proportion : ", round(singlet@tools$vloupe$proportion[1],3)*100,"% \n Type : ", singlet@tools$vloupe$type[1],"\n Isotype : ", singlet@tools$vloupe$igh_c_genes[1], 
                  "\n Heavy : ", singlet@tools$vloupe$V_lourde[1], "/", singlet@tools$vloupe$D_lourde[1], "/", singlet@tools$vloupe$J_lourde[1], 
                  "\n Light : ", singlet@tools$vloupe$V_legere[1], "/", singlet@tools$vloupe$J_legere[1]))
  })
  output$H_clonotype_cell <- renderText({
    sortie()
    return(paste0("Heavy : ", singlet@tools$vloupe$V_lourde[1], "/", singlet@tools$vloupe$D_lourde[1], "/", singlet@tools$vloupe$J_lourde[1]))
  })
  output$L_clonotype_cell <- renderText({
    sortie()
    return(paste0("Light : ", singlet@tools$vloupe$V_legere[1], "/", singlet@tools$vloupe$J_legere[1]))
  })
  
  output$Dynamic_Sub_Spe <- renderUI({
    choice <- list()
    for(i in 1:length(metadata)){
      for(j in 1:length(List[[metadata[i]]])){
        if(metadata[i] %in% input$Subgroup){choice[[List[[metadata[i]]][j]]] <- list(checkboxGroupInput(inputId = paste0("p",List[[metadata[i]]][j]), label = NULL, choices = List[[metadata[i]]][j]))}
      }
      fluidRow()
    }
    return(choice)
  })
  output$Dynamic_Group_Spe <- renderUI({
    choice <- list()
    for(i in 1:length(metadata)){
      for(j in 1:length(List[[metadata[i]]])){
        if(metadata[i] %in% input$Groupes){choice[[List[[metadata[i]]][j]]] <- list(checkboxGroupInput(inputId = List[[metadata[i]]][j], label = NULL, choices = List[[metadata[i]]][j]))}
      }
      fluidRow()
    }
    return(choice)
  })
  output$Dynamic_Split_Spe <- renderUI({
    choice <- list()   
    for(i in 1:length(metadata)){
      for(j in 1:length(List[[metadata[i]]])){
        if(metadata[i] %in% input$Splites){choice[[List[[metadata[i]]][j]]] <- list(checkboxGroupInput(inputId = paste0("s",List[[metadata[i]]][j]), label = NULL, choices = List[[metadata[i]]][j]))}
      }
      fluidRow()
    }
    return(choice)
  })
  
  D3feature <- reactive({
    gene_feature <- input$Dvariables
    singlet2 <- singlet2()
    all_feature <- rownames(as.matrix(singlet2[["RNA"]]@counts))
    plot.data <- FetchData(object = singlet2, vars = c(all_feature, "PC_1", "PC_2", "PC_3", "tSNE_1", "tSNE_2", "tSNE_3", "UMAP_1", "UMAP_2", "UMAP_3"), slot = 'data')
    plot.data$changed <- ifelse(test = plot.data[[gene_feature]] <1, yes = plot.data[[gene_feature]], no = 1)
    plot.data$label <- paste(rownames(plot.data)," - ", plot.data[[gene_feature]], sep="")
    return(plot.data)
  })
  output$Dfeature_pca = renderUI({if(length(input$Dvariables)>0){plotlyOutput("DFeaturePlot_PCA", width = "100%",  height = "800px")}})
  output$Dfeature_umap = renderUI({if(length(input$Dvariables)>0){plotlyOutput("DFeaturePlot_UMAP", width = "100%",  height = "800px")}})
  output$Dfeature_tsne = renderUI({if(length(input$Dvariables)>0){plotlyOutput("DFeaturePlot_TSNE", width = "100%",  height = "800px")}})
  
  output$DFeaturePlot_PCA <- renderPlotly({
    plot_ly(data = D3feature(), x = ~PC_1, y = ~PC_2, z = ~PC_3,color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
            colors = c('darkgreen', 'red'), opacity = .5,type = "scatter3d", mode = "markers", marker = list(size = 3, width=2),text=~label, hoverinfo="text")# %>% layout(title=gene_feature)
  })
  output$DFeaturePlot_UMAP<- renderPlotly({
    plot_ly(data = D3feature(), x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
            colors = c('darkgreen', 'red'), opacity = .5,type = "scatter3d", mode = "markers", marker = list(size = 3, width=2),text=~label, hoverinfo="text")# %>% layout(title=gene_feature)
  })
  output$DFeaturePlot_TSNE<- renderPlotly({
    plot_ly(data = D3feature(), x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3,color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
            colors = c('darkgreen', 'red'), opacity = .5,type = "scatter3d", mode = "markers", marker = list(size = 3, width=2),text=~label, hoverinfo="text") #%>% layout(title=gene_feature)
  })
  
  updateSelectizeInput(session, 'Svariables', choices = feature, server = TRUE)
  updateSelectizeInput(session, 'variables', choices = feature, server = TRUE)
  updateSelectizeInput(session, 'Dvariables', choices = feature, server = TRUE)
  updateSelectizeInput(session, 'DDvariables', choices = feature, server = TRUE)
  
  output$feature_pca = renderUI({if(length(input$variables)>0){plotOutput("FeaturePlot_PCA", width = "100%",  height = "500px")}})
  output$feature_umap = renderUI({if(length(input$variables)>0){plotOutput("FeaturePlot_UMAP", width = "100%",  height = "650px")}})
  output$feature_tsne = renderUI({if(length(input$variables)>0){plotOutput("FeaturePlot_TSNE", width = "100%",  height = "650px")}})
  
  output$DDorder_trajectory = renderUI({if(length(input$DDvariables)>0){plotOutput("order_trajectory", width = "100%",  height = "1000px")}})
  
  #################################################################################################
  ## -- FeaturePlot -- ## 
  output$FeaturePlot_PCA <- renderPlot({FeaturePlot(tokeep(), features = input$variables, reduction = "pca", split.by = split(), pt.size = 1, combine = T, ncol = 2) & NoAxes() & NoLegend() & theme(title = element_text(size=20))})
  output$FeaturePlot_UMAP <- renderPlot({FeaturePlot(tokeep(), features = input$variables, reduction = "umap", split.by = split(), pt.size = 1, combine = T) & NoAxes() & NoLegend()})
  output$FeaturePlot_TSNE <- renderPlot({FeaturePlot(tokeep(), features = input$variables, reduction = "tsne", split.by = split(), pt.size = 1, combine = T) & NoAxes() & NoLegend()})
  output$feature_other <- renderPlotly({FeaturePlot(singlet, features = "CD19", interactive = T)})
  
  ## -- PCA -- ##
  output$PCA <- renderPlot({
    Seurat::DimPlot(object = tokeep(), group.by = group(), split.by = split(), label.size = 0.0, pt.size = 2, reduction = 'pca') & theme(title = element_text(size=20),legend.position = "top",legend.title = element_text(size=10),legend.text = element_text(size=10)
    ) & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 6))) & xlab(label = paste0("PCA 1 : ", round(Seurat::Stdev(singlet[["pca"]])[1],2), " %")) & ylab(label = paste0("PCA 2 : ", round(Seurat::Stdev(singlet[["pca"]])[2],2), " %"))
  })
  output$D_PCA <- renderPlotly({
    plot_ly(data = D3plot(), x = ~PC_1, y = ~PC_2, z = ~PC_3,color = singlet@meta.data[[input$metadata]],type = "scatter3d", mode = "markers",marker = list(size = 3, width=2),
            text=~singlet@meta.data[[input$metadata]], hoverinfo="text") %>%layout(title=input$metadata,scene = list(xaxis = list(title = paste0("PCA 1 : ", round(Seurat::Stdev(singlet[["pca"]])[1],2), " %")),yaxis = list(title = paste0("PCA 2 : ", round(Seurat::Stdev(singlet[["pca"]])[2],2), " %")),zaxis = list(title = paste0("PCA 3 : ", round(Seurat::Stdev(singlet[["pca"]])[3],2), " %"))))
  })
  
  ## -- UMAP -- ##
  output$UMAP <- renderPlot({Seurat::DimPlot(object = tokeep(), group.by = group(), split.by = split(), label.size = 0.0, pt.size = 2, reduction = 'umap') & theme(legend.position = "top") & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4)))})
  output$D_UMAP <- renderPlotly({plot_ly(data = D3plot(), x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,color = singlet@meta.data[[input$metadata]],type = "scatter3d", mode = "markers",marker = list(size = 3, width=2),text=~singlet@meta.data[[input$metadata]], hoverinfo="text") %>% layout(title=input$metadata)})
  
  ## -- TSNE -- ##
  output$TSNE <- renderPlot({Seurat::DimPlot(object = tokeep(), group.by = group(), split.by = split(), label.size = 0.0, pt.size = 2, reduction = 'tsne') & theme(legend.position = "top") & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4)))})
  output$D_TSNE <- renderPlotly({plot_ly(data = D3plot(), x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3,color = singlet@meta.data[[input$metadata]],type =  "scatter3d", mode = "markers",marker = list(size = 3, width=2),text=~singlet@meta.data[[input$metadata]], hoverinfo="text") %>% layout(title=input$metadata)})
  
  ## -- Expression différentielle -- ##
  # RCHOP / Excipient
  output$DE_Heatmap_RE <- renderPlot({
    sortie()
    Seurat::Idents(singlet)<-"Condition"
    Seurat::DoHeatmap(singlet, cells = rownames(singlet@meta.data)[which(singlet@meta.data$Condition==c("Excipient","RCHOP"))], features = rownames(singlet@tools$DE_RE)[1:50], size = 3, assay = 'SCT', slot = "scale.data")
  })
  output$DE_RidgePlot_RE <- renderPlot({
    sortie()
    Seurat::Idents(singlet)<-"Condition"
    Seurat::RidgePlot(singlet, idents = c("Excipient","RCHOP"), features = rownames(singlet@tools$DE_RE)[1:12], ncol = 3, assay = 'RNA', slot = "data") & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
  })
  output$DE_VlnPlot_RE <- renderPlot({
    sortie()
    Seurat::Idents(singlet)<-"Condition"
    Seurat::VlnPlot(singlet, idents = c("Excipient","RCHOP"), features = rownames(singlet@tools$DE_RE)[1:12], sort = T, ncol = 3, assay = 'RNA', slot = "data") & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom")
  })
  output$Linear_RE <- renderPlot({
    sortie()
    result <- rownames(singlet@tools$avg.b.cells_RE[which(singlet@tools$avg.b.cells_RE$CTRL/singlet@tools$avg.b.cells_RE$STIM < 0.75),])
    ggplot(singlet@tools$avg.b.cells_RE, aes(CTRL, STIM)) + geom_point() + ggtitle("B Cells") 
    #Seurat::LabelPoints(plot = p1, points = result, repel = TRUE, xnudge = 0.2, ynudge = 0.5)
  })
  output$DE_DotPlot_RE <- renderPlot({
    sortie()
    Seurat::Idents(singlet)<-"Condition"
    Seurat::DotPlot(singlet, idents = c("Excipient","RCHOP"), features = rownames(singlet@tools$DE_RE)[1:12]) + Seurat::RotatedAxis()
  })
  output$DE_info_RE <- renderPrint({print(sortie()@tools$DE_RE)})
  
  # Pré-greffe / Excipient
  output$DE_Heatmap_PE <- renderPlot({
    sortie()
    Seurat::Idents(singlet)<-"Condition"
    Seurat::DoHeatmap(singlet, cells = rownames(singlet@meta.data)[which(singlet@meta.data$Condition==c("Excipient","Pré-greffe"))], features = rownames(singlet@tools$DE_PE)[1:50], size = 3, assay = 'SCT', slot = "scale.data")
  })
  output$DE_RidgePlot_PE <- renderPlot({
    sortie()
    Seurat::Idents(singlet)<-"Condition"
    Seurat::RidgePlot(singlet, idents = c("Excipient","Pré-greffe"), features = rownames(singlet@tools$DE_PE)[1:12], ncol = 3, assay = 'RNA', slot = "data") & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
  })
  output$DE_VlnPlot_PE <- renderPlot({
    sortie()
    Seurat::Idents(singlet)<-"Condition"
    Seurat::VlnPlot(singlet, idents = c("Excipient","Pré-greffe"), features = rownames(singlet@tools$DE_PE)[1:12], sort = T, ncol = 3, assay = 'RNA', slot = "data") & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom")
  })
  output$Linear_PE <- renderPlot({
    sortie()
    result <- rownames(singlet@tools$avg.b.cells_PE[which(singlet@tools$avg.b.cells_PE$CTRL/singlet@tools$avg.b.cells_PE$STIM > 2),])
    ggplot(singlet@tools$avg.b.cells_PE, aes(CTRL, STIM)) + geom_point() + ggtitle("B Cells") 
    #Seurat::LabelPoints(plot = p1, points = result, repel = TRUE, xnudge = 0.2, ynudge = 0.2)
  })
  output$DE_DotPlot_PE <- renderPlot({
    sortie()
    Seurat::Idents(singlet)<-"Condition"
    Seurat::DotPlot(singlet, idents = c("Excipient","Pré-greffe"), features = rownames(singlet@tools$DE_PE)[1:12]) + Seurat::RotatedAxis()
  })
  output$DE_info_PE <- renderPrint({print(sortie()@tools$DE_PE)})
  
  ## -- Enrichissement -- ##
  #output$Enrichissement <- renderPlot({sortie()@tools$enrichissement_RE})
  output$KEGG <- renderPlot({sortie()@tools$KEGG})
  output$GO_Biological <- renderPlot({sortie()@tools$GO_Biological})
  output$GO_Cellular <- renderPlot({sortie()@tools$GO_Cellular})
  output$GO_Molecular <- renderPlot({sortie()@tools$GO_Molecular})
  
  ## -- Heatmap -- ##
  output$Heatmap <- renderPlot({sortie() ; Seurat::DoHeatmap(singlet, features = heatmap()$gene, group.by = "Condition") + NoLegend()})
  output$Heatmap_feature <- renderPrint({print(sortie()@commands[["FindAllMarkers"]])})
  
  ## -- Mitochondrie Figure -- ##
  output$MT_VlnPlot <- renderPlot({sortie() ; singlet@tools$mitochondrie_all + Seurat::VlnPlot(singlet, features = "nFeature_RNA", group.by = "Condition") +  Seurat::VlnPlot(singlet, features = "nCount_RNA", group.by = "Condition") + Seurat::VlnPlot(singlet, features = "percent.mt", group.by = "Condition")})
  output$MT_FeatureScatter <- renderPlot({sortie() ; Seurat::FeatureScatter(singlet, feature1 = "nCount_RNA", feature2 = "percent.mt") + Seurat::FeatureScatter(singlet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")})
  output$MT_FeatureScatter2 <- renderPlot({sortie() ; Seurat::FeatureScatter(singlet, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "HTO_classification") + Seurat::FeatureScatter(singlet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "HTO_classification")})
  
  ## -- The 50 most highly variable genes -- ##
  output$top50 <- renderPlot({
    sortie()
    Seurat::VariableFeaturePlot(singlet)
    Seurat::LabelPoints(plot = Seurat::VariableFeaturePlot(singlet), points = head(VariableFeatures(singlet), 50), repel = TRUE)
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
    if(length(input$numSelector)==0){rv$hallmark <- sortie()@tools$hallmarks}
    else{rv$hallmark <- input$numSelector}
    rv$order <- input$hallmark_order
  })
  observeEvent(input$Subsets,{
    rv$hallmark <- sortie()@tools$hallmarks
    rv$order <-NULL}, ignoreNULL = FALSE)
  
  output$hallmark_Heatmap <- renderPlot({
    sortie()
    singlet@meta.data$active.idents <- singlet@active.ident
    singlet@tools$meta_variable <- c("seurat_clusters", "Condition", "Greffe", "Phénotype", "Phase", "orig.ident") # "clonotype_id",
    if(is.null(rv$order)){dittoSeq::dittoHeatmap(singlet, genes = NULL, metas = singlet@tools$hallmarks, heatmap.colors = rev(colorblind_vector(50)),annot.by = singlet@tools$meta_variable, cluster_cols = T, fontsize = 12)
    }else{dittoSeq::dittoHeatmap(singlet, genes = NULL, metas = rv$hallmark, heatmap.colors = rev(colorblind_vector(50)),annot.by = singlet@tools$meta_variable, cluster_cols = F, fontsize = 12, order.by = rv$order)}  
  })
  
  output$hallmark_VlnPlot <- renderPlot({dittoSeq::dittoPlot(sortie() , input$hallmark_order_vln, group.by = input$metadata_order_vln, legend.show = FALSE) + theme(title = element_text(size=20), axis.text = element_text(size=15)) +  ylab(label = "Score") })
  output$hallmark_HD <- renderPlot({dittoSeq::dittoScatterHex(sortie(),x.var = input$hallmark_order_X, y.var = input$hallmark_order_Y, do.contour = TRUE, split.by =  input$metadata_order_density) + theme_classic() + scale_fill_gradientn(colors = rev(colorblind_vector(11))) + geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2)})
  output$hallmark_RP <- renderPlot({
    sortie()
    ES2 <- data.frame(singlet[[]], Seurat::Idents(singlet))
    colnames(ES2)[ncol(ES2)] <- "cluster"
    escape::ridgeEnrichment(ES2, gene.set = input$hallmark_order_RP, group = input$metadata_group_RP, facet = input$metadata_facet_RP, add.rug = TRUE)
  })
  #output$hallmark_PCA <- renderPlot({
  #ES2 <- data.frame(singlet[[]], Seurat::Idents(singlet))
  #PCA <- performPCA(enriched = ES2, groups = c("cluster", "SingleR.calls"))
  #pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE)
  #pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = FALSE, facet = "cluster") 
  #masterPCAPlot(ES2, PCx = "PC1", PCy = "PC2", top.contribution = 10)
  #})
  
  ## -- Monocle -- ##
  output$cell_trajectory <- renderPlot({plot_cell_trajectory(data, color_by = input$color_trajectory)})
  output$order_trajectory <- renderPlot({
    data_expressed_genes <-  row.names(subset(fData(data), num_cells_expressed >= 10))
    data_filtered <- data[data_expressed_genes,]
    my_genes <- row.names(subset(fData(data_filtered), gene_short_name %in% input$in8))
    cds_subset <- data_filtered[my_genes,]
    p <- plot_genes_in_pseudotime(cds_subset, color_by = input$color_order)
    p + theme(text = element_text(size=20), legend.title = element_text(color = "black", size = 14),legend.text = element_text(color = "black", size = 14),axis.text=element_text(size=14),axis.title=element_text(size=14))
  })
  
  ## -- Sunburst -- ##
  output$dataTable = renderPlotly({plot_ly(sortie()@tools$sunburst, ids = ~ids ,labels = ~labels, parents = ~parents, values = ~values, marker = list(colors = c( "#BEBADA", "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462")), type = 'sunburst', branchvalues = 'total', hoverinfo = "text", hovertext = paste(singlet@tools$sunburst$labels, ":", round((singlet@tools$sunburst$values/singlet@tools$sunburst$values[1])*100,2),"%", "\nTotal : " ,singlet@tools$sunburst$values))})
  
  ## -- Clonotype -- ##
  output$VDJ_Clonotype = renderPlotly({
    plot_ly(x = sortie()@tools$Clonotype[[1]], y = sortie()@tools$Clonotype[[2]], name = "Info", type = "bar", 
            hovertemplate = paste0('Clonotype : %{x}\n', "Proportion : ", round(singlet@tools$vloupe$proportion[1:5],3)*100,"% \nType : ", singlet@tools$vloupe$type[1:5], "\nIsotype : ", singlet@tools$vloupe$igh_c_genes[1:5], 
                                   "\nHeavy : ", singlet@tools$vloupe$V_lourde[1:5], " / ", singlet@tools$vloupe$D_lourde[1:5], " / ", singlet@tools$vloupe$J_lourde[1:5], "\nLight : ", singlet@tools$vloupe$V_legere[1:5], " / ", singlet@tools$vloupe$J_legere[1:5])
    ) %>% layout(title='Frequencies of the mains clonotypes', yaxis =list(title="Number of cells"))
  })
  
  ## -- VDJ -- ##
  output$V = renderPlotly({plot_ly(x = sortie()@tools$V[[1]], y = sortie()@tools$V[[2]], name = "Clonotype", type = "bar") %>% layout(title='Frequencies V genes : Heavy and Lights chains', xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))})
  output$D = renderPlotly({plot_ly(x = sortie()@tools$D[[1]], y = sortie()@tools$D[[2]], name = "Clonotype", type = "bar") %>% layout(title='Frequencies D genes : Heavy chain',xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))})
  output$J = renderPlotly({plot_ly(x = sortie()@tools$J[[1]], y = sortie()@tools$J[[2]], name = "Clonotype", type = "bar") %>% layout(title='Frequencies J genes : Heavy and Lights chains',xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))})
  
  ## -- Heavy chain -- ##
  output$VDJ_Heavy = renderPlotly({plot_ly(x = sortie()@tools$Heavy[[1]], y = sortie()@tools$Heavy[[2]], name = "Heavy Chain", type = "bar", hovertemplate = paste0("Locus : %{x}\nProportion : ", round(singlet@tools$Heavy[[2]]/sum(singlet@tools$Heavy[[2]]),3)*100,"%")) %>% layout(title='Frequencies Heavy Chain' , xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))})
  output$VDJ_DHeavy = renderPlotly({plot_ly(x = sortie()@tools$Isotype[[1]], y = sortie()@tools$Isotype[[2]], name = "Heavy Chain", type = "bar", hovertemplate = paste0("Locus : %{x}\nProportion : ", round(singlet@tools$Isotype[[2]]/sum(singlet@tools$Isotype[[2]]),3)*100,"%")) %>%layout(title='Frequencies Heavy Chain with details' , xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))  })
  
  ## -- Light chain -- ##
  output$VDJ_Light = renderPlotly({plot_ly(x = sortie()@tools$Light[[1]], y = sortie()@tools$Light[[2]], name = "Light Chain", type = "bar", hovertemplate = paste0("Locus : %{x}\nProportion : ", round(singlet@tools$Light[[2]]/sum(singlet@tools$Light[[2]]),3)*100,"%")) %>% layout(title='Frequencies Light Chain',yaxis =list(title="Number of cells"))})
  output$VDJ_DLight = renderPlotly({plot_ly(x = sortie()@tools$Type[[1]], y = sortie()@tools$Type[[2]], name = "Clonotype", type = "bar", hovertemplate = paste0("Locus : %{x}\nProportion : ", round(singlet@tools$Type[[2]]/sum(singlet@tools$Type[[2]]),3)*100,"%")) %>% layout(title='Frequencies Light Chain with details', yaxis =list(title="Number of cells"))})
  
  ## -- Méta-VDJ -- ##
  output$FL08_VDJ = renderPlotly({load(file = "www/FL08G0293/FL08G0293_VDJ.RData") ; e})
  output$FL09_VDJ = renderPlotly({load(file = "www/FL09C1164/FL09C1164_VDJ.RData") ; e})
  output$FL12_VDJ = renderPlotly({load(file = "www/FL12C1888/FL12C1888_VDJ.RData") ; e})
  output$FL14_VDJ = renderPlotly({load(file = "www/FL140304/FL140304_VDJ.RData") ; e})
  output$FL02_VDJ = renderPlotly({load(file = "www/FL02G095/FL02G095_VDJ.RData") ; e})
  output$FL05_VDJ = renderPlotly({load(file = "www/FL05G0330/FL05G0330_VDJ.RData") ; e})
  
  output$FL08_Meta = renderPlotly({load(file = "www/FL08G0293/FL08G0293_VDJ.RData") ; f})
  output$FL09_Meta = renderPlotly({load(file = "www/FL09C1164/FL09C1164_VDJ.RData") ; f})
  output$FL12_Meta = renderPlotly({load(file = "www/FL12C1888/FL12C1888_VDJ.RData") ; f})
  output$FL14_Meta = renderPlotly({load(file = "www/FL140304/FL140304_VDJ.RData") ; f})
  output$FL02_Meta = renderPlotly({load(file = "www/FL02G095/FL02G095_VDJ.RData") ; f})
  output$FL05_Meta = renderPlotly({load(file = "www/FL05G0330/FL05G0330_VDJ.RData") ; f})
  
  
  ## -- Rapport Patient -- ##
  output$PDF_report <- downloadHandler(
    filename = function() {paste0("rapport_", input$patient,".pdf")},
    content = function(file) {file.copy(paste0("www/",input$patient,"/rapport_", input$patient,".pdf"), file)}
  )
  output$HTML_report <- downloadHandler(
    filename = function() {paste0("rapport_", input$patient,".html")},
    content = function(file) {file.copy(paste0("www/",input$patient,"/rapport_", input$patient,".html"), file)}
  )
  output$Word_report <- downloadHandler(
    filename = function() {paste0("rapport_", input$patient,".docx")},
    content = function(file) {file.copy(paste0("www/",input$patient,"/rapport_", input$patient,".docx"), file)}
  )
})