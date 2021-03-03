######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## LIBRARY                                                                                                                       ########
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
library("tidyverse") # for data manipulation & plots
library("lubridate") # right color
library("BUSpaRse")
library("TENxBUSData")
library("simpleSingleCell") # for scRNA storage & manipulation
library("scater") # for QC control
library("scran") # analysis pipeline
library("Seurat")
library("uwot") # UMAP dim-red
library("DropletUtils") #utility functions for handling single-cell (RNA-seq)
library("AnnotationHub") # ensbl query
library("AnnotationDbi") # ensbl query
library("sctransform") # sc normalization,
library("SingleR")
library("velocyto.R")
library(Matrix)
library(stringr)
library(celldex)
library(scRNAseq)
library(CellSIUS)
require(RColorBrewer)
library(ggplot2)
library(rowr)
library(shiny)
library(shinyjs)
library(grid)
library(gridExtra)
library(shinydashboard)
library(Signac)
library(ape)
library(plotly)
library(lemon)
library(scone)
suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(dittoSeq))

setwd(dir = "/home/boris/Documents/lipinskib/boris/Cellranger/result/")

######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## CREATION DE L'OBJECT SEURAT + DEMULTIPLEXAGE                                                                                  ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ########
patient <- "hFL_180008B"

seurat_object <- function(patient){
  ############################################################################################################################################# 
  ##### --        Préprocessing         -- #####
  ##############################################
  # Matrice mRNA + object seurat
  rwa.mRNA <- Read10X(data.dir = paste0(patient,"/mRNA/",patient,"_CellrangerCount/outs/filtered_feature_bc_matrix/")) 
  cso <- CreateSeuratObject(counts = rwa.mRNA, project = patient, min.cells = 3, min.features = 200)
  umis <- GetAssayData(object = cso, slot = "counts")
  
  # Matrice HTO
  raw.hto <- Read10X(paste0(patient,"/HTO/umi_count/"), gene.column = 1)        
  colnames(raw.hto) <- paste0(colnames(raw.hto),"-1")                          
  hto <- raw.hto[c(1:4),]                                                       # Suppression des séquences unmapped : rownames(raw.hto)
  rownames(hto) <- c("pregreffe","subtilisin","collagenase","actinomycin_D")    
  
  joint.bcs <- intersect(colnames(umis),colnames(hto))                          # Sélection des cellules avec barcode commun HTO / mRNA
  
  umis <- umis[, joint.bcs]                                                     # Sélection des lignes qui correspondent aux cellules en commun
  hto <- as.matrix(hto[, joint.bcs])
  
  # Object seurat : HTOs + mRNA
  hashtag <- CreateSeuratObject(counts = umis, assay = "RNA", project = patient)
  hashtag <- NormalizeData(hashtag)
  hashtag <- FindVariableFeatures(hashtag, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  hashtag <- ScaleData(hashtag,features = VariableFeatures(hashtag))
  
  hashtag[["HTO"]] <- CreateAssayObject(counts = hto)                           # Ajoute des données HTO comme un nouvel assay indépendant du mRNA
  hashtag <- NormalizeData(hashtag, assay = "HTO",normalization.method = "CLR")
  
  # Association : cellules / échantillons
  hashtag <- HTODemux(hashtag, assay = "HTO", positive.quantile = 0.99)
  table(hashtag$HTO_classification.global)                                      # Result
  
  Idents(hashtag)<- 'HTO_classification.global'
  singlet <- subset(hashtag, idents = "Singlet")
  
  
  
  
  #############################################################################################################################################
  ##### -- identification des cellules -- ##### 
  #############################################
  # DB: DatabaseImmuneCellExpressionData()
  hpca.se <- celldex::DatabaseImmuneCellExpressionData()
  results <- SingleR(test = as.SingleCellExperiment(singlet), ref = hpca.se, labels = hpca.se$label.main)
  results.fine <- SingleR(test = as.SingleCellExperiment(singlet), ref = hpca.se, labels = hpca.se$label.fine)
  
  #singlet$SingleR.pruned.calls <- results$pruned.labels
  singlet$SingleR.calls <- results$labels
  #singlet$SingleR.pruned.calls.fine <- results.fine$pruned.labels
  singlet$SingleR.calls.fine <- results.fine$labels
  
  # DB: BlueprintEncodeData()
  #BD <- celldex::BlueprintEncodeData()
  #results.blueprint <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.main)
  #results.blueprint.fine <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.fine)
  
  #singlet$SingleR.pruned.calls.blueprint <- results$pruned.labels
  #singlet$SingleR.calls.blueprint <- results$labels
  #singlet$SingleR.pruned.calls.blueprint.fine <- results.blueprint.fine$pruned.labels
  #singlet$SingleR.calls.blueprint.fine <- results.blueprint.fine$labels
  
  # Récap
  table(results$labels)
  table(results.fine$labels)
  
  
  #############################################################################################################################################
  ##### -- Enrichissement des gènes -- ##### 
  ##########################################
  GS <- getGeneSets(library = "H")
  ES <- enrichIt(obj = singlet, gene.sets = GS, groups = 1000, cores = 12)
  singlet <- AddMetaData(singlet, ES)
  singlet@tools$hallmarks <- names(ES)
  singlet@tools$meta_variable <- c("seurat_clusters", "HTO_maxID", "SingleR.calls", "clonotype_id", "Phase")
  
  #############################################################################################################################################
  ##### --             VDJ             -- ##### 
  #############################################
  # BCR
  bcr <- read.csv(paste0(patient,"/VDJ/",patient,"_CellrangerVDJ/outs/filtered_contig_annotations.csv")) #bcr$barcode <- gsub("-1", "", bcr$barcode)                                   # Remove the -1 at the end of each barcode: VERIFIER SI IL Y A LES TIRETS DANS OBJECT SINGLET!!!
  bcr <- bcr[!duplicated(bcr$barcode), ]                                        # Subsets so only the first line of each barcode is kept,as each entry for given barcode will have same clonotype.
  bcr <- bcr[,c("barcode", "raw_clonotype_id","chain","v_gene","d_gene","j_gene","c_gene","cdr3")]
  names(bcr)[names(bcr) == "raw_clonotype_id"] <- "clonotype_id"
  
  # Clonotypes
  clonotype <- read.csv(paste0(patient,"/VDJ/",patient,"_CellrangerVDJ/outs/clonotypes.csv"))
  
  # Merge
  bcr <- merge(bcr, clonotype[, c("clonotype_id", "cdr3s_aa")])                 # Selection de ceux pour lesquels on a une séquence cdr3
  bcr <- bcr[, c(2,1,3,4,5,6,7,8,9)]                                            # Mettre les codes barres comme nom de la première colonne
  rownames(bcr) <- bcr[,1]
  bcr[,1] <- NULL
  
  # Selection des 4 premiers clonotypes : optionnel
  raw <- rbind(bcr[bcr$clonotype_id == 'clonotype1',],
               bcr[bcr$clonotype_id == 'clonotype2',],
               bcr[bcr$clonotype_id == 'clonotype3',],
               bcr[bcr$clonotype_id == 'clonotype4',])
  
  bcr$clonotype_id <- as.numeric(str_replace(unlist(bcr$clonotype_id,"clonotype",use.names=F), "clonotype", ""))
  
  # Add BCR to metadata
  singlet <- AddMetaData(object=singlet, metadata = bcr)                        #, col.name = colnames(as.data.frame(bcr)))   # Add to the Seurat object's metadata.
  
  return(singlet)
}
singlet <- seurat_object(patient)

#### Extraction des sous pop HTO : pregreffe / actinomycin RCHOP / actinomycin placébo 
Idents(singlet)<-"hash.ID"
pregreffe <- subset(singlet, idents = "pregreffe")
actinomycin <- subset(singlet, idents = "actinomycin-D")
collagenase <- subset(singlet, idents = "collagenase")
subtilisin <- subset(singlet, idents = "subtilisin")
postgreffe <- subset(singlet, idents = c("subtilisin", "actinomycin-D", "collagenase"))




######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## Visualisation                                                                                                                 ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
visualitation <- function(singlet){
  annotations <- read.csv("/home/boris/Bureau/R_project/scShiny/annotation_FindAllMarkers.csv")
  
  ## -- ADD MITOCHONDRIAL ANALYSES -- ## 
  singlet[["percent.mt"]] <- PercentageFeatureSet(singlet, pattern = "^MT-")
  singlet <- subset(singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)            # QC Filter : tester 15%
  
  ## -- ADD PCA TO SEURAT OBJECT -- ## 
  # Pre-processing  
  singlet <- NormalizeData(singlet, normalization.method = "LogNormalize", scale.factor = 10000) 
  singlet <- FindVariableFeatures(singlet, selection.method = "vst", nfeatures = 2000)
  
  singlet <- ScaleData(singlet, features = rownames(singlet))                                               # Scale data :singlet@assays$RNA@scale.data, singlet[["RNA"]]@scale.data
  
  # Add phase cycle information 
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  singlet <- CellCycleScoring(singlet, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)    # Attribution à chaque cellule d'un état de phase
  #singlet <- ScaleData(singlet, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(singlet)) # Correction de la matrice d'expression dépendemment du cycle cellulaire
  
  # Reduction dimension
  singlet <- RunPCA(singlet, features = VariableFeatures(singlet), ndims.print = 1:10, nfeatures.print = 30)
  singlet <- FindNeighbors(singlet, reduction = "pca", dims = 1:40)
  singlet <- FindClusters(singlet, resolution = 0.5)                                                        # head(Idents(singlet), 10) 
  singlet <- RunUMAP(singlet, reduction = "pca", dims = 1:40)
  singlet <- RunTSNE(singlet, reduction = "pca", dims = 1:40)
  
  ## -- Ohter PCA analysis -- ## 
  singlet <- JackStraw(singlet, num.replicate = 100, reduction = "pca")
  singlet <- ScoreJackStraw(singlet, dims = 1:20)
  
  singlet@commands[["FindAllMarkers"]] <- FindAllMarkers(singlet, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  #singlet@commands[["FindAllMarkers"]] <- merge(singlet@commands[["FindAllMarkers"]], annotations, by.x="gene", by.y="gene_name")
  #singlet@commands[["FindAllMarkers"]] %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  
  singlet@meta.data$nFeature_HTO <- NULL
  singlet@meta.data$HTO_secondID <- NULL
  singlet@meta.data$HTO_classification <- NULL
  singlet@meta.data$RNA_snn_res.0.5 <- NULL
  singlet@meta.data$hash.ID <- NULL
  
  return(singlet)
}
singlet <- visualitation(singlet)

save(singlet, file = paste0("/home/boris/Documents/analyse/singlet_",patient,".RData"))
write.csv(singlet@meta.data, file = paste0("/home/boris/Documents/analyse/jupyter/metadata_matrix_",patient,".csv"))
write.csv(as.matrix(singlet[["RNA"]]@counts), file = paste0("/home/boris/Documents/analyse/jupyter/count_matrix_", patient,".csv"))

load(file = paste0("/home/boris/Documents/analyse/singlet_",patient,".RData"))

#save(singlet, file = paste0(patient,"/R/singlet_",patient,".RData"))
#write.csv(singlet@meta.data, file = paste0(patient,"/R/metadata_matrix_",patient,".csv"))
#write.csv(as.matrix(singlet[["RNA"]]@counts), file = paste0(patient,"/R/count_matrix_", patient,".csv"))


######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## ESCAPE - Enrichissement de gène                                                                                               ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

#The Heatmap
singlet@meta.data$active.idents <- singlet@active.ident
dittoHeatmap(singlet, genes = NULL, metas = singlet@tools$hallmarks, heatmap.colors = rev(colorblind_vector(50)),
             annot.by = singlet@tools$meta_variable, cluster_cols = TRUE, fontsize = 7)

dittoHeatmap(singlet, genes = NULL, metas = names(ES), heatmap.colors = rev(colorblind_vector(50)),
             annot.by = meta_variable, cluster_cols = TRUE, fontsize = 7)

dittoHeatmap(singlet, genes = NULL, metas = c("HALLMARK_APOPTOSIS", "HALLMARK_DNA_REPAIR", "HALLMARK_P53_PATHWAY"), 
             heatmap.colors = rev(colorblind_vector(50)), annot.by = meta_variable, cluster_cols = TRUE, fontsize = 7)

#The Violin Plot
dittoPlot(singlet, "HALLMARK_DNA_REPAIR", group.by = "SingleR.calls") + scale_fill_manual(values = colorblind_vector(5))

#Hex Density Enrichment Plots
dittoScatterHex(singlet,x.var = "HALLMARK_DNA_REPAIR", y.var = "HALLMARK_MTORC1_SIGNALING", do.contour = TRUE) + theme_classic() + 
  scale_fill_gradientn(colors = rev(colorblind_vector(11))) + geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2)  

dittoScatterHex(singlet, x.var = "HALLMARK_DNA_REPAIR", y.var = "HALLMARK_MTORC1_SIGNALING", do.contour = TRUE, split.by = "SingleR.calls") + 
  theme_classic() + scale_fill_gradientn(colors = rev(colorblind_vector(11))) + geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2) 

# Enrichment along a Ridge Plot
ES2 <- data.frame(singlet[[]], Idents(singlet))
colnames(ES2)[ncol(ES2)] <- "cluster"
ridgeEnrichment(ES2, gene.set = "HALLMARK_DNA_REPAIR", group = "SingleR.calls", add.rug = TRUE)
ridgeEnrichment(ES2, gene.set = "HALLMARK_DNA_REPAIR", group = "cluster", facet = "SingleR.calls", add.rug = TRUE)

# The Split Violin Plot
splitEnrichment(ES2, split = "SingleR.calls", gene.set = "HALLMARK_DNA_REPAIR")
splitEnrichment(ES2, x.axis = "cluster", split = "SingleR.calls", gene.set = "HALLMARK_DNA_REPAIR")

# Expanded Analysis
ES2 <- data.frame(singlet[[]], Idents(singlet))
PCA <- performPCA(enriched = ES2, groups = c("cluster", "SingleR.calls"))
pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE)
pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = FALSE, facet = "cluster") 
masterPCAPlot(ES2, PCx = "PC1", PCy = "PC2", top.contribution = 10)

#Signficance
output <- getSignificance(ES2, group = "cluster", fit = "linear.model")

######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## To test                                                                                                                       ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#annotations_finale <- singlet@meta.data[which(singlet@meta.data$gene_name %in% top10),] 
#singlet <- BuildClusterTree(singlet)
#Tool(singlet,'BuildClusterTree')
#PlotClusterTree(singlet)

#singlet <- SCTransform(singlet)
#Stdev(singlet[["pca"]])

sce.singlet <- as.SingleCellExperiment(singlet)
colData(sce.singlet)

sce.singlet@metadata <- singlet@meta.data
metadata(sce.singlet)$which_qc

fluidigm <- as.SingleCellExperiment(singlet)
#Add QC for Mitochondrial genes
fluidigm <- addPerCellQC(fluidigm)
colData(fluidigm)
fluidigm <- addPerFeatureQC(fluidigm)
rowData(fluidigm)

View(fluidigm)


######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## 3D Visualisation                                                                                                              ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
# cell plotting : Embeddings(object = singlet, reduction = "umap") --> Visualize what headings are called so that you can extract them to form a dataframe
singlet2 <- singlet
singlet2 <- RunUMAP(singlet2, reduction = "pca", dims = 1:40, n.components = 3L)
singlet2 <- RunTSNE(singlet2, reduction = "pca", dims = 1:40, dim.embed = 3)
meta_variable = c("orig.ident", "HTO_maxID", "SingleR.calls", "clonotype_id","chain","v_gene", "d_gene", "j_gene","c_gene", "cdr3", "Phase", "old.ident", "seurat_clusters")
plot.data <- FetchData(object = singlet2, vars = c(meta_variable, "PC_1", "PC_2", "PC_3", "tSNE_1", "tSNE_2", "tSNE_3", "UMAP_1", "UMAP_2", "UMAP_3"))
plot.data$label <- paste(rownames(plot.data))

plot_ly(data = plot.data, x = ~PC_1, y = ~PC_2, z = ~PC_3, 
        color = ~seurat_clusters, colors = RColorBrewer::brewer.pal(length(levels(singlet@meta.data$seurat_clusters)),"Spectral"),
        type = "scatter3d", mode = "markers", 
        marker = list(size = 3, width=2),
        text=~label, hoverinfo="text")

plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
        color = ~seurat_clusters, colors = RColorBrewer::brewer.pal(length(levels(singlet@meta.data$seurat_clusters)),"Spectral"),
        type = "scatter3d", mode = "markers", 
        marker = list(size = 3, width=2),
        text=~label, hoverinfo="text")

plot_ly(data = plot.data, x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
        color = ~seurat_clusters, colors = RColorBrewer::brewer.pal(length(levels(singlet@meta.data$seurat_clusters)),"Spectral"),
        type = "scatter3d", mode = "markers", 
        marker = list(size = 3, width=2),
        text=~label, hoverinfo="text")


gene_feature <- "CD19"
all_feature <- rownames(as.matrix(singlet2[["RNA"]]@counts))
plot.data <- FetchData(object = singlet2, vars = c(all_feature, "PC_1", "PC_2", "PC_3", "tSNE_1", "tSNE_2", "tSNE_3", "UMAP_1", "UMAP_2", "UMAP_3"), slot = 'data')
plot.data$changed <- ifelse(test = plot.data[[gene_feature]] <1, yes = plot.data[[gene_feature]], no = 1)
plot.data$label <- paste(rownames(plot.data)," - ", plot.data[[gene_feature]], sep="")

plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        colors = c('darkgreen', 'red'), opacity = .5,
        type = "scatter3d", mode = "markers",
        marker = list(size = 3, width=2), 
        text=~label, hoverinfo="text"
) %>% layout(title=gene_feature)

plot_ly(data = plot.data, x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
        color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        colors = c('darkgreen', 'red'), opacity = .5,
        type = "scatter3d", mode = "markers",
        marker = list(size = 3, width=2), 
        text=~label, hoverinfo="text"
) %>% layout(title=gene_feature)

plot_ly(data = plot.data, x = ~PC_1, y = ~PC_2, z = ~PC_3, 
        color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        colors = c('darkgreen', 'red'), opacity = .5,
        type = "scatter3d", mode = "markers",
        marker = list(size = 3, width=2), 
        text=~label, hoverinfo="text"
) %>% layout(title=gene_feature)


Cutoff <- quantile(plot.data[,gene_feature], probs = .95)
plot.data$"ExprCutoff" <- ifelse(test = plot.data[,gene_feature] < Cutoff, yes = plot.data[,gene_feature], no = Cutoff)
plot.data$label <- paste(rownames(plot.data)," - ", plot.data[,gene_feature], sep="")

plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, name = gene_feature, # Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
        color = ~ExprCutoff, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        colors = c('darkgrey', 'red'), opacity = .5,
        type = "scatter3d", mode = "markers", marker = list(size = 1), 
        text=~label,hoverinfo="text"
) %>% layout(title=gene_feature)




######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## Summary metadata                                                                                                              ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
result=list()
summary <- list(singlet@meta.data)
names(summary) <- patient
resume=list()
for (elm in 1:length(names(summary[[patient]]))) {
  if (names(summary[[patient]])[elm] %in% meta_variable) {
    resume[[names(summary[[1]])[elm]]] = as.list(as.vector(table(summary[[1]][[elm]])))
    names(resume[[names(summary[[1]])[elm]]]) = rownames(as.matrix(table(summary[[1]][[elm]])))
  }
}
result[[patient]] <- resume
save(result, file = "/home/boris/Documents/analyse/shiny/result.RData")
load(file = "/home/boris/Documents/analyse/shiny/result.RData")




######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## Analyses                                                                                                                      ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
# Main visualisation : TSNE single and multi label
DimPlot(singlet, reduction = "tsne", label = TRUE)
BD <- "SingleR.calls" #"SingleR.calls.fine", "SingleR.calls.blueprint", "SingleR.calls.blueprint.fine"


Idents(singlet)<-"HTO_maxID"
tokeep <- levels(Idents(singlet))
tokeep <- tokeep[tokeep %in% "pregreffe"]
singlet2 <- subset(singlet, idents = tokeep)
Idents(singlet2)<-"seurat_clusters"


plots <- TSNEPlot(object = singlet2, group.by = c("HTO_maxID", BD), split.by = "Phase", label.size = 0.0, pt.size = 1)
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 1)))

# Heatmap de correlation
top10 <- singlet@commands[["FindAllMarkers"]] %>% group_by("HTO_maxID") %>% top_n(n = 100, wt = avg_log2FC)
DoHeatmap(singlet, features = top10$gene, group.by = "HTO_maxID") + NoLegend()

# PCA
DimPlot(singlet, reduction = "pca")
print(singlet[["pca"]], dims = 1:10, nfeatures = 10)
ElbowPlot(singlet, ndims = 50, reduction = "pca")
DimHeatmap(singlet, dims = 1:10, cells = 100, balanced = TRUE)
VizDimLoadings(singlet, dims = 1:5, reduction = "pca")
JackStrawPlot(singlet, dims = 1:15)

# UMAP
DimPlot(singlet, reduction = "umap", label = TRUE)
DimPlot(singlet, reduction = "umap", group.by = "HTO_maxID")

# TSNE
DimPlot(singlet, reduction = "tsne", group.by = "HTO_maxID")
DimPlot(singlet, reduction = "tsne", group.by = BD)
DimPlot(singlet, reduction = "tsne", group.by = "Phase")
TSNEPlot(object = singlet, group.by = "clonotype_id", label.size = 0.0, pt.size = 1)

# Mitochondrie 
VlnPlot(singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "HTO_classification")
VlnPlot(singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # Visualize QC metrics as a violin plot
FeatureScatter(singlet, feature1 = "nCount_RNA", feature2 = "percent.mt") + FeatureScatter(singlet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 <- FeatureScatter(singlet, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "HTO_classification")
plot2 <- FeatureScatter(singlet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "HTO_classification")
plot1 + plot2

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(singlet), 10)
VariableFeaturePlot(singlet)
LabelPoints(plot = VariableFeaturePlot(singlet), points = top10, repel = TRUE) # plot variable features with and without labels

annotations_finale <- annotations[which(annotations$gene_name %in% top10),] 

# Visualize the distribution of cell cycle markers across
RidgePlot(singlet, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Verification Expression marqueurs B_T
feature <-c("MS4A1", "CD19", "CD79A", "CD79B")
FeaturePlot(singlet, features = feature, cols = RColorBrewer::brewer.pal(6,"YlOrBr"), combine = TRUE, blend = TRUE) & NoAxes() & NoLegend()
FeaturePlot(singlet, features = feature, split.by = "Phase")
FeaturePlot(singlet, features = c("MS4A1", "CD19", "CD79A", "CD79B"), reduction='umap')
FeaturePlot(singlet, features = c("MS4A1", "CD19", "CD79A", "CD79B"), reduction='tsne')
VlnPlot(singlet, features = c("MS4A1", "CD19", "CD79A", "CD79B"), group.by = "HTO_classification")
FeaturePlot(singlet, features = c("CD3E", "CD3D", "CD56", "CD4", "CD2", "CD25", "CD62L", "CD197")) # Expression marqueurs T cells
VlnPlot(singlet, features = c("CD3E", "CD3D", "CD56", "CD4", "CD2", "CD25", "CD62L", "CD197"), group.by = "HTO_classification")




######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## VDJ - Métadata                                                                                                                ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
## ATTENTION pour que installation scRepertoire ok il faut d'abord désintaller ggplot2 et Seurat, installer scRepertoire puis re installer ggplot2 et Seurat
suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))
library(scRepertoire)

# Chargement des données VDJ
patient1 <- "hFL_180008B"
csv1 <- list(read.csv(paste0(patient1,"/VDJ/",patient1,"_CellrangerVDJ/outs/filtered_contig_annotations.csv"), stringsAsFactors = FALSE))
patient2 <- "hFL_130337"
csv2 <- list(read.csv(paste0(patient2,"/VDJ/",patient2,"_CellrangerVDJ/outs/filtered_contig_annotations.csv"), stringsAsFactors = FALSE))

combined <- combineBCR(c(csv1,csv2), samples = c(patient1, patient2), ID = c("BCR","BCR"))

quantContig(combined, cloneCall="gene+nt", scale = TRUE)
quantContig_output <- quantContig(combined, cloneCall="gene+nt", scale = TRUE, exportTable = TRUE)
quantContig_output




######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## CellSIUS - Cell Subtype Identification from Upregulated gene Sets                                                             ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
singlet_matrix <- as.matrix(singlet[["RNA"]]@counts)
CellSIUS.out <- CellSIUS(singlet_matrix, as.character(Idents(singlet)))
Result_List = CellSIUS_GetResults(CellSIUS.out=CellSIUS.out)

my_tsne <- Embeddings(singlet[["tsne"]])
d = CellSIUS_plot(coord = my_tsne, CellSIUS.out = CellSIUS.out)
d + scale_color_brewer(palette = "Blues")

Final_Clusters = CellSIUS_final_cluster_assignment(CellSIUS.out=CellSIUS.out, group_id=Idents(singlet), min_n_genes = 3)
table(Final_Clusters)
p <- qplot(my_tsne[,1],my_tsne[,2],xlab="tSNE1",ylab="tSNE2",color=as.factor(Final_Clusters))
p + theme(legend.position = "none")




######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## GENOTYPAGE - Long Ranger                                                                                                      ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#vawk '{print $1,$2}' variants.vcf > SNV.loci.txt
#sed -i 's/\s/:/g' SNV.loci.txt 
setwd(dir = "/home/boris/Documents/lipinskib/")

# Construct the data.frame
barcodes <- read.table("/home/boris/Documents/lipinskib/flinovo/datas/datas_analysees/hFL_180008B/hFL_180008B_CellrangerCount/outs/filtered_feature_bc_matrix/barcodes.tsv", header = F) #read in the cell barcodes output by Cell Ranger
snps <- read.table("/home/boris/Documents/lipinskib/boris/longranger/result/wgs/hFL_180008B/outs/SNV.loci.txt", header = F) # Construct the final table to add to the Seurat object

# charge raw data
snv_matrix <- readMM("/home/boris/Documents/lipinskib/boris/longranger/result/vartex/hFL_180008B/matrix.mtx") # Read in the sparse genotype matrix
colnames(snv_matrix) <- barcodes$V1
row.names(snv_matrix) <- snps$V1 #View(head(snv_matrix[1:10]))

snv_matrix <- as.data.frame(as.matrix(t(snv_matrix))) # convert the matrix to a dataframe #View(head(snv_matrix[1:10]))
#sub_snv_matrix <- snv_matrix[98640:98825] #BCL2 target

colnames(snv_matrix)
taget_name <-"chr18:63319316" #or colnames(snv_matrix)[1]
target <- data.frame(snv_matrix[taget_name])

row.names(target) <- barcodes$V1 ; colnames(target) <- taget_name

target[taget_name] <- str_replace(as.character(target[taget_name]), "0", "No Call") # No reads detected
target[taget_name] <- str_replace(as.character(target[taget_name]), "1", "ref/ref") # Only ref detected
target[taget_name] <- str_replace(as.character(target[taget_name]), "2", "alt/alt") # Only alt detected
target[taget_name] <- str_replace(as.character(target[taget_name]), "3", "alt/ref") # Both alleles detected

target <- t(target)
target <- target[ , which(colnames(target) %in% colnames(singlet))] 

singlet <- AddMetaData(object = singlet, metadata = as.data.frame(target), col.name = colnames(as.data.frame(target)))
TSNEPlot(object = singlet, group.by = "target", label.size = 0.0, pt.size = 1)
#do.label = T, colors.use = c("azure2","black", "yellow","red"), plot.title = "chr18_63319316", do.return = T, plot.order = c("alt/alt" ,"alt/ref", "ref/ref","No Call"))
#tableau <- data.frame(c(1:638),c(1:638),c(1:638),c(1:638)); for (i in 1:length(snv_matrix)){tableau[i,] <- table(factor(as.vector(snv_matrix[i]), levels = c(0:3)))} ;colnames(tableau) <- c('0','1','2','3')