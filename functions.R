#.libPaths( c( .libPaths(), '/home/boris/R/x86_64-pc-linux-gnu-library/4.1/') )
library(Seurat)
library(dittoSeq)
library(plotly)
library(shiny)
library(shinybusy)
library(shinyWidgets)
library(shinydashboard)
library(dashboardthemes)
library(shinyjs)

## -- Worflow -- ##
metadata <- function(singlet){  
  ##### -- Identification des cellules -- #####
  BD <- celldex::BlueprintEncodeData()
  results.blueprint <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.main)
  results.blueprint.fine <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.fine)
  singlet$Phénotype <- results.blueprint$labels 
  singlet$Phénotype.fine <- results.blueprint.fine$labels 
  
  
  ##### -- Enrichissement des gènes -- ##### 
  GS <- getGeneSets(library = "H")
  ES <- enrichIt(obj = singlet, gene.sets = GS, groups = 1000, cores = 12)
  names(ES) <- str_replace_all(names(ES), "HALLMARK_", "")
  singlet <- AddMetaData(singlet, ES)
  singlet@tools$hallmarks <- names(ES)
  
  
  ##### -- GoT  -- ##### 
  load(file=paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient, "/R/", patient, "_GoT.Rdata"))
  singlet <- AddMetaData(object = singlet, metadata = GOT)
  
  
  ##### -- Phase cycle -- #### 
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  singlet <- CellCycleScoring(singlet, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)    # Attribution à chaque cellule d'un état de phase
  #singlet <- ScaleData(singlet, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(singlet)) # Correction de la matrice d'expression dépendemment du cycle cellulaire
  Idents(singlet)<-"seurat_clusters"
  
  
  ##### -- VDJ -- ##### 
  # Cellranger
  bcr <- read.csv(paste0(patient,"/VDJ/CellrangerVDJ/outs/filtered_contig_annotations.csv"))  # bcr$barcode <- gsub("-1", "", bcr$barcode)                                   # Remove the -1 at the end of each barcode: VERIFIER SI IL Y A LES TIRETS DANS OBJECT SINGLET!!!
  bcr <- bcr[!duplicated(bcr$barcode), ]                                                                  # Subsets so only the first line of each barcode is kept,as each entry for given barcode will have same clonotype.
  bcr <- bcr[,c("barcode", "raw_clonotype_id","chain","v_gene","d_gene","j_gene","c_gene","cdr3")] ; names(bcr)[names(bcr) == "raw_clonotype_id"] <- "clonotype_id"
  clonotype <- read.csv(paste0(patient,"/VDJ/CellrangerVDJ/outs/clonotypes.csv"))             # Clonotypes
  bcr <- merge(bcr, clonotype[, c("clonotype_id", "cdr3s_aa")])                 # Selection de ceux pour lesquels on a une séquence cdr3
  bcr <- bcr[, c(2,1,3,4,5,6,7,8,9)]                                            # Mettre les codes barres comme nom de la première colonne
  rownames(bcr) <- bcr[,1]
  bcr[,1] <- NULL
  bcr$clonotype_id <- as.numeric(str_replace(unlist(bcr$clonotype_id,"clonotype",use.names=F), "clonotype", ""))
  singlet <- AddMetaData(object=singlet, metadata = bcr)                        #, col.name = colnames(as.data.frame(bcr)))   # Add to the Seurat object's metadata.
  
  # Figure
  load(file=paste0("/home/boris/Documents/lipinskib/flinovo/result/",patient,"/R/VDJ.RData"))
  singlet@tools$Clonotype <- Clonotype
  singlet@tools$Type <- Type
  singlet@tools$Isotype <- Isotype
  singlet@tools$V <- V
  singlet@tools$D <- D
  singlet@tools$J <- J
  singlet@tools$Heavy <- Heavy
  singlet@tools$Light <- Light
  singlet@tools$vloupe <- vloupe
  
  
  ##### -- Ohter -- #####
  singlet@tools$meta_variable <- c("seurat_clusters", "Condition", "Greffe", "Phénotype", "orig.ident", "Phase")#, "clonotype_id" "BCL2_K22K", "BCL2_L23L", "CD79B_Y696H")# c("seurat_clusters", "Condition", "Phénotype", "Phase", "K29Q", "L37M", "M11I", "pGln45")
  singlet@meta.data$nFeature_HTO <- NULL
  singlet@meta.data$HTO_secondID <- NULL
  singlet@meta.data$HTO_classification <- NULL
  singlet@meta.data$RNA_snn_res.0.5 <- NULL
  singlet@meta.data$hash.ID <- NULL
  singlet@meta.data$Greffe <- as.character(singlet@meta.data$Condition)   # Dissociation Pré/Post greffe
  singlet@meta.data[singlet@meta.data$Condition == "RCHOP", "Greffe" ] <- "Post-Greffe"
  singlet@meta.data[singlet@meta.data$Condition == "Excipient", "Greffe" ] <- "Post-Greffe"
  singlet@meta.data[singlet@meta.data$Condition == "Pré-greffe", "Greffe" ] <- "Pré-Greffe"
  #colnames(singlet@meta.data)[which(colnames(singlet@meta.data)=="SingleR.calls")] <- "Phénotype"
  
  
  ##### -- Sunburst -- #### 
  make_sunburst_data <- function(met){
    df <- data.frame(ids=character(), labels=character(), parents=character(), values=integer(),stringsAsFactors=FALSE)
    
    df0 <- met %>% mutate(path = paste(patient,  sep=";")) %>% dplyr::select(path, value)
    df1 <- met %>% mutate(path = paste(patient, etat,  sep=";")) %>% dplyr::select(path, value)
    df2 <- met %>% mutate(path = paste(patient, etat, condition,  sep=";")) %>% dplyr::select(path, value)
    df3 <- met %>% mutate(path = paste(patient, etat, condition, phénotype, sep=";")) %>% dplyr::select(path, value)
    df4 <- met %>% mutate(path = paste(patient, etat, condition, phénotype, sub, sep=";")) %>% dplyr::select(path, value)
    df5 <- met %>% mutate(path = paste(patient, etat, condition, phénotype, sub, phase, sep=";")) %>% dplyr::select(path, value)
    
    for (dfX in c(df0,df1,df2,df3,df4,df5)) {
      xnames <- row.names(table(dfX))
      max <- 0
      
      for (rows in xnames) {
        orga = strsplit(rows,";")
        if (max < length(orga[[1]])){max = length(orga[[1]])}}
      
      for (i in 1:max){
        for (rows in xnames) {
          orga = strsplit(rows,";")
          if ((length(orga[[1]]))==i){
            c = ""
            c2 = ""
            for (j in 1:i){c = paste0(c,orga[[1]][j],"-")}
            for (k in 1:i-1){c2 = paste0(c2,orga[[1]][k],"-")}
            c = str_sub(c,1,-2)
            c2 = str_sub(c2,2,-2)
            df[rows,"labels"]= orga[[1]][i]
            df[rows,"ids"]=c
            df[rows,"parents"]=c2
          }
        }
      }
    }
    df <- df[-2,]
    df$values <- c(table(df0)[,1],table(df1)[,1],table(df2)[,1],table(df3)[,1],table(df4)[,1],table(df5)[,1])
    return(df)
  }
  met <- data.frame(
    patient = singlet@meta.data$orig.ident, 
    etat = singlet@meta.data$Greffe, 
    condition = singlet@meta.data$Condition, 
    phénotype = singlet@meta.data$Phénotype, 
    sub = singlet@meta.data$Phénotype.fine, 
    phase = singlet@meta.data$Phase,
    value = rep(1, length(singlet@meta.data$orig.ident)), stringsAsFactors = FALSE
  )
  singlet@tools$sunburst <- make_sunburst_data(met)
  
  
  ##### -- DE -- ##### 
  Idents(singlet)<-"Condition"
  singlet[["RNA"]]@counts <- as.matrix(singlet[["RNA"]]@counts)+1
  singlet@tools$DE_RE <- FindMarkers(singlet, slot = "counts", ident.1 = "RCHOP", ident.2 = "Excipient", test.use = "DESeq2")
  singlet@tools$DE_PE <- FindMarkers(singlet, slot = "counts", ident.1 = "Pré-greffe", ident.2 = "Excipient", test.use = "DESeq2")
  singlet[["RNA"]]@counts <- as.matrix(singlet[["RNA"]]@counts)-1
  Idents(singlet)<-"seurat_clusters"
  
  
  ##### -- Linear correlation -- ##### 
  linear_correlation <- function(ctrl, stim){
    ctrl$stim <- "CTRL" ; stim$stim <- "STIM"
    immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20, k.filter=80)
    immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20, k.weight=80)
    DefaultAssay(immune.combined) <- "integrated"
    immune.combined <- ScaleData(immune.combined, verbose = FALSE)
    immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
    immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
    immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
    immune.combined <- FindClusters(immune.combined, resolution = 0.5)
    
    Idents(immune.combined)<-"Phénotype" ; DefaultAssay(immune.combined) <- "RNA"
    nk.markers <- FindConservedMarkers(immune.combined, ident.1 = "B-cells", grouping.var = "stim", verbose = FALSE)
    
    b.cells <- subset(immune.combined, idents = "B-cells") ; Idents(b.cells) <- "stim"
    return(as.data.frame(log1p(AverageExpression(b.cells, verbose = FALSE)$RNA)))
  }
  Excipient <- seurat_subset(singlet, "Condition", "Excipient")
  RCHOP <- seurat_subset(singlet, "Condition", "RCHOP")
  Pregreffe <- seurat_subset(singlet, "Condition", "Pré-greffe")
  
  singlet@tools$avg.b.cells_RE <- linear_correlation(Excipient, RCHOP)
  singlet@tools$avg.b.cells_PE <- linear_correlation(Pregreffe, Excipient)
  
  
  return(singlet)
}
seurat_object <- function(patient){
  ##### -- Préprocessing -- #####
  ##### -- Matrice mRNA -- #####
  rwa.mRNA <- Read10X(data.dir = paste0(patient,"/mRNA/CellrangerCount/outs/filtered_feature_bc_matrix/")) 
  cso <- CreateSeuratObject(counts = rwa.mRNA, project = patient, min.cells = 3, min.features = 200)
  umis <- GetAssayData(object = cso, slot = "counts")
  
  ##### -- Matrice HTO -- #####
  raw.hto <- Read10X(paste0(patient,"/HTO/result/umi_count/"), gene.column = 1)        
  colnames(raw.hto) <- paste0(colnames(raw.hto),"-1")                          
  hto <- raw.hto[c(1:3),]                                                       # Suppression des séquences unmapped : rownames(raw.hto)
  rownames(hto) <- c("Pré-greffe","Excipient","RCHOP")    
  joint.bcs <- intersect(colnames(umis),colnames(hto))                          # Sélection des cellules avec barcode commun HTO / mRNA
  umis <- umis[, joint.bcs]                                                     # Sélection des lignes qui correspondent aux cellules en commun
  hto <- as.matrix(hto[, joint.bcs])
  
  ##### -- Object seurat -- #####
  hashtag <- CreateSeuratObject(counts = umis, assay = "RNA", project = patient)
  hashtag <- NormalizeData(hashtag)
  hashtag <- FindVariableFeatures(hashtag, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  hashtag <- ScaleData(hashtag,features = VariableFeatures(hashtag))
  hashtag[["HTO"]] <- CreateAssayObject(counts = hto)                           # Ajoute des données HTO comme un nouvel assay indépendant du mRNA
  hashtag <- NormalizeData(hashtag, assay = "HTO", normalization.method = "CLR")
  hashtag <- HTODemux(hashtag, assay = "HTO", positive.quantile = 0.99)         # Association : cellules / échantillons
  Idents(hashtag)<- 'HTO_classification.global'
  singlet <- subset(hashtag, idents = "Singlet") 
  colnames(singlet@meta.data)[which(colnames(singlet@meta.data)=="HTO_maxID")] <- "Condition"
  
  # VlnPlot 1
  singlet[["percent.mt"]] <- PercentageFeatureSet(singlet, pattern = "^MT-")
  singlet@tools$mitochondrie_all <- VlnPlot(singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Condition")
  
  return(singlet)
}
visualisation <- function(singlet){
  ## -- Pre-processing  -- ##
  singlet <- NormalizeData(singlet, normalization.method = "LogNormalize", scale.factor = 10000) 
  singlet <- FindVariableFeatures(singlet, selection.method = "vst", nfeatures = 2000)
  singlet <- ScaleData(singlet, features = rownames(singlet))                                               # Scale data :singlet@assays$RNA@scale.data, singlet[["RNA"]]@scale.data
  singlet <- RunPCA(singlet, features = VariableFeatures(singlet), ndims.print = 1:10, nfeatures.print = 30)# Reduction dimension
  singlet <- FindNeighbors(singlet, reduction = "pca", dims = 1:40, compute.SNN = T)
  singlet <- FindClusters(singlet, resolution = 0.5)                                                        # head(Idents(singlet), 10) 
  singlet <- RunUMAP(singlet, reduction = "pca", dims = 1:40)
  singlet <- RunTSNE(singlet, reduction = "pca", dims = 1:40)
  
  ## -- FindAllMarkers -- ## 
  singlet@commands[["FindAllMarkers"]] <- FindAllMarkers(singlet, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  
  return(singlet)
}
processing <- function(patient){
  singlet <- seurat_object(patient)
  singlet <- visualisation(singlet)
  singlet <- metadata(singlet)
  return(singlet)
}

## -- Merge -- ##
metadata_merge <- function(singlet){
  for (metadata in colnames(singlet@meta.data)) {
    Idents(singlet)<- metadata
    singlet<- StashIdent(singlet, save.name = metadata)
  }
  return(singlet)
}
numeric_merge <- function(singlet){
  for (metadata in colnames(singlet@meta.data)) {
    if (metadata %in% singlet@tools$hallmarks){
      singlet@meta.data[[metadata]]<- as.numeric(singlet@meta.data[[metadata]])
    }
  }
  return(singlet)
}

## -- App -- ##
seurat_subset <- function(singlet, item, sub_item){
  Idents(singlet) <- item
  sub_singlet <- subset(singlet, idents = sub_item)
  sub_singlet <- visualisation(sub_singlet)
  return(sub_singlet)
}
gene_subset <- function(singlet, gene, expression){
  expr <- FetchData(object = singlet, vars = gene)
  sub_singlet <- singlet[, which(x = expr > expression)]
  sub_singlet <- visualisation(sub_singlet)
  return(sub_singlet)
}
QC_subset <- function(singlet, maximum_sub, percent_mt_sub){
  if(maximum_sub==""){maximum_sub = 10000}
  if(percent_mt_sub==""){percent_mt_sub = 50}
  
  ## -- QC filtres-- ##
  Idents(singlet)<-"seurat_clusters"
  sub_singlet <- subset(singlet, subset = nFeature_RNA > 200 & nFeature_RNA < maximum_sub & percent.mt < percent_mt_sub)            # QC Filter : tester 15%
  sub_singlet <- visualisation(sub_singlet)
  
  return(sub_singlet)
}


## -- Other -- ## 

#library(dplyr)
#library(SingleR)
#library(stringr)
#library(celldex)
#library(monocle)
#library(escape)