library(Seurat)
library(celldex)
library(SingleR)
library(escape)
library(dittoSeq)
library(stringr)
library(dplyr)

metadata <- function(singlet){  
  ##### -- identification des cellules -- ##### 
  # DB: DatabaseImmuneCellExpressionData()
  #hpca.se <- celldex::DatabaseImmuneCellExpressionData()
  #results <- SingleR(test = as.SingleCellExperiment(singlet), ref = hpca.se, labels = hpca.se$label.main)
  #results.fine <- SingleR(test = as.SingleCellExperiment(singlet), ref = hpca.se, labels = hpca.se$label.fine)
  
  #singlet$SingleR.pruned.calls <- results$pruned.labels
  #singlet$SingleR.calls <- results$labels
  #singlet$SingleR.pruned.calls.fine <- results.fine$pruned.labels
  #singlet$SingleR.calls.fine <- results.fine$labels
  
  # DB: BlueprintEncodeData()
  BD <- celldex::BlueprintEncodeData()
  results.blueprint <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.main)
  results.blueprint.fine <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.fine)
  
  #singlet$SingleR.pruned.calls.blueprint <- results$pruned.labels
  singlet$SingleR.calls <- results.blueprint$labels
  #singlet$SingleR.pruned.calls.blueprint.fine <- results.blueprint.fine$pruned.labels
  singlet$SingleR.calls.fine <- results.blueprint.fine$labels
  
  # Récap
  table(results.blueprint$labels)
  table(results.blueprint.fine$labels)
  
  
  ##### -- Enrichissement des gènes -- ##### 
  GS <- getGeneSets(library = "H")
  ES <- enrichIt(obj = singlet, gene.sets = GS, groups = 1000, cores = 12)
  names(ES) <- str_replace_all(names(ES), "HALLMARK_", "")
  singlet <- AddMetaData(singlet, ES)
  singlet@tools$hallmarks <- names(ES)
  singlet@tools$meta_variable <- c("seurat_clusters", "Condition", "Greffe", "Phénotype", "clonotype_id", "Phase", "BCL2_K22K", "BCL2_L23L", "CD79B_Y696H")# c("seurat_clusters", "Condition", "Phénotype", "Phase", "K29Q", "L37M", "M11I", "pGln45")
  
  
  ##### -- VDJ -- ##### 
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
  raw <- rbind(bcr[bcr$clonotype_id == 'clonotype1',],bcr[bcr$clonotype_id == 'clonotype2',],bcr[bcr$clonotype_id == 'clonotype3',],bcr[bcr$clonotype_id == 'clonotype4',])
  bcr$clonotype_id <- as.numeric(str_replace(unlist(bcr$clonotype_id,"clonotype",use.names=F), "clonotype", ""))
  
  # Add BCR to metadata
  singlet <- AddMetaData(object=singlet, metadata = bcr)                        #, col.name = colnames(as.data.frame(bcr)))   # Add to the Seurat object's metadata.
  
  
  ##### -- GoT  -- ##### 
  load(file=paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient, "/GOT/result/", patient, "_GoT.Rdata"))
  singlet <- AddMetaData(object = singlet, metadata = GOT)
  
  
  ##### -- Phase cycle -- #### 
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  singlet <- CellCycleScoring(singlet, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)    # Attribution à chaque cellule d'un état de phase
  #singlet <- ScaleData(singlet, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(singlet)) # Correction de la matrice d'expression dépendemment du cycle cellulaire
  Idents(singlet)<-"seurat_clusters"
  
  ## -- Ohter -- ## 
  singlet@meta.data$nFeature_HTO <- NULL
  singlet@meta.data$HTO_secondID <- NULL
  singlet@meta.data$HTO_classification <- NULL
  singlet@meta.data$RNA_snn_res.0.5 <- NULL
  singlet@meta.data$hash.ID <- NULL
  singlet@meta.data$Greffe <- as.character(singlet@meta.data$Condition)   # Dissociation Pré/Post greffe
  singlet@meta.data[singlet@meta.data$Condition == "RCHOP", "Greffe" ] <- "Post-Greffe"
  singlet@meta.data[singlet@meta.data$Condition == "Excipient", "Greffe" ] <- "Post-Greffe"
  singlet@meta.data[singlet@meta.data$Condition == "Pré-greffe", "Greffe" ] <- "Pré-Greffe"
  colnames(singlet@meta.data)[which(colnames(singlet@meta.data)=="SingleR.calls")] <- "Phénotype"
  
  
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
    sub = singlet@meta.data$SingleR.calls.fine, 
    phase = singlet@meta.data$Phase,
    value = rep(1, length(singlet@meta.data$orig.ident)), stringsAsFactors = FALSE
  )
  singlet@tools$sunburst <- make_sunburst_data(met)
  
  ##### -- VDJ -- #### 
  load(file=paste0("/home/boris/Documents/lipinskib/flinovo/result/",patient,"/R/VDJ.RData"))
  
  singlet@tools$Clonotype <- Clonotype
  singlet@tools$Type <- Type
  singlet@tools$Isotype <- Isotype
  singlet@tools$V <- V
  singlet@tools$D <- D
  singlet@tools$J <- J
  singlet@tools$vloupe <- vloupe
  
  return(singlet)
}
seurat_object <- function(patient){
  ##### -- Préprocessing -- #####
  # Matrice mRNA + object seurat
  rwa.mRNA <- Read10X(data.dir = paste0(patient,"/mRNA/",patient,"_CellrangerCount/outs/filtered_feature_bc_matrix/")) 
  cso <- CreateSeuratObject(counts = rwa.mRNA, project = patient, min.cells = 3, min.features = 200)
  umis <- GetAssayData(object = cso, slot = "counts")
  
  # Matrice HTO
  raw.hto <- Read10X(paste0(patient,"/HTO/umi_count/"), gene.column = 1)        
  colnames(raw.hto) <- paste0(colnames(raw.hto),"-1")                          
  hto <- raw.hto[c(1:3),]                                                       # Suppression des séquences unmapped : rownames(raw.hto)
  rownames(hto) <- c("Pré-greffe","Excipient","RCHOP")    
  
  joint.bcs <- intersect(colnames(umis),colnames(hto))                          # Sélection des cellules avec barcode commun HTO / mRNA
  
  umis <- umis[, joint.bcs]                                                     # Sélection des lignes qui correspondent aux cellules en commun
  hto <- as.matrix(hto[, joint.bcs])
  
  # Object seurat : HTOs + mRNA
  hashtag <- CreateSeuratObject(counts = umis, assay = "RNA", project = patient)
  hashtag <- NormalizeData(hashtag)
  hashtag <- FindVariableFeatures(hashtag, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  hashtag <- ScaleData(hashtag,features = VariableFeatures(hashtag))
  
  hashtag[["HTO"]] <- CreateAssayObject(counts = hto)                           # Ajoute des données HTO comme un nouvel assay indépendant du mRNA
  hashtag <- NormalizeData(hashtag, assay = "HTO", normalization.method = "CLR")
  
  # Association : cellules / échantillons
  hashtag <- HTODemux(hashtag, assay = "HTO", positive.quantile = 0.99)
  table(hashtag$HTO_maxID)  
  table(hashtag$HTO_classification.global)                                      # Result
  
  Idents(hashtag)<- 'HTO_classification.global'
  
  singlet <- subset(hashtag, idents = "Singlet")
  colnames(singlet@meta.data)[which(colnames(singlet@meta.data)=="HTO_maxID")] <- "Condition"
  
  singlet[["percent.mt"]] <- PercentageFeatureSet(singlet, pattern = "^MT-")
  singlet@tools$mitochondrie_all <- VlnPlot(singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Condition")
  
  return(singlet)
}
visualitation <- function(singlet, maximum=2500, percent_mt=5){
  ## -- QC filtres-- ##
  Idents(singlet)<-"seurat_clusters"
  singlet <- subset(singlet, subset = nFeature_RNA > 200 & nFeature_RNA < maximum & percent.mt < percent_mt)            # QC Filter : tester 15%
  
  ## -- Pre-processing  -- ##
  singlet <- NormalizeData(singlet, normalization.method = "LogNormalize", scale.factor = 10000) 
  singlet <- FindVariableFeatures(singlet, selection.method = "vst", nfeatures = 2000)
  singlet <- ScaleData(singlet, features = rownames(singlet))                                               # Scale data :singlet@assays$RNA@scale.data, singlet[["RNA"]]@scale.data
  singlet <- RunPCA(singlet, features = VariableFeatures(singlet), ndims.print = 1:10, nfeatures.print = 30)# Reduction dimension
  singlet <- FindNeighbors(singlet, reduction = "pca", dims = 1:40, compute.SNN = T)
  singlet <- FindClusters(singlet, resolution = 0.5)                                                        # head(Idents(singlet), 10) 
  #singlet <- RunSPCA(singlet, features = VariableFeatures(singlet), ndims.print = 1:10, nfeatures.print = 30, graph = singlet@graphs[["RNA_snn"]])# Reduction dimension
  singlet <- RunUMAP(singlet, reduction = "pca", dims = 1:40)
  singlet <- RunTSNE(singlet, reduction = "pca", dims = 1:40)
  
  ## -- Ohter PCA analysis -- ## 
  #singlet <- JackStraw(singlet, num.replicate = 100, reduction = "pca")
  #singlet <- ScoreJackStraw(singlet, dims = 1:20)
  
  ## -- FindAllMarkers -- ## 
  #singlet@commands[["FindAllMarkers"]] <- FindAllMarkers(singlet, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  #singlet@commands[["FindAllMarkers"]] <- merge(singlet@commands[["FindAllMarkers"]], annotations, by.x="gene", by.y="gene_name")
  #singlet@commands[["FindAllMarkers"]] %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  
  return(singlet)
}
processing <- function(patient){
  singlet <- seurat_object(patient)
  singlet <- visualitation(singlet)
  singlet <- metadata(singlet)
  return(singlet)
}
seurat_subset <- function(singlet, item, sub_item, maximum_sub=2500, percent_mt_sub=5, subset=NULL, expression=NULL){
  Idents(singlet) <- item
  #if(is.null(subset)){sub_singlet <- subset(singlet, idents = sub_item)}else{sub_singlet <- subset(singlet, idents = sub_item, subset= subsets > expression)}
  sub_singlet <- subset(singlet, idents = sub_item)
  sub_singlet <- visualitation(sub_singlet, maximum=maximum_sub, percent_mt=percent_mt_sub)
  return(sub_singlet)
}
gene_subset <- function(singlet, gene, expression, maximum_sub=2500, percent_mt_sub=5){
  expr <- FetchData(object = singlet, vars = gene)
  sub_singlet <- singlet[, which(x = expr > expression)]
  sub_singlet <- visualitation(sub_singlet, maximum=maximum_sub, percent_mt=percent_mt_sub)
  return(sub_singlet)
}



seurat_object_namanm <- function(i){
  ##### -- Préprocessing -- #####
  # Matrice mRNA + object seurat
  rwa.mRNA <- Read10X(data.dir = paste0("/home/boris/Documents/lipinskib/out_nam/",i,"/outs/filtered_feature_bc_matrix/")) 
  cso <- CreateSeuratObject(counts = rwa.mRNA, project = as.character(i), min.cells = 3, min.features = 200)
  umis <- GetAssayData(object = cso, slot = "counts")
  
  # Matrice HTO
  #raw.hto <- Read10X(paste0(patient,"/HTO/umi_count/"), gene.column = 1)        
  #colnames(raw.hto) <- paste0(colnames(raw.hto),"-1")                          
  #hto <- raw.hto[c(1:3),]                                                       # Suppression des séquences unmapped : rownames(raw.hto)
  #rownames(hto) <- c("Pré-greffe","Excipient","RCHOP")    
  
  #joint.bcs <- intersect(colnames(umis),colnames(hto))                          # Sélection des cellules avec barcode commun HTO / mRNA
  
  #umis <- umis[, joint.bcs]                                                     # Sélection des lignes qui correspondent aux cellules en commun
  #hto <- as.matrix(hto[, joint.bcs])
  
  # Object seurat : HTOs + mRNA
  hashtag <- CreateSeuratObject(counts = umis, assay = "RNA", project = as.character(i))
  hashtag <- NormalizeData(hashtag)
  hashtag <- FindVariableFeatures(hashtag, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  hashtag <- ScaleData(hashtag,features = VariableFeatures(hashtag))
  
  #hashtag[["HTO"]] <- CreateAssayObject(counts = hto)                           # Ajoute des données HTO comme un nouvel assay indépendant du mRNA
  #hashtag <- NormalizeData(hashtag, assay = "HTO", normalization.method = "CLR")
  
  # Association : cellules / échantillons
  #hashtag <- HTODemux(hashtag, assay = "HTO", positive.quantile = 0.99)
  #table(hashtag$HTO_maxID)  
  #table(hashtag$HTO_classification.global)                                      # Result
  
  Idents(hashtag)<- 'HTO_classification.global'
  
  singlet <- hashtag
  #colnames(singlet@meta.data)[which(colnames(singlet@meta.data)=="HTO_maxID")] <- "Condition"
  
  singlet[["percent.mt"]] <- PercentageFeatureSet(singlet, pattern = "^MT-")
  singlet@tools$mitochondrie_all <- VlnPlot(singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  return(singlet)
}
metadata_namnam <- function(singlet){  
  singlet <- all
  # DB: BlueprintEncodeData()
  BD <- celldex::BlueprintEncodeData()
  results.blueprint <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.main)
  results.blueprint.fine <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.fine)
  
  #singlet$SingleR.pruned.calls.blueprint <- results$pruned.labels
  singlet$SingleR.calls <- results.blueprint$labels
  #singlet$SingleR.pruned.calls.blueprint.fine <- results.blueprint.fine$pruned.labels
  singlet$SingleR.calls.fine <- results.blueprint.fine$labels
  
  
  ##### -- Enrichissement des gènes -- ##### 
  #GS <- getGeneSets(library = "H")
  #ES <- enrichIt(obj = singlet, gene.sets = GS, groups = 1000, cores = 12)
  #names(ES) <- str_replace_all(names(ES), "HALLMARK_", "")
  #singlet <- AddMetaData(singlet, ES)
  #singlet@tools$hallmarks <- names(ES)
  #singlet@tools$meta_variable <- c("seurat_clusters", "Condition", "Greffe", "Phénotype", "clonotype_id", "Phase", "BCL2_K22K", "BCL2_L23L", "CD79B_Y696H")# c("seurat_clusters", "Condition", "Phénotype", "Phase", "K29Q", "L37M", "M11I", "pGln45")
  
  
  ##### -- Phase cycle -- #### 
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  singlet <- CellCycleScoring(singlet, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)    # Attribution à chaque cellule d'un état de phase
  #singlet <- ScaleData(singlet, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(singlet)) # Correction de la matrice d'expression dépendemment du cycle cellulaire
  Idents(singlet)<-"seurat_clusters"
  
  ## -- Ohter -- ## 
  singlet@meta.data$nFeature_HTO <- NULL
  singlet@meta.data$HTO_secondID <- NULL
  singlet@meta.data$HTO_classification <- NULL
  singlet@meta.data$RNA_snn_res.0.5 <- NULL
  singlet@meta.data$hash.ID <- NULL
  colnames(singlet@meta.data)[which(colnames(singlet@meta.data)=="SingleR.calls")] <- "Phénotype"
  
  all<- singlet
  
  ##### -- Sunburst -- #### 
  make_sunburst_data <- function(met){
    df <- data.frame(ids=character(), labels=character(), parents=character(), values=integer(),stringsAsFactors=FALSE)
    
    df1 <- met %>% mutate(path = paste(patient,  sep=";")) %>% dplyr::select(path, value)
    df3 <- met %>% mutate(path = paste(patient, phénotype, sep=";")) %>% dplyr::select(path, value)
    df4 <- met %>% mutate(path = paste(patient, phénotype, sub, sep=";")) %>% dplyr::select(path, value)
    df5 <- met %>% mutate(path = paste(patient, phénotype, sub, phase, sep=";")) %>% dplyr::select(path, value)
    
    for (dfX in c(df1,df3,df4,df5)) {
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
    df$values <- c(table(df1)[,1],table(df3)[,1],table(df4)[,1],table(df5)[,1])
    return(df)
  }
  met <- data.frame(
    patient = singlet@meta.data$orig.ident, 
    phénotype = singlet@meta.data$Phénotype, 
    sub = singlet@meta.data$SingleR.calls.fine, 
    phase = singlet@meta.data$Phase,
    value = rep(1, length(singlet@meta.data$orig.ident)), stringsAsFactors = FALSE
  )
  singlet@tools$sunburst <- make_sunburst_data(met)
  
  
  return(singlet)
}



