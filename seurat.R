library(celldex)
library(SingleR)
library(escape)
library(stringr)
library(dplyr)
library(future)
library(DESeq2)
library(cowplot)
source(file = "/home/boris/Bureau/scShiny/functions.R")

## -- Loading -- ## 
setwd(dir = "/home/boris/Documents/lipinskib/flinovo/result/")
siege <- c("FL12C1888","FL09C1164","FL02G095","FL05G0330", "FL140304","FL08G0293") #'all' 
patient <- siege[1]


## -- Workflow -- ##
#Full
singlet <- processing(patient)
save(singlet, file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
#save(all, file = paste0("/home/boris/Documents/lipinskib/flinovo/result/",patient,"/R/singlet_", patient,".RData"))

#app
load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
all@assays[["SCT"]]@scale.data <- subset(all@assays[["SCT"]]@scale.data, rownames(all@assays[["SCT"]]@scale.data) %in% c(rownames(all@tools$DE_PE)[1:50], rownames(all@tools$DE_RE)[1:50]))
singlet <- DietSeurat(all, counts = FALSE, data = T, scale.data = T,features = NULL, assays = NULL, dimreducs = c("pca","umap",'tsne'), graphs = NULL )
all@assays[["HTO"]] <- list()
all@assays[["RNA"]] <- list()
save(singlet, file = paste0("/home/boris/Bureau/scShiny/www/", patient,".RData"))

PCAPlot(object = singlet, label.size = 0.0, pt.size = 2) & theme(title = element_text(size=20),legend.position = "top",legend.title = element_text(size=10),legend.text = element_text(size=10)) & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 6))) & xlab(label = paste0("PCA 1 : ", round(Stdev(singlet[["pca"]])[1],2), " %")) & ylab(label = paste0("PCA 2 : ", round(Stdev(singlet[["pca"]])[2],2), " %"))
Idents(singlet)<-"Condition" ; DoHeatmap(singlet, cells = rownames(singlet@meta.data)[which(singlet@meta.data$Condition==c("Excipient","RCHOP"))], features = rownames(singlet@tools$DE_RE)[1:50], size = 3)



## -- Sub sample -- ##
patient <- "FL05G0330"
load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
all <- seurat_subset(singlet,"Condition", c("RCHOP","Pré-greffe")) ; save(all, file = paste0("/home/boris/Documents/analyse/singlet_", patient,"_PG_RCHOP.RData"))
all <- seurat_subset(singlet,"Condition", c("Excipient","Pré-greffe")) ; save(all, file = paste0("/home/boris/Documents/analyse/singlet_", patient,"_PG_EX.RData"))
all <- seurat_subset(singlet,"Condition", c("RCHOP","Excipient")) ; save(all, file = paste0("/home/boris/Documents/analyse/singlet_", patient,"_EX_RCHOP.RData"))



## -- Séparation Pré-greffe/Excipient -- ##
patient <- "FL05G0330"
load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
Idents(all) <- "Condition"
all <- subset(all, idents = c('Pré-greffe','Excipient'))
all <- subset(all, idents = 'B-cells')
all <- visualisation(all)

PCAPlot(object = all, label.size = 0.0, pt.size = 2) & 
  theme(title = element_text(size=20),legend.position = "top",legend.title = element_text(size=10),legend.text = element_text(size=10)) &guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 6))) & 
  xlab(label = paste0("PCA 1 : ", round(Stdev(all[["pca"]])[1],2), " %")) &ylab(label = paste0("PCA 2 : ", round(Stdev(all[["pca"]])[2],2), " %"))




## -- DE intersect -- ##
DE <- list()
for (patient in siege) {
  load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
  DE[[patient]]$DE_RE <- all@tools$DE_RE[which(all@tools$DE_RE$p_val_adj < 0.05),]
  DE[[patient]]$DE_PE <- all@tools$DE_PE[which(all@tools$DE_PE$p_val_adj < 0.05),]
}
liste = read.table("/home/boris/Bureau/dissoc.txt", header = T)
liste <- liste$x
X <- 500
DE_RE <- Reduce(intersect, list(rownames(DE[["FL140304"]]$DE_RE)[1:X],rownames(DE[["FL12C1888"]]$DE_RE)[1:X],rownames(DE[["FL09C1164"]]$DE_RE)[1:X],rownames(DE[["FL08G0293"]]$DE_RE)[1:X],rownames(DE[["FL02G095"]]$DE_RE)[1:X],rownames(DE[["FL05G0330"]]$DE_RE)[1:X]))
DE_PE <- Reduce(intersect, list(rownames(DE[["FL140304"]]$DE_PE)[1:X],rownames(DE[["FL12C1888"]]$DE_PE)[1:X],rownames(DE[["FL09C1164"]]$DE_PE)[1:X],rownames(DE[["FL08G0293"]]$DE_PE)[1:X],rownames(DE[["FL02G095"]]$DE_PE)[1:X],rownames(DE[["FL05G0330"]]$DE_PE)[1:X]))

DE_RE_candidat <- setdiff(DE_RE,liste)
DE_PE_candidat <- setdiff(DE_PE,liste)

DE_PE_FC <- data.frame(matrix(ncol = 40, nrow = 6))
rownames(DE_PE_FC) <- siege ; colnames(DE_PE_FC) <- DE_PE_candidat
for (patient in siege) {DE_PE_FC[patient,] <- DE[[patient]]$DE_PE[rownames(DE[[patient]]$DE_PE) %in% DE_PE_candidat,]$avg_log2FC}

write.table(DE_PE_FC, "/home/boris/Bureau/DE_PE_FC.txt")
write.table(DE_PE_candidat_FC, "/home/boris/Bureau/DE_PE_FC.txt")



## -- Méta patient -- ##
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[1],".RData")) ; P1 <- all ; rownames(P1@meta.data) <- paste0("FL14_",rownames(P1@meta.data)) ; P1 <- metadata_merge(P1)
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[2],".RData")) ; P2 <- all ; rownames(P2@meta.data) <- paste0("FL12_",rownames(P2@meta.data)) ; P2 <- metadata_merge(P2)
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[3],".RData")) ; P3 <- all ; rownames(P3@meta.data) <- paste0("FL09_",rownames(P3@meta.data)) ; P3 <- metadata_merge(P3)
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[4],".RData")) ; P4 <- all ; rownames(P4@meta.data) <- paste0("FL08_",rownames(P4@meta.data)) ; P4 <- metadata_merge(P4)
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[5],".RData")) ; P5 <- all ; rownames(P5@meta.data) <- paste0("FL02_",rownames(P5@meta.data)) ; P5 <- metadata_merge(P5)
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[6],".RData")) ; P6 <- all ; rownames(P6@meta.data) <- paste0("FL05_",rownames(P6@meta.data)) ; P6 <- metadata_merge(P6)

all <- merge(P1, y = c(P2,P3,P4,P5,P6), add.cell.ids = c("FL14","FL12","FL09","FL08","FL02","FL05"), project = "FL")
save(all, file = paste0("/home/boris/Documents/analyse/singlet_all.RData"))
load(file = paste0("/home/boris/Documents/analyse/singlet_all.RData"))
all <- visualisation(all)
all@tools$hallmarks <- c("ADIPOGENESIS","ALLOGRAFT_REJECTION","ANDROGEN_RESPONSE","ANGIOGENESIS","APICAL_JUNCTION","APICAL_SURFACE","APOPTOSIS","BILE_ACID_METABOLISM","CHOLESTEROL_HOMEOSTASIS","COAGULATION","COMPLEMENT","DNA_REPAIR","E2F_TARGETS","EPITHELIAL_MESENCHYMAL_TRANSITION","ESTROGEN_RESPONSE_EARLY","ESTROGEN_RESPONSE_LATE","FATTY_ACID_METABOLISM","G2M_CHECKPOINT","GLYCOLYSIS","HEDGEHOG_SIGNALING","HEME_METABOLISM","HYPOXIA","IL2_STAT5_SIGNALING","IL6_JAK_STAT3_SIGNALING","INFLAMMATORY_RESPONSE","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","KRAS_SIGNALING_DN","KRAS_SIGNALING_UP","MITOTIC_SPINDLE","MTORC1_SIGNALING","MYC_TARGETS_V1","MYC_TARGETS_V2","MYOGENESIS","NOTCH_SIGNALING","OXIDATIVE_PHOSPHORYLATION","P53_PATHWAY","PANCREAS_BETA_CELLS","PEROXISOME","PI3K_AKT_MTOR_SIGNALING","PROTEIN_SECRETION","REACTIVE_OXYGEN_SPECIES_PATHWAY","SPERMATOGENESIS","TGF_BETA_SIGNALING","TNFA_SIGNALING_VIA_NFKB","UNFOLDED_PROTEIN_RESPONSE","UV_RESPONSE_DN","UV_RESPONSE_UP","WNT_BETA_CATENIN_SIGNALING","XENOBIOTIC_METABOLISM")
all@tools$meta_variable <- c("seurat_clusters", "Condition", "Greffe", "Phénotype", "clonotype_id", "Phase", "old.ident", "BCL2_L23L","BCL2_K22K","CD79B_Y696H","EZH2_A682G_1","EZH2_A682G_2","EZH2_A692V_1","EZH2_A692V_2","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S","EZH2_A692V","EZH2_A682G")#, "BCL2_K22K", "BCL2_L23L", "CD79B_Y696H")# c("seurat_clusters", "Condition", "Phénotype", "Phase", "K29Q", "L37M", "M11I", "pGln45")
all <- numeric_merge(all)
save(all, file = paste0("/home/boris/Documents/analyse/singlet_all.RData"))

Idents(all)<-"orig.ident"
PCAPlot(object = all, group.by = "Phénotype", label.size = 0.0, pt.size = 2) & theme(title = element_text(size=20),legend.position = "top",legend.title = element_text(size=10),legend.text = element_text(size=10)) & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 6))) & xlab(label = paste0("PCA 1 : ", round(Stdev(all[["pca"]])[1],2), " %")) & ylab(label = paste0("PCA 2 : ", round(Stdev(all[["pca"]])[2],2), " %"))




## -- Réponse : RP/ RC -- ##
load(file = paste0("/home/boris/Documents/analyse/singlet_all.RData"))
all <- seurat_subset(all, "Condition", c("Pré-greffe"))
all@meta.data[all@meta.data$orig.ident == "FL140304", "Reponse"] <- "RP"
all@meta.data[all@meta.data$orig.ident == "FL12C1888", "Reponse"] <- "RP"
all@meta.data[all@meta.data$orig.ident == "FL08G0293", "Reponse"] <- "RP"

all@meta.data[all@meta.data$orig.ident == "FL09C1164", "Reponse"] <- "RC"
all@meta.data[all@meta.data$orig.ident == "FL02G095", "Reponse"] <- "RC"
all@meta.data[all@meta.data$orig.ident == "FL05G0330", "Reponse"] <- "RC"

Idents(all)<-"Reponse" ; all[["RNA"]]@counts <- as.matrix(all[["RNA"]]@counts)+1
all@tools$DE_R <- FindMarkers(all, slot = "counts", ident.1 = "RP", ident.2 = "RC", test.use = "DESeq2", max.cells.per.ident = 2000)
all[["RNA"]]@counts <- as.matrix(all[["RNA"]]@counts)-1 ; Idents(all)<-"seurat_clusters"
save(all, file = "/home/boris/Documents/analyse/singlet_all_PG.RData")





## -- TEST -- ##
liste = as.matrix(rownames(all@meta.data)[which(all@meta.data$APOPTOSIS>0.07 & all@meta.data$Phénotype=="B-cells" & all@meta.data$Condition ==c("RCHOP"))])
liste2 = as.matrix(rownames(all@meta.data)[which(all@meta.data$APOPTOSIS>0.13 & all@meta.data$Phénotype=="B-cells" & all@meta.data$Condition ==c("Excipient"))])
sub_all <- subset(all, cells = c(liste,liste2))

singlet <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000) 
singlet <- FindVariableFeatures(singlet, selection.method = "vst", nfeatures = 2000)
singlet <- ScaleData(singlet, features = rownames(singlet))                                               # Scale data :singlet@assays$RNA@scale.data, singlet[["RNA"]]@scale.data
singlet <- RunPCA(singlet, features = VariableFeatures(singlet), ndims.print = 1:10, nfeatures.print = 30)# Reduction dimension
singlet <- FindNeighbors(singlet, reduction = "pca", dims = 1:40, compute.SNN = T)
singlet <- FindClusters(singlet, resolution = 0.5)                                                        # head(Idents(singlet), 10) 
singlet <- RunUMAP(singlet, reduction = "pca", dims = 1:40)

plots <- UMAPPlot(object = singlet, group.by = "Condition", split.by = NULL, label.size = 0.0, pt.size = 2)
plots & theme(title = element_text(size=20),legend.position = "top",legend.title = element_text(size=10),legend.text = element_text(size=10)
) & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 6))) & xlab(label = paste0("PCA 1 : ", round(Stdev(singlet[["pca"]])[1],2), " %")) & ylab(label = paste0("PCA 2 : ", round(Stdev(singlet[["pca"]])[2],2), " %"))

liste = as.matrix(rownames(all@meta.data)[which(all@meta.data$APOPTOSIS>0.18 & all@meta.data$Phénotype=="B-cells" & all@meta.data$Condition ==c("RCHOP"))])
liste2 = as.matrix(rownames(all@meta.data)[which(all@meta.data$APOPTOSIS>0.18 & all@meta.data$Phénotype=="B-cells" & all@meta.data$Condition ==c("Excipient"))])

all <- seurat_subset(all, "Phénotype", )
apop <- as.vector(read.table("/home/boris/Bureau/gene_apop.txt"))
DoHeatmap(singlet, group.by = "Condition", feature = apop$V1)






## -- TEST : DEenrichRPlot -- ##
patient <- "FL05G0330"
load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
DEenrichRPlot(all, slot = "counts", ident.1 = "Excipient", ident.2 = "RCHOP", test.use = "DESeq2", max.cells.per.ident = 2000, balanced=T, p.val.cutoff=0.05, return.gene.list)


## -- TEST : DE appareillé : RCHOP/EXP + PG : RC/RP -- ##



## -- TEST : DE RCHOP : apop+ / apop- -- ##
patient <- "FL05G0330"
load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
liste = as.matrix(rownames(all@meta.data)[which(all@meta.data$APOPTOSIS>0.07 & all@meta.data$Phénotype=="B-cells" & all@meta.data$Condition ==c("RCHOP"))])
liste2 = as.matrix(rownames(all@meta.data)[which(all@meta.data$APOPTOSIS>0.13 & all@meta.data$Phénotype=="B-cells" & all@meta.data$Condition ==c("Excipient"))])
sub_all <- subset(all, cells = c(liste,liste2))