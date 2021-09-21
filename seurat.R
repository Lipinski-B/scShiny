library(celldex)
library(SingleR)
library(escape)
library(stringr)
library(dplyr)
library(future)
library(DESeq2)
library(cowplot)
library(Seurat)
library(dittoSeq)
library(DT)

library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE

source(file = "/home/boris/Bureau/scShiny/functions.R")

## -- Loading -- ## 
setwd(dir = "/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/")
siege <- c("FL09C1164","FL02G095","FL05G0330","FL140304", "FL08G0293", "FL12C1888") #'all'
patient <- siege[5]
load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))


singlet@assays$RNA@meta.features <- data.frame()
save(singlet, file = paste0("/home/boris/Bureau/scShiny/www/", patient,"/", patient,".RData"))

################################################################################################################################################################################################################################################################################################################################################
## -- Workflow -- ##
singlet <- processing(patient) ; diet(singlet)
#save(singlet, file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
#save(singlet, file = paste0("/home/boris/Documents/lipinskib/flinovo/result/",patient,"/R/singlet_", patient,".RData"))
rm(singlet) ; gc() ; gc() ; gc()


## -- Sub sample -- ##
patient <- "FL05G0330"
load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
singlet <- seurat_subset(singlet,"Condition", c("RCHOP","Pré-greffe")) ; save(singlet, file = paste0("/home/boris/Documents/analyse/singlet_", patient,"_PG_RCHOP.RData"))
singlet <- seurat_subset(singlet,"Condition", c("Excipient","Pré-greffe")) ; save(singlet, file = paste0("/home/boris/Documents/analyse/singlet_", patient,"_PG_EX.RData"))
singlet <- seurat_subset(singlet,"Condition", c("RCHOP","Excipient")) ; save(singlet, file = paste0("/home/boris/Documents/analyse/singlet_", patient,"_EX_RCHOP.RData"))



################################################################################################################################################################################################################################################################################################################################################
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

## -- TEST : DE appareillé : RCHOP/EXP + PG : RC/RP -- ##
## test 1 = Réponse : RP/RC -- ##
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

## test 2 = Réponse : RP/RC -- ##




################################################################################################################################################################################################################################################################################################################################################
## -- DE intersect -- ##
liste = read.table("/home/boris/Bureau/dissoc.txt", header = T) ; liste <- liste$x
DE <- list()
for (patient in siege) {
  load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
  DE[[patient]]$DE_RE <- singlet@tools$DE_RE[which(singlet@tools$DE_RE$p_val_adj < 0.05),]
  #DE[[patient]]$DE_PE <- all@tools$DE_PE[which(all@tools$DE_PE$p_val_adj < 0.05),]
}
X <- 500
DE_RE <- Reduce(intersect, list(rownames(DE[["FL140304"]]$DE_RE)[1:X],rownames(DE[["FL12C1888"]]$DE_RE)[1:X],rownames(DE[["FL08G0293"]]$DE_RE)[1:X]))
DE_PE <- Reduce(intersect, list(rownames(DE[["FL140304"]]$DE_PE)[1:X],rownames(DE[["FL12C1888"]]$DE_PE)[1:X],rownames(DE[["FL09C1164"]]$DE_PE)[1:X],rownames(DE[["FL08G0293"]]$DE_PE)[1:X],rownames(DE[["FL02G095"]]$DE_PE)[1:X],rownames(DE[["FL05G0330"]]$DE_PE)[1:X]))

DE_RE <- Reduce(intersect, list(rownames(DE[["FL05G0330"]]$DE_RE)[1:X],rownames(DE[["FL02G095"]]$DE_RE)[1:X],rownames(DE[["FL09C1164"]]$DE_RE)[1:X]))

DE_RE_candidat <- setdiff(DE_RE,liste)
DE_PE_candidat <- setdiff(DE_PE,liste)

DE_PE_FC <- data.frame(matrix(ncol = 40, nrow = 6))
rownames(DE_PE_FC) <- siege ; colnames(DE_PE_FC) <- DE_PE_candidat
for (patient in siege) {DE_PE_FC[patient,] <- DE[[patient]]$DE_PE[rownames(DE[[patient]]$DE_PE) %in% DE_PE_candidat,]$avg_log2FC}

write.table(DE_PE_FC, "/home/boris/Bureau/DE_PE_FC.txt")
write.table(DE_PE_candidat_FC, "/home/boris/Bureau/DE_PE_FC.txt")





################################################################################################################################################################################################################################################################################################################################################
## -- TEST : DE RCHOP : apop+ / apop- -- ##
patient <- "FL05G0330"
load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
seuil = 0.13
liste = as.matrix(rownames(singlet@meta.data)[which(singlet@meta.data$APOPTOSIS<seuil & singlet@meta.data$Phénotype=="B-cells" & singlet@meta.data$Condition =="RCHOP")])
liste2 = as.matrix(rownames(singlet@meta.data)[which(singlet@meta.data$APOPTOSIS>seuil & singlet@meta.data$Phénotype=="B-cells" & singlet@meta.data$Condition =="RCHOP")])
singlet@meta.data[singlet@meta.data$Condition == "RCHOP" & singlet@meta.data$orig.ident %in% liste, "Apop" ] <- "+"
singlet@meta.data[singlet@meta.data$Condition == "RCHOP" & singlet@meta.data$orig.ident %in% liste2, "Apop" ] <- "-"
