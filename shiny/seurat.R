library(Seurat)
library(celldex)
library(SingleR)
library(escape)
library(dittoSeq)
library(stringr)
library(dplyr)
library(future)
library(DESeq2)

## -- Loading -- ## 
source(file = "/home/boris/Bureau/scShiny/shiny/functions.R")
setwd(dir = "/home/boris/Documents/lipinskib/flinovo/result/")
siege <- c("FL140304","FL12C1888","FL09C1164","FL08G0293") #'all'
patient <- siege[2]
load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData")) 
save(all, file = paste0("/home/boris/Documents/lipinskib/flinovo/result/",patient,"/R/singlet_", patient,".RData"))




## -- Workflow -- ## 
#all <- processing(patient)
#save(all, file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))

DE <- list()
for (patient in siege) {
  load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData")) 
  DE[[patient]]$DE_RE <- all@tools$DE_RE
  DE[[patient]]$DE_PE <- all@tools$DE_PE
}


Reduce(intersect, list(rownames(DE[["FL140304"]]$DE_RE)[1:50],
                       rownames(DE[["FL12C1888"]]$DE_RE)[1:50],
                       rownames(DE[["FL09C1164"]]$DE_RE)[1:50],
                       rownames(DE[["FL08G0293"]]$DE_RE)[1:50]))

Reduce(intersect, list(rownames(DE[["FL140304"]]$DE_PE)[1:50],
                       rownames(DE[["FL12C1888"]]$DE_PE)[1:50],
                       rownames(DE[["FL09C1164"]]$DE_PE)[1:50],
                       rownames(DE[["FL08G0293"]]$DE_PE)[1:50]))

#FeaturePlot(all2, features = features)
#FeaturePlot(all2, features = c("MS4A1", "CD79A"), blend = TRUE)

# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
DotPlot(pbmc3k.final, features = features, split.by = "groups") + RotatedAxis()



## -- Méta patient -- ##
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[1],".RData")) ; P1 <- all ; rownames(P1@meta.data) <- paste0("FL14_",rownames(P1@meta.data)) ; P1 <- metadata_merge(P1)
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[2],".RData")) ; P2 <- all ; rownames(P2@meta.data) <- paste0("FL12_",rownames(P2@meta.data)) ; P2 <- metadata_merge(P2)
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[3],".RData")) ; P3 <- all ; rownames(P3@meta.data) <- paste0("FL09_",rownames(P3@meta.data)) ; P3 <- metadata_merge(P3)
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[4],".RData")) ; P4 <- all ; rownames(P4@meta.data) <- paste0("FL08_",rownames(P4@meta.data)) ; P4 <- metadata_merge(P4)

all <- merge(P1, y = c(P2,P3,P4), add.cell.ids = c("FL14","FL12","FL09","FL08"), project = "FL")
all <- visualisation(all)
all@tools$hallmarks <- c("ADIPOGENESIS","ALLOGRAFT_REJECTION","ANDROGEN_RESPONSE","ANGIOGENESIS","APICAL_JUNCTION","APICAL_SURFACE","APOPTOSIS","BILE_ACID_METABOLISM","CHOLESTEROL_HOMEOSTASIS","COAGULATION","COMPLEMENT","DNA_REPAIR","E2F_TARGETS",
                             "EPITHELIAL_MESENCHYMAL_TRANSITION","ESTROGEN_RESPONSE_EARLY","ESTROGEN_RESPONSE_LATE","FATTY_ACID_METABOLISM","G2M_CHECKPOINT","GLYCOLYSIS","HEDGEHOG_SIGNALING","HEME_METABOLISM","HYPOXIA","IL2_STAT5_SIGNALING","IL6_JAK_STAT3_SIGNALING",
                             "INFLAMMATORY_RESPONSE","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","KRAS_SIGNALING_DN","KRAS_SIGNALING_UP","MITOTIC_SPINDLE","MTORC1_SIGNALING","MYC_TARGETS_V1","MYC_TARGETS_V2","MYOGENESIS","NOTCH_SIGNALING",
                             "OXIDATIVE_PHOSPHORYLATION","P53_PATHWAY","PANCREAS_BETA_CELLS","PEROXISOME","PI3K_AKT_MTOR_SIGNALING","PROTEIN_SECRETION","REACTIVE_OXYGEN_SPECIES_PATHWAY","SPERMATOGENESIS",
                             "TGF_BETA_SIGNALING","TNFA_SIGNALING_VIA_NFKB","UNFOLDED_PROTEIN_RESPONSE","UV_RESPONSE_DN","UV_RESPONSE_UP","WNT_BETA_CATENIN_SIGNALING","XENOBIOTIC_METABOLISM")
all@tools$meta_variable <- c("seurat_clusters", "Condition", "Greffe", "Phénotype", "clonotype_id", "Phase", "old.ident", "BCL2_L23L","BCL2_K22K","CD79B_Y696H","EZH2_A682G_1","EZH2_A682G_2","EZH2_A692V_1","EZH2_A692V_2","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S","EZH2_A692V","EZH2_A682G")#, "BCL2_K22K", "BCL2_L23L", "CD79B_Y696H")# c("seurat_clusters", "Condition", "Phénotype", "Phase", "K29Q", "L37M", "M11I", "pGln45")
all <- numeric_merge(all)
save(all, file = paste0("/home/boris/Documents/analyse/singlet_all.RData"))

## -- Méta subset -- ##
all <- seurat_subset(all, "Condition", c("Excipient","Pré-greffe"))
GS <- getGeneSets(library = "H")
ES <- enrichIt(obj = all, gene.sets = GS, groups = 400, cores = 10)
names(ES) <- str_replace_all(names(ES), "HALLMARK_", "")
all <- AddMetaData(all, ES)
all@commands[["FindAllMarkers"]] <- FindAllMarkers(all, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)