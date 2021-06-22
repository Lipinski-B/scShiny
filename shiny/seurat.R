library(Seurat)
library(celldex)
library(SingleR)
library(escape)
library(dittoSeq)
library(stringr)
library(dplyr)

source(file = "/home/boris/Bureau/scShiny/shiny/functions.R")
setwd(dir = "/home/boris/Documents/lipinskib/flinovo/result/")

## -- Loading -- ## 
siege <- c("FL140304","FL12C1888","FL09C1164","FL08G0293")
patient <- siege[2]
#load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))

## -- Workflow -- ## 
all <- processing(patient)
save(all, file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))

## -- Méta patient -- ##
P1 <- all ; rownames(P1@meta.data) <- paste0("FL14_",rownames(P1@meta.data))
P2 <- all ; rownames(P2@meta.data) <- paste0("FL12_",rownames(P2@meta.data))
all <- merge(P1, y = P2, add.cell.ids = c("FL12", "FL14"), project = "FL")

all <- visualitation(all)
all <- metadata(all)
save(all, file = paste0("/home/boris/Documents/analyse/singlet_FULL.RData"))







all <- seurat_object_namanm(1)
all <- visualitation(all)
all <- metadata_namnam(all)

save(all, file = paste0("/home/boris/Documents/analyse/1_namnam.RData"))

all <- seurat_object_namanm(2)
all <- visualitation(all)
all <- metadata_namnam(all)

save(all, file = paste0("/home/boris/Documents/analyse/2_namnam.RData"))

all <- seurat_object_namanm(3)
all <- visualitation(all)
all <- metadata_namnam(all)

save(all, file = paste0("/home/boris/Documents/analyse/3_namnam.RData"))

