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
patient <- siege[4]
#load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))


all <- seurat_object(patient)
saveRDS(all, file = paste0("/home/boris/Bureau/celloracle/singlet_", patient,".Rds"))

## -- Workflow -- ## 
all <- processing(patient)
save(all, file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))


## -- Méta patient -- ##
P1 <- all ; rownames(P1@meta.data) <- paste0("FL14_",rownames(P1@meta.data))
P2 <- all ; rownames(P2@meta.data) <- paste0("FL12_",rownames(P2@meta.data))
P3 <- all ; rownames(P3@meta.data) <- paste0("FL09_",rownames(P3@meta.data))
all <- merge(P1, y = c(P2,P3), add.cell.ids = c("FL14", "FL12","FL09"), project = "FL")

all <- visualitation(all)
all <- metadata(all)
save(all, file = paste0("/home/boris/Documents/analyse/singlet_RP.RData"))

