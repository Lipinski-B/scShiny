.libPaths( c( .libPaths(), '/home/boris/R/x86_64-pc-linux-gnu-library/4.0') )
library(Seurat)
library(monocle)
library(escape)
library(dittoSeq)
library(plotly)
library(dplyr)
library(shiny)
library(shinyjs)
library(shinybusy)
library(shinyWidgets)
library(shinydashboard)
library(dashboardthemes)

source(file = "functions.R")
load(file = paste0("/home/boris/Documents/analyse/singlet_FL09C1164.RData"))

feature <- row.names(as.matrix(all[["RNA"]]@counts))
#all@meta.data$SingleR.calls.fine<- all@meta.data$d_gene <- all@meta.data[70] <- all@meta.data$v_gene <- all@meta.data$j_gene <- all@meta.data$c_gene <- all@meta.data$old.ident <- all@meta.data$BCL2_L23L <- all@meta.data$BCL2_K22K <- all@meta.data$CD79B_Y696H <- all@meta.data$EZH2_A682G_1 <- all@meta.data$EZH2_A682G_2  <- all@meta.data$EZH2_A692V_1  <- all@meta.data$EZH2_A692V_2  <- all@meta.data$EZH2_Y646S  <- all@meta.data$EZH2_Y646N  <- all@meta.data$EZH2_Y646H  <- all@meta.data$EZH2_Y646F  <- all@meta.data$EZH2_Y646C <- NULL
metadata <- c() ; List <- list()
for(i in 1:length(colnames(all@meta.data))){if (length(levels(as.factor(all@meta.data[[i]]))) > 1 && length(levels(as.factor(all@meta.data[[i]]))) < 600 && is.numeric(levels(as.factor(all@meta.data[[1]])))==F ){metadata <- c(metadata, colnames(all@meta.data)[i])}}
for(i in 1:length(metadata) ){List[[metadata[i]]] <- levels(as.factor(all@meta.data[[metadata[i]]]))}
singlet <- all

singlet@tools$meta_variable <- c("seurat_clusters","Condition","Greffe","Phénotype","clonotype_id","Phase","CD79B_Y696H")

#FL09 : réponse modéré

#annotations <- read.csv("/home/boris/Bureau/scShiny/annotation_FindAllMarkers.csv")
#Idents(singlet) <- "Phénotype"
#singlet@commands[["FindAllMarkers"]] <- FindAllMarkers(singlet, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
#invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)), detach, character.only = TRUE, unload = TRUE))

#IGHG1 ET SCREPERTOIRE, meta patient, rapport discussion, immunosect, fiche récap manon (docker, singularity, schiny, nextflow)