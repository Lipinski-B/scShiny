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
load(file = paste0("/home/boris/Documents/analyse/singlet_FL12C1888.RData"))

metadata <- c()
for(i in 1:length(colnames(all@meta.data))){if (length(levels(as.factor(all@meta.data[[i]]))) > 1 && length(levels(as.factor(all@meta.data[[i]]))) < 25 && is.numeric(levels(as.factor(all@meta.data[[1]])))==F ){metadata <- c(metadata, colnames(all@meta.data)[i])}}

List <- list()
for(i in 1:length(metadata) ){List[[metadata[i]]] <- levels(as.factor(all@meta.data[[metadata[i]]]))}

singlet <- all

#Idents(singlet) <- "Phénotype"
#singlet@commands[["FindAllMarkers"]] <- FindAllMarkers(singlet, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

#singlet@tools$meta_variable <- c("seurat_clusters", "Condition", "Phénotype", "Phase")#, "K29Q", "L37M", "M11I", "pGln45")

#invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)), detach, character.only = TRUE, unload = TRUE))


#library("uwot") # UMAP dim-red
#library("DropletUtils") #utility functions for handling single-cell (RNA-seq)
#library("AnnotationHub") # ensbl query
#library("AnnotationDbi") # ensbl query
#library("sctransform") # sc normalization,
#library("SingleR")
#library("velocyto.R")
#library(Matrix)
#library(stringr)
#library(celldex)
#library(scRNAseq)
#library(CellSIUS)
#require(RColorBrewer)
#library(ggplot2)
#library(rowr)
#library(reticulate)
#library(cowplot)
#library(clues)
#library(DT)
#library(scone)
#library(Signac)
#library(ape)
#library("tidyverse") # for data manipulation & plots
#library("lubridate") # right color
#library("BUSpaRse")
#library("TENxBUSData")
#library("simpleSingleCell") # for scRNA storage & manipulation
#library("scater") # for QC control
#library("scran") # analysis pipeline
#library(ggplot2)
#library(gridExtra)

#meta_variable <- c("seurat_clusters", "HTO_maxID", "Greffe", "SingleR.calls", "clonotype_id","chain", "v_gene", "d_gene", "j_gene","c_gene", "cdr3", "Phase")
#meta_variable <- c("seurat_clusters", "Condition", "Phénotype", "Phase", "K29Q", "L37M", "M11I", "pGln45")
#annotations <- read.csv("/home/boris/Bureau/scShiny/annotation_FindAllMarkers.csv")
#hallmark = all@tools$hallmarks
