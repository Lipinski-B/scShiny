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
library(ggplot2)
library(gridExtra)
library(shinydashboard)
library(Signac)
library(ape)
library(plotly)
library(scone)
suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(dittoSeq))
library(shinyWidgets)
library(monocle)
library(cowplot)
library(clues)
library(dplyr)
library(DT)
library(dashboardthemes)
library(reticulate)

load(file = paste0("/home/boris/Documents/analyse/singlet_FL12C1888.RData"))

meta_variable <- c("seurat_clusters", "HTO_maxID", "Greffe", "SingleR.calls", "clonotype_id", "Phase")
#meta_variable <- c("seurat_clusters", "HTO_maxID", "Greffe", "SingleR.calls", "clonotype_id","chain", "v_gene", "d_gene", "j_gene","c_gene", "cdr3", "Phase")
annotations <- read.csv("/home/boris/Bureau/scShiny/annotation_FindAllMarkers.csv")
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
hallmark = all@tools$hallmarks

metadata <- c()
for(i in 1:length(colnames(all@meta.data))){
if (length(levels(as.factor(all@meta.data[[i]]))) > 1 && length(levels(as.factor(all@meta.data[[i]]))) < 25 && is.numeric(levels(as.factor(all@meta.data[[1]])))==F ){
    metadata <- c(metadata, colnames(all@meta.data)[i])
  }
}

List <- list()
for(i in 1:length(metadata) ){List[[metadata[i]]] <- levels(as.factor(all@meta.data[[metadata[i]]]))}

singlet <- all