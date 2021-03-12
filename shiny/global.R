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

load(file = "/home/boris/Documents/analyse/singlet_hFL_180008B.RData")
meta_variable <- c("seurat_clusters", "HTO_maxID", "Greffe", "SingleR.calls", "clonotype_id", "Phase")
#meta_variable <- c("seurat_clusters", "HTO_maxID", "Greffe", "SingleR.calls", "clonotype_id","chain", "v_gene", "d_gene", "j_gene","c_gene", "cdr3", "Phase")
annotations <- read.csv("/home/boris/Bureau/scShiny/annotation_FindAllMarkers.csv")
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
hallmark = singlet@tools$hallmarks
