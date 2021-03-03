######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## LIBRARY                                                                                                                       ########
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
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

load(file = "/home/boris/Documents/analyse/singlet_hFL_180008B.RData")