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
load(file = "/home/boris/Documents/analyse/singlet_FL08G0293.RData")

feature <- row.names(as.matrix(all[["RNA"]]@counts)); metadata <- c() ; List <- list()
for(i in 1:length(colnames(all@meta.data))){if (length(levels(as.factor(all@meta.data[[i]]))) > 1 && length(levels(as.factor(all@meta.data[[i]]))) < 600 && is.numeric(levels(as.factor(all@meta.data[[1]])))==F ){metadata <- c(metadata, colnames(all@meta.data)[i])}}
for(i in 1:length(metadata) ){List[[metadata[i]]] <- levels(as.factor(all@meta.data[[metadata[i]]]))}
all@tools$meta_variable <- c("seurat_clusters","Condition","Greffe","Phénotype","clonotype_id","Phase")   #c("dd","a") #si objet namnam
singlet <- all

