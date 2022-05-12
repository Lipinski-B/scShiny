#' Global parameters
#'
#'
#'

library(plotly)
library(shinyWidgets)
library(shinydashboard)

#load(file = "datasets/Patient/All_All_FL08G0293.RData")
load(file = "/home/boris/Bureau/Flinovo/result/analyse_meta/All_Post-greffe_Other.RData")
#load(file = "/home/boris/Bureau/Flinovo/result/analyse_patient/FL08G0293/All_All_FL08G0293.RData")
#load(file = paste0("/home/boris/Bureau/Flinovo/result/analyse_meta/All_Post-greffe_All.RData"))

Seurat::Idents(singlet)<-'Phénotype.fine'
singlet <- subset(singlet, idents = names(table(singlet@meta.data$Phénotype.fine)[which(table(singlet@meta.data$Phénotype.fine) > 10)]))
Seurat::Idents(singlet)<-'Sample'