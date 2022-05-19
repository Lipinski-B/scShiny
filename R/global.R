#' Global parameters
#'
#'
#'

library(plotly)
library(shinyWidgets)
library(shinydashboard)

load(file = "/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_meta/All_Post-greffe_Other.RData")
Seurat::Idents(singlet)<-'Phénotype.fine'
singlet <- subset(singlet, idents = names(table(singlet@meta.data$Phénotype.fine)[which(table(singlet@meta.data$Phénotype.fine) > 10)]))
Seurat::Idents(singlet)<-'Sample'