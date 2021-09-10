source(file = "functions.R")
load(file = "www/FL08G0293.RData")
#load(file = "/home/boris/Documents/analyse/FL08G0293.RData") 

feature <- row.names(as.matrix(all[["RNA"]]@data)); metadata <- c() ; List <- list()
for(i in 1:length(colnames(all@meta.data))){if (length(levels(as.factor(all@meta.data[[i]]))) > 1 && length(levels(as.factor(all@meta.data[[i]]))) < 600 && is.numeric(levels(as.factor(all@meta.data[[1]])))==F ){metadata <- c(metadata, colnames(all@meta.data)[i])}}
for(i in 1:length(metadata) ){List[[metadata[i]]] <- levels(as.factor(all@meta.data[[metadata[i]]]))}
singlet <- all