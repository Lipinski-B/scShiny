library(celldex)
library(SingleR)
library(escape)
library(stringr)
library(dplyr)
library(future)
library(DESeq2)
library(cowplot)
library(Seurat)
library(dittoSeq)

library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE

source(file = "/home/boris/Bureau/scShiny/functions.R")

## -- Loading -- ## 
setwd(dir = "/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/")
siege <- c("FL09C1164","FL02G095","FL05G0330","FL140304", "FL08G0293", "FL12C1888") #'all'
patient <- siege[5]
load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))

################################################################################################################################################################################################################################################################################################################################################
## -- Workflow -- ##
singlet <- processing(patient) ; diet(singlet)
#save(singlet, file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
#save(singlet, file = paste0("/home/boris/Documents/lipinskib/flinovo/result/",patient,"/R/singlet_", patient,".RData"))
rm(singlet) ; gc() ; gc() ; gc()


## -- Sub sample -- ##
patient <- "FL05G0330"
load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
singlet <- seurat_subset(singlet,"Condition", c("RCHOP","Pré-greffe")) ; save(singlet, file = paste0("/home/boris/Documents/analyse/singlet_", patient,"_PG_RCHOP.RData"))
singlet <- seurat_subset(singlet,"Condition", c("Excipient","Pré-greffe")) ; save(singlet, file = paste0("/home/boris/Documents/analyse/singlet_", patient,"_PG_EX.RData"))
singlet <- seurat_subset(singlet,"Condition", c("RCHOP","Excipient")) ; save(singlet, file = paste0("/home/boris/Documents/analyse/singlet_", patient,"_EX_RCHOP.RData"))


for (patient in siege) {
  load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
  Idents(singlet)<-"Condition"
  singlet@tools$KEGG <- DEenrichRPlot(singlet, assay ="RNA", ident.1 = "RCHOP", ident.2 = "Excipient", test.use = "negbinom" ,max.cells.per.ident = Inf, balanced=T, p.val.cutoff=0.05, return.gene.list=T, num.pathway = 15, enrich.database = "KEGG_2021_Human", max.genes = Inf, logfc.threshold = 0.1)
  singlet@tools$GO_Biological <- DEenrichRPlot(singlet, assay ="RNA", ident.1 = "RCHOP", ident.2 = "Excipient", test.use = "negbinom", max.cells.per.ident = Inf, balanced=T, p.val.cutoff=0.05, return.gene.list=T, num.pathway = 15, enrich.database = "GO_Biological_Process_2021", max.genes = Inf, logfc.threshold = 0.1)
  singlet@tools$GO_Cellular <- DEenrichRPlot(singlet, assay ="RNA", ident.1 = "RCHOP", ident.2 = "Excipient", test.use = "negbinom", max.cells.per.ident = Inf, balanced=T, p.val.cutoff=0.05, return.gene.list=T, num.pathway = 15, enrich.database = "GO_Cellular_Component_2021", max.genes = Inf, logfc.threshold = 0.1)
  singlet@tools$GO_Molecular <- DEenrichRPlot(singlet, assay ="RNA", ident.1 = "RCHOP", ident.2 = "Excipient", test.use = "negbinom", max.cells.per.ident = Inf, balanced=T, p.val.cutoff=0.05, return.gene.list=T, num.pathway = 15, enrich.database = "GO_Molecular_Function_2021", max.genes = Inf, logfc.threshold = 0.1)
  Idents(singlet)<-"seurat_clusters"
  diet(singlet)
}




################################################################################################################################################################################################################################################################################################################################################
## -- Méta patient -- ##
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[1],".RData")) ; P1 <- all ; rownames(P1@meta.data) <- paste0("FL14_",rownames(P1@meta.data)) ; P1 <- metadata_merge(P1)
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[2],".RData")) ; P2 <- all ; rownames(P2@meta.data) <- paste0("FL12_",rownames(P2@meta.data)) ; P2 <- metadata_merge(P2)
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[3],".RData")) ; P3 <- all ; rownames(P3@meta.data) <- paste0("FL09_",rownames(P3@meta.data)) ; P3 <- metadata_merge(P3)
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[4],".RData")) ; P4 <- all ; rownames(P4@meta.data) <- paste0("FL08_",rownames(P4@meta.data)) ; P4 <- metadata_merge(P4)
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[5],".RData")) ; P5 <- all ; rownames(P5@meta.data) <- paste0("FL02_",rownames(P5@meta.data)) ; P5 <- metadata_merge(P5)
load(file = paste0("/home/boris/Documents/analyse/singlet_", siege[6],".RData")) ; P6 <- all ; rownames(P6@meta.data) <- paste0("FL05_",rownames(P6@meta.data)) ; P6 <- metadata_merge(P6)

all <- merge(P1, y = c(P2,P3,P4,P5,P6), add.cell.ids = c("FL14","FL12","FL09","FL08","FL02","FL05"), project = "FL")
save(all, file = paste0("/home/boris/Documents/analyse/singlet_all.RData"))
load(file = paste0("/home/boris/Documents/analyse/singlet_all.RData"))
all <- visualisation(all)
all@tools$hallmarks <- c("ADIPOGENESIS","ALLOGRAFT_REJECTION","ANDROGEN_RESPONSE","ANGIOGENESIS","APICAL_JUNCTION","APICAL_SURFACE","APOPTOSIS","BILE_ACID_METABOLISM","CHOLESTEROL_HOMEOSTASIS","COAGULATION","COMPLEMENT","DNA_REPAIR","E2F_TARGETS","EPITHELIAL_MESENCHYMAL_TRANSITION","ESTROGEN_RESPONSE_EARLY","ESTROGEN_RESPONSE_LATE","FATTY_ACID_METABOLISM","G2M_CHECKPOINT","GLYCOLYSIS","HEDGEHOG_SIGNALING","HEME_METABOLISM","HYPOXIA","IL2_STAT5_SIGNALING","IL6_JAK_STAT3_SIGNALING","INFLAMMATORY_RESPONSE","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","KRAS_SIGNALING_DN","KRAS_SIGNALING_UP","MITOTIC_SPINDLE","MTORC1_SIGNALING","MYC_TARGETS_V1","MYC_TARGETS_V2","MYOGENESIS","NOTCH_SIGNALING","OXIDATIVE_PHOSPHORYLATION","P53_PATHWAY","PANCREAS_BETA_CELLS","PEROXISOME","PI3K_AKT_MTOR_SIGNALING","PROTEIN_SECRETION","REACTIVE_OXYGEN_SPECIES_PATHWAY","SPERMATOGENESIS","TGF_BETA_SIGNALING","TNFA_SIGNALING_VIA_NFKB","UNFOLDED_PROTEIN_RESPONSE","UV_RESPONSE_DN","UV_RESPONSE_UP","WNT_BETA_CATENIN_SIGNALING","XENOBIOTIC_METABOLISM")
all@tools$meta_variable <- c("seurat_clusters", "Condition", "Greffe", "Phénotype", "clonotype_id", "Phase", "old.ident", "BCL2_L23L","BCL2_K22K","CD79B_Y696H","EZH2_A682G_1","EZH2_A682G_2","EZH2_A692V_1","EZH2_A692V_2","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S","EZH2_A692V","EZH2_A682G")#, "BCL2_K22K", "BCL2_L23L", "CD79B_Y696H")# c("seurat_clusters", "Condition", "Phénotype", "Phase", "K29Q", "L37M", "M11I", "pGln45")
all <- numeric_merge(all)
save(all, file = paste0("/home/boris/Documents/analyse/singlet_all.RData"))

## -- TEST : DE appareillé : RCHOP/EXP + PG : RC/RP -- ##
## test 1 = Réponse : RP/RC -- ##
load(file = paste0("/home/boris/Documents/analyse/singlet_all.RData"))
all <- seurat_subset(all, "Condition", c("Pré-greffe"))
all@meta.data[all@meta.data$orig.ident == "FL140304", "Reponse"] <- "RP"
all@meta.data[all@meta.data$orig.ident == "FL12C1888", "Reponse"] <- "RP"
all@meta.data[all@meta.data$orig.ident == "FL08G0293", "Reponse"] <- "RP"

all@meta.data[all@meta.data$orig.ident == "FL09C1164", "Reponse"] <- "RC"
all@meta.data[all@meta.data$orig.ident == "FL02G095", "Reponse"] <- "RC"
all@meta.data[all@meta.data$orig.ident == "FL05G0330", "Reponse"] <- "RC"

Idents(all)<-"Reponse" ; all[["RNA"]]@counts <- as.matrix(all[["RNA"]]@counts)+1
all@tools$DE_R <- FindMarkers(all, slot = "counts", ident.1 = "RP", ident.2 = "RC", test.use = "DESeq2", max.cells.per.ident = 2000)
all[["RNA"]]@counts <- as.matrix(all[["RNA"]]@counts)-1 ; Idents(all)<-"seurat_clusters"
save(all, file = "/home/boris/Documents/analyse/singlet_all_PG.RData")

## test 2 = Réponse : RP/RC -- ##




################################################################################################################################################################################################################################################################################################################################################
## -- DE intersect -- ##
DE <- list()
for (patient in siege) {
  load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
  DE[[patient]]$DE_RE <- all@tools$DE_RE[which(all@tools$DE_RE$p_val_adj < 0.05),]
  DE[[patient]]$DE_PE <- all@tools$DE_PE[which(all@tools$DE_PE$p_val_adj < 0.05),]
}
liste = read.table("/home/boris/Bureau/dissoc.txt", header = T)
liste <- liste$x
X <- 500
DE_RE <- Reduce(intersect, list(rownames(DE[["FL140304"]]$DE_RE)[1:X],rownames(DE[["FL12C1888"]]$DE_RE)[1:X],rownames(DE[["FL09C1164"]]$DE_RE)[1:X],rownames(DE[["FL08G0293"]]$DE_RE)[1:X],rownames(DE[["FL02G095"]]$DE_RE)[1:X],rownames(DE[["FL05G0330"]]$DE_RE)[1:X]))
DE_PE <- Reduce(intersect, list(rownames(DE[["FL140304"]]$DE_PE)[1:X],rownames(DE[["FL12C1888"]]$DE_PE)[1:X],rownames(DE[["FL09C1164"]]$DE_PE)[1:X],rownames(DE[["FL08G0293"]]$DE_PE)[1:X],rownames(DE[["FL02G095"]]$DE_PE)[1:X],rownames(DE[["FL05G0330"]]$DE_PE)[1:X]))

DE_RE_candidat <- setdiff(DE_RE,liste)
DE_PE_candidat <- setdiff(DE_PE,liste)

DE_PE_FC <- data.frame(matrix(ncol = 40, nrow = 6))
rownames(DE_PE_FC) <- siege ; colnames(DE_PE_FC) <- DE_PE_candidat
for (patient in siege) {DE_PE_FC[patient,] <- DE[[patient]]$DE_PE[rownames(DE[[patient]]$DE_PE) %in% DE_PE_candidat,]$avg_log2FC}

write.table(DE_PE_FC, "/home/boris/Bureau/DE_PE_FC.txt")
write.table(DE_PE_candidat_FC, "/home/boris/Bureau/DE_PE_FC.txt")






################################################################################################################################################################################################################################################################################################################################################
## -- TEST : DEenrichRPlot -- ##
ER <- list()
for (patient in siege) {
  load(file = paste0("/home/boris/Bureau/scShiny/www/", patient,"/", patient,".RData"))
  ER[[patient]]$KEGG$pos <- singlet@tools$KEGG$pos[1]
  ER[[patient]]$GO_Biological$pos <- singlet@tools$GO_Biological$pos[1]
  ER[[patient]]$GO_Cellular$pos <- singlet@tools$GO_Cellular$pos[1]
  ER[[patient]]$GO_Molecular$pos <- singlet@tools$GO_Molecular$pos[1]
    
  ER[[patient]]$KEGG$neg <- singlet@tools$KEGG$neg[1]
  ER[[patient]]$GO_Biological$neg <- singlet@tools$GO_Biological$neg[1]
  ER[[patient]]$GO_Cellular$neg <- singlet@tools$GO_Cellular$neg[1]
  ER[[patient]]$GO_Molecular$neg <- singlet@tools$GO_Molecular$neg[1]
  
  rm(singlet) ; gc() ; gc() ; gc()
}

ER_result <- list()
ER_result[["Complete"]]$GO_B$UP <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Biological$pos[1], ER[["FL12C1888"]]$GO_Biological$pos[1], ER[["FL08G0293"]]$GO_Biological$pos[1])) ))>0, Reduce(intersect, list(ER[["FL140304"]]$GO_Biological$pos[1], ER[["FL12C1888"]]$GO_Biological$pos[1], ER[["FL08G0293"]]$GO_Biological$pos[1])), NA)
ER_result[["Complete"]]$GO_B$DOWN <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Biological$neg[1], ER[["FL12C1888"]]$GO_Biological$neg[1], ER[["FL08G0293"]]$GO_Biological$neg[1])) ))>0,Reduce(intersect, list(ER[["FL140304"]]$GO_Biological$neg[1], ER[["FL12C1888"]]$GO_Biological$neg[1], ER[["FL08G0293"]]$GO_Biological$neg[1])), NA)
ER_result[["Complete"]]$GO_M$UP <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Molecular$pos[1], ER[["FL12C1888"]]$GO_Molecular$pos[1], ER[["FL08G0293"]]$GO_Molecular$pos[1])) ))>0, Reduce(intersect, list(ER[["FL140304"]]$GO_Molecular$pos[1], ER[["FL12C1888"]]$GO_Molecular$pos[1], ER[["FL08G0293"]]$GO_Molecular$pos[1])) , NA)
ER_result[["Complete"]]$GO_M$DOWN <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Molecular$neg[1], ER[["FL12C1888"]]$GO_Molecular$neg[1], ER[["FL08G0293"]]$GO_Molecular$neg[1])) ))>0, Reduce(intersect, list(ER[["FL140304"]]$GO_Molecular$neg[1], ER[["FL12C1888"]]$GO_Molecular$neg[1], ER[["FL08G0293"]]$GO_Molecular$neg[1])), NA)
ER_result[["Complete"]]$GO_C$UP <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Cellular$pos[1], ER[["FL12C1888"]]$GO_Cellular$pos[1], ER[["FL08G0293"]]$GO_Cellular$pos[1])) ))>0, Reduce(intersect, list(ER[["FL140304"]]$GO_Cellular$pos[1], ER[["FL12C1888"]]$GO_Cellular$pos[1], ER[["FL08G0293"]]$GO_Cellular$pos[1])), NA)
ER_result[["Complete"]]$GO_C$DOWN <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Cellular$neg[1], ER[["FL12C1888"]]$GO_Cellular$neg[1], ER[["FL08G0293"]]$GO_Cellular$neg[1])) ))>0, Reduce(intersect, list(ER[["FL140304"]]$GO_Cellular$neg[1], ER[["FL12C1888"]]$GO_Cellular$neg[1], ER[["FL08G0293"]]$GO_Cellular$neg[1])), NA)
ER_result[["Complete"]]$KEGG$UP <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$KEGG$pos[1], ER[["FL12C1888"]]$KEGG$pos[1], ER[["FL08G0293"]]$KEGG$pos[1])) ))>0, Reduce(intersect, list(ER[["FL140304"]]$KEGG$pos[1], ER[["FL12C1888"]]$KEGG$pos[1], ER[["FL08G0293"]]$KEGG$pos[1])), NA)
ER_result[["Complete"]]$KEGG$DOWN <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$KEGG$neg[1], ER[["FL12C1888"]]$KEGG$neg[1], ER[["FL08G0293"]]$KEGG$neg[1])) ))>0, Reduce(intersect, list(ER[["FL140304"]]$KEGG$neg[1], ER[["FL12C1888"]]$KEGG$neg[1], ER[["FL08G0293"]]$KEGG$neg[1])), NA)


ER_result[["Partiel"]]$GO_B$UP <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Biological$pos[1], ER[["FL02G095"]]$GO_Biological$pos[1], ER[["FL05G0330"]]$GO_Biological$pos[1])) ))>0, Reduce(intersect, list(ER[["FL09C1164"]]$GO_Biological$pos[1], ER[["FL02G095"]]$GO_Biological$pos[1], ER[["FL05G0330"]]$GO_Biological$pos[1])), NA)
ER_result[["Partiel"]]$GO_B$DOWN <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Biological$neg[1], ER[["FL02G095"]]$GO_Biological$neg[1], ER[["FL05G0330"]]$GO_Biological$neg[1])) ))>0, Reduce(intersect, list(ER[["FL09C1164"]]$GO_Biological$neg[1], ER[["FL02G095"]]$GO_Biological$neg[1], ER[["FL05G0330"]]$GO_Biological$neg[1])), NA)
ER_result[["Partiel"]]$GO_M$UP <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Molecular$pos[1], ER[["FL02G095"]]$GO_Molecular$pos[1], ER[["FL05G0330"]]$GO_Molecular$pos[1])) ))>0, Reduce(intersect, list(ER[["FL09C1164"]]$GO_Molecular$pos[1], ER[["FL02G095"]]$GO_Molecular$pos[1], ER[["FL05G0330"]]$GO_Molecular$pos[1])), NA)
ER_result[["Partiel"]]$GO_M$DOWN <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Molecular$neg[1], ER[["FL02G095"]]$GO_Molecular$neg[1], ER[["FL05G0330"]]$GO_Molecular$neg[1])) ))>0, Reduce(intersect, list(ER[["FL09C1164"]]$GO_Molecular$neg[1], ER[["FL02G095"]]$GO_Molecular$neg[1], ER[["FL05G0330"]]$GO_Molecular$neg[1])) , NA)
ER_result[["Partiel"]]$GO_C$UP <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Cellular$pos[1], ER[["FL02G095"]]$GO_Cellular$pos[1], ER[["FL05G0330"]]$GO_Cellular$pos[1])) ))>0, Reduce(intersect, list(ER[["FL09C1164"]]$GO_Cellular$pos[1], ER[["FL02G095"]]$GO_Cellular$pos[1], ER[["FL05G0330"]]$GO_Cellular$pos[1])), NA)
ER_result[["Partiel"]]$GO_C$DOWN <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Cellular$neg[1], ER[["FL02G095"]]$GO_Cellular$neg[1], ER[["FL05G0330"]]$GO_Cellular$neg[1])) ))>0, Reduce(intersect, list(ER[["FL09C1164"]]$GO_Cellular$neg[1], ER[["FL02G095"]]$GO_Cellular$neg[1], ER[["FL05G0330"]]$GO_Cellular$neg[1])), NA)
ER_result[["Partiel"]]$KEGG$UP <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$KEGG$pos[1], ER[["FL02G095"]]$KEGG$pos[1], ER[["FL05G0330"]]$KEGG$pos[1])) ))>0, Reduce(intersect, list(ER[["FL09C1164"]]$KEGG$pos[1], ER[["FL02G095"]]$KEGG$pos[1], ER[["FL05G0330"]]$KEGG$pos[1])), NA)
ER_result[["Partiel"]]$KEGG$DOWN <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$KEGG$neg[1], ER[["FL02G095"]]$KEGG$neg[1], ER[["FL05G0330"]]$KEGG$neg[1])) ))>0, Reduce(intersect, list(ER[["FL09C1164"]]$KEGG$neg[1], ER[["FL02G095"]]$KEGG$neg[1], ER[["FL05G0330"]]$KEGG$neg[1])), NA)


df <- tibble(Groupe = c("Complete", "Complete", "Complete", "Complete", "Partiel", "Partiel", "Partiel", "Partiel"),
             Enrichissement = c("GO : Biological Process", "GO : Molecular Function", "GO : Cellular Component", "KEGG", "GO : Biological", "GO : Molecular", "GO : Cellular", "KEGG"),
             UP = c(paste(unlist(ER_result[["Complete"]]$GO_B$P), collapse = " ; "), paste(unlist(ER_result[["Complete"]]$GO_M$P), collapse = " ; "), paste(unlist(ER_result[["Complete"]]$GO_C$P), collapse = " ; "), paste(unlist(ER_result[["Complete"]]$KEGG$P), collapse = " ; "),
                    paste(unlist(ER_result[["Partiel"]]$GO_B$P), collapse = " ; "),  paste(unlist(ER_result[["Partiel"]]$GO_M$P), collapse = " ; "), paste(unlist(ER_result[["Partiel"]]$GO_C$P), collapse = " ; "), paste(unlist(ER_result[["Partiel"]]$KEGG$P), collapse = " ; ")),
             DOWN = c(paste(unlist(ER_result[["Complete"]]$GO_B$N), collapse = " ; "), paste(unlist(ER_result[["Complete"]]$GO_M$N), collapse = " ; "), paste(unlist(ER_result[["Complete"]]$GO_C$N), collapse = " ; "), paste(unlist(ER_result[["Complete"]]$KEGG$N), collapse = " ; "),
                      paste(unlist(ER_result[["Partiel"]]$GO_B$N), collapse = " ; "),  paste(unlist(ER_result[["Partiel"]]$GO_M$N), collapse = " ; "), paste(unlist(ER_result[["Partiel"]]$GO_C$N), collapse = " ; "), paste(unlist(ER_result[["Partiel"]]$KEGG$N), collapse = " ; ")))
      
df %>% datatable(rownames = FALSE)


met <- data.frame(
  Cas = c(rep("", length(unlist(ER_result)))),
  Enrichissement = c(rep("", length(unlist(ER_result)))),
  Mode = c(rep("", length(unlist(ER_result)))),
  Name = c(rep("", length(unlist(ER_result)))),
  value = rep(1, 35), stringsAsFactors = FALSE
)

i=1
for (cas in names(ER_result)) {
  for (enrichissement in names(ER_result[[cas]])) {
     for (mode in names(ER_result[[cas]][[enrichissement]])) {
       for (name in unlist(ER_result[[cas]][[enrichissement]][[mode]][])) {
         met[i,"Cas"]=cas
         met[i,"Enrichissement"]=enrichissement
         met[i,"Mode"]=mode
         met[i,"Name"]=name
         i=i+1
       }
     } 
  }
}


df <- data.frame(ids=character(), labels=character(), parents=character(), values=integer(),stringsAsFactors=FALSE)
df0 <- met %>% mutate(path = paste(Cas,  sep=";")) %>% dplyr::select(path, value)
df1 <- met %>% mutate(path = paste(Cas, Enrichissement,  sep=";")) %>% dplyr::select(path, value)
df2 <- met %>% mutate(path = paste(Cas, Enrichissement, Mode,  sep=";")) %>% dplyr::select(path, value)
df3 <- met %>% mutate(path = paste(Cas, Enrichissement, Mode, Name, sep=";")) %>% dplyr::select(path, value)

for (dfX in c(df0,df1,df2,df3)) {
  xnames <- row.names(table(dfX))
  max <- 0
  
  for (rows in xnames) {
    orga = strsplit(rows,";")
    if (max < length(orga[[1]])){max = length(orga[[1]])}}
  
  for (i in 1:max){
    for (rows in xnames) {
      orga = strsplit(rows,";")
      if ((length(orga[[1]]))==i){
        c = ""
        c2 = ""
        for (j in 1:i){c = paste0(c,orga[[1]][j],"-")}
        for (k in 1:i-1){c2 = paste0(c2,orga[[1]][k],"-")}
        c = str_sub(c,1,-2)
        c2 = str_sub(c2,2,-2)
        df[rows,"labels"]= orga[[1]][i]
        df[rows,"ids"]=c
        df[rows,"parents"]=c2
      }
    }
  }
}
df <- df[-3,]
df$values <- c(table(df0)[,1],table(df1)[,1],table(df2)[,1],table(df3)[,1])

plot_ly(df, ids = ~ids ,labels = ~labels, parents = ~parents, values = ~values, branchvalues = 'total', type = 'sunburst')




################################################################################################################################################################################################################################################################################################################################################
## -- TEST : DE RCHOP : apop+ / apop- -- ##
patient <- "FL05G0330"
load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
seuil = 0.13
liste = as.matrix(rownames(all@meta.data)[which(all@meta.data$APOPTOSIS<seuil & all@meta.data$Phénotype=="B-cells" & all@meta.data$Condition =="RCHOP")])
liste2 = as.matrix(rownames(all@meta.data)[which(all@meta.data$APOPTOSIS>seuil & all@meta.data$Phénotype=="B-cells" & all@meta.data$Condition =="RCHOP")])
singlet@meta.data[singlet@meta.data$Condition == "RCHOP", "Apop" ] <- "+"
singlet@meta.data[singlet@meta.data$Condition == "RCHOP", "Apop" ] <- "-"