suppressMessages(library(scRepertoire))
suppressMessages(library(Seurat))
setwd(dir = "/home/boris/Documents/lipinskib/flinovo/result/")

load(file = "/home/boris/Documents/analyse/singlet_FL12C1888.RData")

patient1 <- "FL12C1888"
csv1 <- list(read.csv(paste0(patient1,"/VDJ/",patient1,"_CellrangerVDJ/outs/filtered_contig_annotations.csv"), stringsAsFactors = FALSE))
patient2 <- "FL140304"
csv2 <- list(read.csv(paste0(patient2,"/VDJ/",patient2,"_CellrangerVDJ/outs/filtered_contig_annotations.csv"), stringsAsFactors = FALSE))
patient3 <- "hFL_180008B"
csv3 <- list(read.csv("/home/boris/Documents/lipinskib/flinovo/result/hFL_180008B/hFL_180008B_CellrangerVDJ/outs/filtered_contig_annotations.csv", stringsAsFactors = FALSE))
patient4 <- "hFL_130337"
csv4 <- list(read.csv("/home/boris/Documents/lipinskib/boris/Cellranger/result/hFL_130337/VDJ/hFL_130337_CellrangerVDJ/outs/filtered_contig_annotations.csv", stringsAsFactors = FALSE))

combined <- combineBCR(c(csv1,csv2,csv3,csv4), samples = c(patient1, patient2, patient3, patient4), ID = c("BCR","BCR","BCR","BCR"))

quantContig(combined, cloneCall="gene+nt", scale = TRUE)
quantContig_output <- quantContig(combined, cloneCall="gene+nt", scale = TRUE, exportTable = TRUE)
quantContig_output

abundanceContig(combined, cloneCall = "gene", scale = F)
abundanceContig(combined, cloneCall = "gene", exportTable = T)
lengthContig(combined, cloneCall="aa", chains = "combined") 
lengthContig(combined, cloneCall="nt", chains = "single") 
compareClonotypes(combined, numbers = 10, samples = c("FL12C1888_BCR", "FL140304_BCR","hFL_180008B_BCR","hFL_130337_BCR"),cloneCall="aa", graph = "alluvial")
vizVgenes(combined, TCR="TRA", facet.x = "sample", facet.y = "ID")
clonalHomeostasis(combined, cloneCall = "gene")
clonalHomeostasis(combined, cloneCall = "aa")
clonalProportion(combined, cloneCall = "gene") 
clonalProportion(combined, cloneCall = "nt") 
clonalOverlap(combined, cloneCall = "gene+nt", method = "morisita")
clonesizeDistribution(combined, cloneCall = "gene+nt", method="ward.D2")
clonalDiversity(combined, cloneCall = "gene", group = "samples")#,n.boots = 100)
clonalDiversity(combined, cloneCall = "gene", group = "ID")



combined <- combineBCR(c(csv1,csv2), samples = c(patient1, patient2), ID = c("BCR","BCR"))
load(file = "/home/boris/Documents/analyse/singlet_FULL.RData")

seurat <- combineExpression(combined, all, 
                            cloneCall="gene+nt",# groupBy = "sample",# proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500), filterNA = T)


colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", 
                                        "#C6FDEC", "#7AC5FF", "#0348A6"))
DimPlot(all, group.by = "Type") + NoLegend() +
  scale_color_manual(values=colorblind_vector(2))


DimPlot(all, label = T) + NoLegend()
table(Idents(all))








library(stringr)
combined <- combineBCR(csv2, samples = patient2, ID = "BCR")
load(file = "/home/boris/Documents/analyse/singlet_FL140304.RData")

combined$FL140304_BCR$barcode <- str_replace(combined$FL140304_BCR$barcode, "FL140304_BCR_", "")
seurat <- combineExpression(combined, all,cloneCall="gene+nt", groupBy = "sample",# proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500), filterNA = T)

DimPlot(all, label = T, reduction = "pca") + NoLegend()
DimPlot(seurat, group.by = "cloneType",reduction = "pca") + NoLegend() + scale_color_manual(values= colorRampPalette(c("#FF4B20", "#FFB433","#C6FDEC", "#7AC5FF", "#0348A6"))(4))

slot(seurat, "meta.data")$cloneType <- factor(slot(seurat, "meta.data")$cloneType, 
                                              levels = c("Hyperexpanded (100 < X <= 500)", 
                                                         "Large (20 < X <= 100)", 
                                                         "Medium (5 < X <= 20)", 
                                                         "Small (1 < X <= 5)", 
                                                         "Single (0 < X <= 1)", NA))
DimPlot(seurat, group.by = "cloneType", reduction = "pca") + scale_color_manual(values = colorblind_vector(5), na.value="grey")

clonalOverlay(seurat, reduction = "umap",freq.cutpoint = 30, bins = 10, facet = "Patient") + guides(color = FALSE)

