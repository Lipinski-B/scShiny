library(monocle)
library(ggplot2)
library(cowplot)
library(clues)
library(dplyr)

############################ Conversion entre Seurat V3 object et Monocle 2.0 ############################
importCDS <- function (otherCDS, seurat_scale=F, import_all = FALSE) {
  if (class(otherCDS)[1] == "Seurat") {
    requireNamespace("Seurat")
    if (!seurat_scale) {
      data <- otherCDS@assays$RNA@counts
    } else {
      data <- otherCDS@assays$RNA@scale.data
    }
    if (class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }
    pd <- tryCatch({
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, error = function(e) {
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)
      message("This Seurat object doesn't provide any meta data")
      pd
    })
    if (length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]
    }
    fData <- data.frame(gene_short_name = row.names(data), 
                        row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    #lowerDetectionLimit <- otherCDS@is.expr
    if (all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
      expr <- "negbinomial.size"
    }
    else if (any(data < 0)) {
      expressionFamily <- uninormal()
      expr <- "unimormal"
    }
    else {
      expressionFamily <- tobit()
      expr <- "tobit"
    }
    print(paste0("expressionFamily ",expr))
    # valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                                  #lowerDetectionLimit = lowerDetectionLimit,
                                  expressionFamily = expressionFamily)
    if (import_all) {
      if ("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
      }
      else {
        mist_list <- otherCDS
      }
    }
    else {
      mist_list <- list()
    }
    if ("var.genes" %in% slotNames(otherCDS)) {
      var.genes <- setOrderingFilter(monocle_cds, otherCDS@var.genes)
    }
    monocle_cds@auxClusteringData$seurat <- mist_list
  }
  else if (class(otherCDS)[1] == "SCESet") {
    requireNamespace("scater")
    message("Converting the exprs data in log scale back to original scale ...")
    data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
    fd <- otherCDS@featureData
    pd <- otherCDS@phenoData
    experimentData = otherCDS@experimentData
    if ("is.expr" %in% slotNames(otherCDS)) 
      lowerDetectionLimit <- otherCDS@is.expr
    else lowerDetectionLimit <- 1
    if (all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    }
    else if (any(data < 0)) {
      expressionFamily <- uninormal()
    }
    else {
      expressionFamily <- tobit()
    }
    if (import_all) {
      mist_list <- otherCDS
    }
    else {
      mist_list <- list()
    }
    monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                                  lowerDetectionLimit = lowerDetectionLimit, expressionFamily = expressionFamily)
    monocle_cds@auxOrderingData$scran <- mist_list
  }
  else {
    stop("the object type you want to export to is not supported yet")
  }
  return(monocle_cds)
} # Importer object seurat (sous forme de matrix) dans monocle: CellDataSet (CDS)
data <- importCDS(singlet)
dim(pData(data))                          # Check out the phenotypic data (nombre de cells + nombre metadatas de l'object CDS)
head(pData(data))

data <- estimateSizeFactors(data)         # Normalization 
data <- estimateDispersions(data)
data <- detectGenes(data, min_expr = 0.1) # Filtering low-quality cells (gêne exprimé si min 1 compte: min_expr = 0.1 ) 

head(fData(data))                         # Visualisation: nombre de cells qui expriment un gêne spé: dans fData
summary(fData(data)$num_cells_expressed)
sum((exprs(data['CD19',])))               # Si on veut regarder expr d'un gene de façon individuelle

head(pData(data))                         # Visualisation: nombre de gênes exprimés dans une cell: dans pData (stored in phenoData)
summary(pData(data)$num_genes_expressed)

############################ Standardise to Z-distribution ############################
x <- pData(data)$num_genes_expressed
x_1 <- (x - mean(x)) / sd(x)
summary(x_1)
ggplot(data.frame(x = x_1), aes(x)) +geom_histogram(bins = 50) +geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')

pData(data)$UMI <- Matrix::colSums(exprs(data)) # add UMI
head(pData(data))
ggplot(pData(data), aes(num_genes_expressed, UMI)) + geom_point() #UMI / gene


############################ Clustering cells WITHOUT markers genes############################
######### 1ere façon #########
disp_table <- dispersionTable(data)     # dispersionTable() function calculates the mean and dispersion values 
head(disp_table, n= 10)
table(disp_table$mean_expression>=0.1)  # Select genes, which have a mean expression >= 0.1, to use in the clustering step

unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
data <- setOrderingFilter(data, unsup_clustering_genes$gene_id)
plot_ordering_genes(data)

######### 2eme façon #########
plot_pc_variance_explained(data, max_components = 50, return_all = FALSE) # Elbow Plot (long)

## Include X dimensions: tSNE
data <- reduceDimension(data, max_components = 2, reduction_method = 'tSNE', verbose = TRUE)
data <- clusterCells(data)
plot_cell_clusters(data)
my_cluster_dim <- pData(data)$Cluster


data <- reduceDimension(data, max_components = 2, num_dim = 10, reduction_method = 'tSNE', verbose = TRUE)
data <- clusterCells(data, num_clusters = 15)
plot_cell_clusters(data)
my_cluster_dim_10 <- pData(data)$Cluster

adjustedRand(as.numeric(my_cluster_dim_10), as.numeric(my_cluster_dim))


############################## Expression differentielle ##############################
## I’ll create a vector that indicates whether a single was clustered in “Cluster 1” or not to identify genes differentially expressed in cluster 1 versus the other clusters.
my_vector <- rep('no', nrow(pData(data))) ## Create vector of no's
my_vector[pData(data)$Cluster == 1] <- rep('yes', sum(pData(data)$Cluster == 1)) ## Change status to yes if the cell was in cluster 1
pData(data)$test <- my_vector ## Add vector to phenoData
head(pData(data))

## Differential expression analysis on the subset of genes where the mean expression (across cells) was >= 0.1
length(unsup_clustering_genes$gene_id)
de_cluster_one <- differentialGeneTest(data[unsup_clustering_genes$gene_id,],fullModelFormulaStr = '~test',cores = 8)
dim(de_cluster_one)

de_cluster_one %>% arrange(qval) %>% head() ## Order by q-value
plot_genes_jitter(data['SERTAD1',], grouping = "Cluster")

## Add another column to the phenotypic data where TRUE means cluster membership in clusters 1, 2, 7 or 10
pData(data)$my_colour <- pData(data)$Cluster == 1 | pData(data)$Cluster == 2 | pData(data)$Cluster == 7 | pData(data)$Cluster == 10
plot_cell_clusters(data, color_by = 'my_colour')


############################## Classifying cells by type ##############################
## Identifier gênes intérêt
MS4A1_id <- row.names(subset(fData(data), gene_short_name == "MS4A1"))
CXCR4_id <- row.names(subset(fData(data),gene_short_name == "CXCR4"))
CD83_id <- row.names(subset(fData(data),gene_short_name == "CD83"))

cth <- newCellTypeHierarchy() ## Ajouter type cellulaire en fonction expression gênes intéret 
cth <- addCellType(cth, "LZ cells", classify_func =function(x) { x[MS4A1_id,] >= 1 & x[CXCR4_id,] < 1 & x[CD83_id,] > 1})
cth <- addCellType(cth, "DZ cells", classify_func = function(x){ x[MS4A1_id,] >= 1 & x[CXCR4_id,] > 1 & x[CD83_id,] < 1})

data <- classifyCells(data, cth, 0.1) ## Classification des cellules
table(pData(data)$CellType)

## Camembert types cellulaires
pie <- ggplot(pData(data),aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

plot_cell_clusters(data, 1, 2, color = "CellType") ## Représentation tSNE fonction type cellulaire
plot_cell_clusters(data, 1, 2, color = "CellType", markers = "CD19") ## Représentation tSNE fonction gênes types
plot_cell_clusters(data, 1, 2, color = "Cluster") +facet_wrap(~CellType) ## Représentation tSNE fonction types cellulaires + clusters


############################## Constructing Single Cell Trajectories ##############################
## Trajectory step 1: choose genes that define a cell's progress
expressed_genes <- row.names(subset(fData(data), num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(data[expressed_genes,],fullModelFormulaStr = "~CellType")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
data <- setOrderingFilter(data, ordering_genes)
plot_ordering_genes(data)

## Trajectory step 2: reduce data dimensionality (Long)
data <- reduceDimension(data, max_components = 2,method = 'DDRTree')

## Trajectory step 3: order cells along the trajectory
data <- orderCells(data)
plot_cell_trajectory(data, color_by = "CellType")
plot_cell_trajectory(data, color_by = "State")

### Ne fonctionne pas ??? 
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$CellType)[,"0"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))}
  else{return(1)}
}

### Ne fonctionne pas 
data <- orderCells(data, root_state = GM_state(data))

### (Suite)
plot_cell_trajectory(data, color_by = "Pseudotime")
plot_cell_trajectory(data, color_by = "State") + facet_wrap(~State, nrow = 1)

data_filtered <- data[expressed_genes,]
my_genes <- row.names(subset(fData(data_filtered),gene_short_name %in% "CD3E"))
cds_subset <- data_filtered[my_genes,]

plot_genes_in_pseudotime(cds_subset, color_by = "CellType")


















### A partir de là c'est plus la partie du blog où il clusterise lui en fontion de cluster
### Je me suis pas vraiment penchée dessus, parce que je sais pas trop si c'est pertinent
############################## Constructing Single Cell Trajectories ##############################
### For the trajectory analysis in this post, I will use only a subset of all cells (cluster 1) and genes that are expressed in at least 10 cells.
expressed_genes <- row.names(subset(fData(data), num_cells_expressed >= 10))

# my_colour is TRUE if a cell belongs to either cluster 3, 4, 7 or 10
data_subset <- data[expressed_genes, pData(data)$my_colour]
data_subset


# keeping only genes expressed in greater than 5% of cells
data_subset <- detectGenes(data_subset, min_expr = 0.1)
fData(data_subset)$use_for_ordering <- fData(data_subset)$num_cells_expressed > 0.05 * ncol(data_subset)

# how many genes are used?
table(fData(data_subset)$use_for_ordering)

plot_pc_variance_explained(data_subset, return_all = FALSE)

data_subset <- reduceDimension(data_subset, max_components = 2,norm_method = 'log',num_dim = 10,reduction_method = 'tSNE',verbose = TRUE)
data_subset <- clusterCells(data_subset, verbose = FALSE)
plot_rho_delta(data_subset, rho_threshold = 2, delta_threshold = 10)

data_subset <- clusterCells(data_subset,rho_threshold = 2,delta_threshold = 10,skip_rho_sigma = T,verbose = FALSE)
table(pData(data_subset)$Cluster)
plot_cell_clusters(data_subset)

clustering_DEG_genes <- differentialGeneTest(data_subset,fullModelFormulaStr = '~Cluster',cores = 8)
dim(clustering_DEG_genes)

clustering_DEG_genes %>% arrange(qval) %>% head()

my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
data_subset <- setOrderingFilter(data_subset, ordering_genes = my_ordering_genes)
data_subset <- reduceDimension(data_subset, method = 'DDRTree')

# the warnings were for use of deprecated code
data_subset <- orderCells(data_subset)

plot_cell_trajectory(data_subset, color_by = "Cluster")

# pseudotime is now a column in the phenotypic data as well as the cell state
head(pData(data_subset))

my_pseudotime_de <- differentialGeneTest(data_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 8)
my_pseudotime_de %>% arrange(qval) %>% head()

# save the top 6 genes
my_pseudotime_de %>% arrange(qval) %>% head(6) %>% select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene <- my_pseudotime_gene$gene_short_name

plot_genes_in_pseudotime(data_subset[my_pseudotime_gene,])

# cluster the top 50 genes that vary as a function of pseudotime
my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(status) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster$status

my_pseudotime_cluster <- plot_pseudotime_heatmap(data_subset[gene_to_cluster,], num_clusters = 3, cores = 8, show_rownames = TRUE, return_heatmap = TRUE)