#Initial analysis of MCF-7 vs. Normal Breast data

#__________________________________________________________________________________________________________________________________________________________________________________________________
integrated.data <- readRDS(file = "C:/Users/jainar/Documents/integrated_data_NormalizeData_BeforeFindClusters")

integrated.data <- FindClusters(object = integrated.data, resolution = 0.5, dims.use = 1:40)

integrated.data <- RunUMAP(object = integrated.data, reduction = "pca", dims = 1:40)


#Find shared variable features between normal and MCF7 data and take out cell cycle genes

variable_features_normal <- integrated.data@assays$integrated@var.features

variable_features_MCF7 <- read.csv("C:/Users/jainar/Documents/MCF7variable_features.tsv", header = FALSE, sep = ",")

variable_features_MCF7 <- as.vector(variable_features_MCF7$V1)

shared_markers <- intersect(variable_features_MCF7, variable_features_normal)

cell_cycle_genes <- read.table("C:/Users/jainar/Documents/regev_lab_cell_cycle_genes.txt", header = FALSE)

cell_cycle_genes <- as.vector(cell_cycle_genes$V1)

shared_markers_nocellcycle <- setdiff(shared_markers, cell_cycle_genes)

# find markers for every cluster compared to all remaining cells, report only the positive ones
integrated.markers <- FindAllMarkers(integrated.data, only.pos = TRUE, logfc.threshold = 0.15)

#________________________________________________________________________________________________________
saveRDS(shared_markers_nocellcycle, "C:/Users/jainar/Documents/shared_markers_nocellcycle")

shared_markers_nocellcycle <- readRDS("C:/Users/jainar/Documents/shared_markers_nocellcycle") 

saveRDS(integrated.markers, "C:/Users/jainar/Documents/integrated.markers.res0.5")

integrated.markers <- readRDS("C:/Users/jainar/Documents/integrated.markers.res0.5")


#find cluster designations for each marker 

geneclus_list <- NULL
for (gene in shared_markers_nocellcycle) {
  if (gene %in% integrated.markers$gene) {
    index <- which(integrated.markers$gene == gene)
    gene_clusters <- gene
    for (i in index) {
      gene_clusters <- c(gene_clusters, as.numeric(integrated.markers[i,]$cluster)-1)
    }
    geneclus_list <- c(geneclus_list, list(gene_clusters))
  }
}

lapply(geneclus_list, function(x) write.table(data.frame(x), 'C:/Users/jainar/Documents/shared_markers_cluster_designations.tsv', append= T, sep='\t' , row.names = F, col.names = F))

#__________________________________________________________________________________________________________________________________________________
#Read in excel files containing excel file that assigns shared markers to MCF-7 clusters and normal breast clusters

#Determine number of markers each MCF7 cluster shares with normal breast clusters 


MCF7_cluster_genes <- read.csv(file = "C:/Users/jainar/Documents/MCF-7 genes for each cluster.csv", header = TRUE, sep = ",")

Normal_cluster_genes <- read.csv(file = "C:/Users/jainar/Documents/Normal breast - genes for each cluster.csv", header = TRUE, sep = ",")

common_clusters <- matrix(nrow = 15, ncol = 6)

rownames(common_clusters) <- colnames(Normal_cluster_genes)

colnames(common_clusters) <- colnames(MCF7_cluster_genes)

for(cluster in 0:14) {
  cluster_genes <- as.vector(Normal_cluster_genes[,cluster+1])
  cluster_genes <- cluster_genes[cluster_genes != ""]
  for(cluster2 in 0:5) {
    cluster_genes2 <- as.vector(MCF7_cluster_genes[, cluster2 + 1])
    cluster_genes2 <- cluster_genes2[cluster_genes2 != ""]
    num_shared <- as.numeric(length(intersect(cluster_genes, cluster_genes2)))
    common_clusters[cluster + 1, cluster2 + 1] <- num_shared
  }
}

write.csv(common_clusters, file = "C:/Users/jainar/Documents/corresponding_clusters.csv", row.names = TRUE, col.names = TRUE)

#___________________________________________________________________
# Correlation between MCF7 data and Normal Breast data (making heatmap)

MCF7 <- readRDS("C:/Users/jainar/Documents/cc_reduced_time.rds")

DimPlot(MCF7)

Idents(object = MCF7) <- 'RNA_snn_res.0.3'

healthy <- readRDS("C:/Users/jainar/Documents/integrated.data.res0.5")

variable_features_MCF7 <- read.csv("C:/Users/jainar/Documents/variable_features.tsv", header = FALSE, sep = ",")

variable_features_MCF7 <- as.vector(variable_features_MCF7$V1)

variable_features_healthy <- healthy@assays$integrated@var.features

combined_variable_features <- c(variable_features_healthy, variable_features_MCF7)

unique_variable_features <- unique(combined_variable_features)

features <- healthy@assays$RNA@data@Dimnames[[1]]
intersect_healthy <- intersect(features, unique_variable_features)

features2 <- MCF7@assays$RNA@data@Dimnames[[1]]
intersect_MCF7 <- intersect(features2, unique_variable_features)

genes.use <- intersect(intersect_healthy, intersect_MCF7)


MCF7_df <- matrix(nrow = length(genes.use), ncol = 6)
rownames(MCF7_df) <- genes.use
colnames(MCF7_df) <- paste("MCF7cluster.", 0:5, sep = "")

healthy_df <- matrix(nrow = length(genes.use), ncol = 15)
rownames(healthy_df) <- genes.use
colnames(healthy_df) <- paste("healthy.cluster", 0:14, sep="")

for(cluster in 0:14) {
  cells.use <- WhichCells(object = healthy, ident = cluster)
  expr <- GetAssayData(object = healthy, assay = "RNA", slot = "data")[genes.use, cells.use]
  values <- as.data.frame(rowSums(as.matrix(expr)))[,1]
  healthy_df[, cluster+1] <- values
}


for(cluster in 0:5) {
  cells.use <- WhichCells(object = MCF7, ident = cluster)
  expr <- GetAssayData(object = MCF7, assay = "RNA", slot = "data")[genes.use, cells.use]
  values <- as.data.frame(rowSums(as.matrix(expr)))[,1]
  MCF7_df[, cluster+1] <- values
}

saveRDS(MCF7_df, "C:/Users/jainar/Documents/MCF7_df")
saveRDS(healthy_df, "C:/Users/jainar/Documents/healthy_df")

#____________________________________________________________________________________________________
#correlation heatmap

combined <- cbind(MCF7_df, healthy_df)
cor.mtx <- cor(combined, method = "spearman")
cor.mtx2 <- cor(x = MCF7_df, y = healthy_df, method = "spearman")
saveRDS(cor.mtx2, "C:/Users/jainar/Documents/cor.mtx")
saveRDS(cor.mtx, "C:/Users/jainar/Documents/Arjun High School/High School Science Research/OHSU Internship - CEDAR/2. Normal Breast vs. MCF7 Analysis/1. Initial analysis/cor.mtx.rds")
pheatmap(cor.mtx2)
library(RColorBrewer)

cor(MCF7_df, method = "spearman")

readRDS(cor.)

# 
# #___________________________________________________________________________________________________________________
# 
# #Myoepithelial signature in healthy cells
# 
# DoHeatmap(healthy, c("TAGLN", "MYLK", "ACTA2", "TP63", "ACTG2", "VIM", "CALD1", "KRT5", "MYL9"))
# 
# cluster.averages <- AverageExpression(healthy, assays = "RNA", features = c("TAGLN", "MYLK", "ACTA2", "ACTG2", "VIM", "CALD1", "KRT5", "MYL9"))
# 
# cluster.averages <- cluster.averages$RNA
# 
# average_expression_per_cluster <- colMeans(cluster.averages)
# 
# average_expression_per_cluster

#_______________________________________________________________________________________________________________






