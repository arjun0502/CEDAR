# Arjun Jain
# 10 July 2019
# Analysis of Normal Breast Epithelium 10x scRNA-seq Data with NormalizeData() and not SCTransform()
  # Data: 4 raw counts matrices for 4 individuals (ind4-7)
  #~24000 cells

#Load Seurat

library(Seurat)

# Download counts matrices
ind4_counts <- read.table(file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/data/GSM3099846_Ind4_Expression_Matrix.txt", sep = "\t", header = T)
ind5_counts <- read.table(file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/data/GSM3099847_Ind5_Expression_Matrix.txt", sep = "\t", header = T)
ind6_counts <- read.table(file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/data/GSM3099848_Ind6_Expression_Matrix.txt", sep = "\t", header = T)
ind7_counts <- read.table(file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/data/GSM3099849_Ind7_Expression_Matrix.txt", sep = "\t", header = T)

# For each individual count matrices, make the genes the rowames and delete the first column
rownames(ind4_counts) <- ind4_counts$X
ind4_counts <- ind4_counts[,-1]

rownames(ind5_counts) <- ind5_counts$X
ind5_counts <- ind5_counts[,-1]

rownames(ind6_counts) <- ind6_counts$X
ind6_counts <- ind6_counts[,-1]

rownames(ind7_counts) <- ind7_counts$X
ind7_counts <- ind7_counts[,-1]

# Determine dimensions (no. of genes x no. of cells) for each dataset
#dim(ind4_counts)  #33694 genes and 4116 cells
#dim(ind5_counts)  #33694 genes and 6964 cells
#dim(ind6_counts) #33694 genes and 6014 cells
#dim(ind7_counts) #33694 genes and 7552 cells

# Create Seurat objects for each individual data

# Filter data when creating seurat objects
# only include features detected in at least 5 cells
# only include cells where at least 200 features are detected

ind4_seurat <- CreateSeuratObject(counts = ind4_counts, min.cells = 5, min.features = 200)

ind5_seurat <- CreateSeuratObject(counts = ind5_counts, min.cells = 5, min.features = 200)

ind6_seurat <- CreateSeuratObject(counts = ind6_counts, min.cells = 5, min.features = 200)

ind7_seurat <- CreateSeuratObject(counts = ind7_counts, min.cells = 5, min.features = 200)
#_____________________________________________________________________________________________________________________________________
# QUALITY CONTROL

## Store mitochondrial percentage in meta data
ind4_seurat[["percent.mt"]] <- PercentageFeatureSet(ind4_seurat, pattern = "^MT-") 
ind5_seurat[["percent.mt"]] <- PercentageFeatureSet(ind5_seurat, pattern = "^MT-") 
ind6_seurat[["percent.mt"]] <- PercentageFeatureSet(ind6_seurat, pattern = "^MT-") 
ind7_seurat[["percent.mt"]] <- PercentageFeatureSet(ind7_seurat, pattern = "^MT-") 

## Visualize QC metrics

VlnPlot(ind4_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2)
plot1 <- FeatureScatter(ind4_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ind4_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(ind5_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(ind5_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ind5_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(ind6_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(ind6_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ind6_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(ind7_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(ind7_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ind7_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

## subset counts matrices based on violin and feature plots
ind4_seurat <- subset(ind4_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 13)
ind5_seurat <- subset(ind5_seurat, subset = nFeature_RNA > 800 & nFeature_RNA < 5477 & percent.mt < 7)
ind6_seurat <- subset(ind6_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 5)
ind7_seurat <- subset(ind7_seurat, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 8)

## Check dimensions ofsubsetted Seurat objects
dim(ind4_seurat) #16067 genes and 4017 cells
dim(ind5_seurat) #18969 genes and 6835 cells
dim(ind6_seurat) #17890 genes and 5889 cells
dim(ind7_seurat) #18298 genes and 7477 cells

## Save seurat objects as RDS files
saveRDS(ind4_seurat, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/individual seurat objects with QC/ind4_seurat_withQC.rds")
saveRDS(ind5_seurat, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/individual seurat objects with QC/ind5_seurat_withQC.rds")
saveRDS(ind6_seurat, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/individual seurat objects with QC/ind6_seurat_withQC.rds")
saveRDS(ind7_seurat, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/individual seurat objects with QC/ind7_seurat_withQC.rds")
#______________________________________________________________________________________________________________________

#Read in RDS Seurat objects

ind4_seurat <- readRDS("C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/individual seurat objects with QC/ind4_seurat_withQC.rds")
ind5_seurat <- readRDS("C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/individual seurat objects with QC/ind5_seurat_withQC.rds")
ind6_seurat <- readRDS("C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/individual seurat objects with QC/ind6_seurat_withQC.rds")
ind7_seurat <- readRDS("C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/individual seurat objects with QC/ind7_seurat_withQC.rds")

#Normalize Seurat objects and then run FindVariableFeatures() on them
ind4_seurat <- NormalizeData(ind4_seurat, verbose = FALSE)
ind5_seurat <- NormalizeData(ind5_seurat, verbose = FALSE)
ind6_seurat <- NormalizeData(ind6_seurat, verbose = FALSE)
ind7_seurat <- NormalizeData(ind7_seurat, verbose = FALSE)

ind4_seurat <- FindVariableFeatures(ind4_seurat, selection.method = "vst")
ind5_seurat <- FindVariableFeatures(ind5_seurat, selection.method = "vst")
ind6_seurat <- FindVariableFeatures(ind6_seurat, selection.method = "vst")
ind7_seurat <- FindVariableFeatures(ind7_seurat, selection.method = "vst")

## Save normalized seurat objects as RDS files
saveRDS(ind4_seurat, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/seurat objects/ind4_seurat_NormalizeData.rds")
saveRDS(ind5_seurat, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/seurat objects/ind5_seurat_NormalizeData.rds")
saveRDS(ind6_seurat, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/seurat objects/ind6_seurat_NormalizeData.rds")
saveRDS(ind7_seurat, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/seurat objects/ind7_seurat_NormalizeData.rds")

#___________________________________________________________________________________________________________________________
ind4_seurat <- readRDS("C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/seurat objects/normal normalization/ind4_seurat_NormalizeData.rds")
ind5_seurat <- readRDS("C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/seurat objects/normal normalization/ind5_seurat_NormalizeData.rds")
ind6_seurat <- readRDS("C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/seurat objects/normal normalization/ind6_seurat_NormalizeData.rds")
ind7_seurat <- readRDS("C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/seurat objects/normal normalization/ind7_seurat_NormalizeData.rds")

# Integrating Seurat objects
integrate.list <- objects()
integrate.list$ind4 <- ind4_seurat
integrate.list$ind5 <- ind5_seurat
integrate.list$ind6 <- ind6_seurat
integrate.list$ind7 <- ind7_seurat

reference.list <- integrate.list[c("ind4", "ind5", "ind6", "ind7")]
integrate.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:40)
integrated.data <- IntegrateData(anchorset = integrate.anchors, dims = 1:40)

saveRDS(integrated.data, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/seurat objects/integrated_data_NormalizeData")
#_________________________________________________________________________________________________________________________________________
integrated.data <- readRDS(file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/seurat objects/integrated_data_NormalizeData")

DefaultAssay(integrated.data) <- "integrated"

integrated.data <- ScaleData(object = integrated.data, verbose = FALSE)
integrated.data <- RunPCA(object = integrated.data, npcs = 40, verbose = FALSE)

integrated.data <- FindNeighbors(object = integrated.data, dims = 1:40)

saveRDS(integrated.data, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/seurat objects/integrated_data_NormalizeData_BeforeFindClusters")

#______________________________________________________________________________________________________________________________________________________
integrated.data <- readRDS(file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/seurat objects/integrated_data_NormalizeData_BeforeFindClusters")

integrated.data <- FindClusters(object = integrated.data, dims.use = 1:40, resolution = 0.01)

integrated.data <- RunUMAP(integrated.data, reduction = "pca", dims = 1:40)

DimPlot(integrated.data, reduction = "umap")

saveRDS(integrated.data, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/seurat objects/NDintegrated_data_res0.1")

rm(integrated.data)

#______________________________________________________________________________________________________________________________________________________________________________
integrated.data <- readRDS(file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/seurat objects/integrated_data_NormalizeData_BeforeFindClusters")

integrated.data <- FindClusters(object = integrated.data, dims.use = 1:40, resolution = 0.5)

integrated.data <- RunUMAP(integrated.data, reduction = "pca", dims = 1:40)

DimPlot(integrated.data, reduction = "umap", label = TRUE)

saveRDS(integrated.data, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/seurat objects/NDintegrated_data_res0.7")

#____________________________________________________________________________________________________________________________________

res_0.1 <- readRDS(file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/seurat objects/NDintegrated_data_res0.1")
 
res_0.1.markers <- FindAllMarkers(res_0.1, only.pos = T)

saveRDS(res_0.1.markers, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/markers/res0.1_markers")

FeaturePlot(res_0.1, "KRT14", min.cutoff = "q9")

FeaturePlot(res_0.1, "KRT18", min.cutoff = "q9")

FeaturePlot(res_0.1, "SLPI", min.cutoff = "q9")

FeaturePlot(res_0.1, "SYTL2", min.cutoff = "q9")

DoHeatmap(res_0.1, features = c("KRT14", "KRT5", "ACTA2", "MYLK", "TP63", "SLPI", "PROM1", "KRT19", "ANKRD30A", "SYTL2"))


#_____________________________________________________________________________________________

res_0.5 <- readRDS(file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/seurat objects/NDintegrated_data_res0.5")

res_0.5.markers <- FindAllMarkers(res_0.5, only.pos = T)

saveRDS(res_0.5.markers, file = "C:/Users/jainar/Documents/Analysis of Normal Breast Epithelium 10x Data/NormalizeData/markers/reS_0.5.markers")


