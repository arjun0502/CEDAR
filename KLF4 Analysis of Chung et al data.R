library(Seurat)
library(TITAN)
#library(devtools)
#install_github("JuliusCampbell/TITAN")

scBC_SO <- readRDS(file = "C:/Users/jainar/Documents/Arjun High School/High School Science Research/OHSU Internship - CEDAR/KLF4 Analysis of Chung et al dataset/scBC_SO.rds")
DimPlot(scBC_SO, reduction = "umap", label = TRUE)
Idents(scBC_SO)

MCF7_KLF4_var_genes <- read.table(file = "C:/Users/jainar/Documents/Arjun High School/High School Science Research/OHSU Internship - CEDAR/KLF4 Analysis of Chung et al dataset/Model_MCF7_KLF4_5000Variable_top50_genes_topics.txt", header = T)
T47D_KLF4_var_genes <- read.table(file = "C:/Users/jainar/Documents/Arjun High School/High School Science Research/OHSU Internship - CEDAR/KLF4 Analysis of Chung et al dataset/Model_T47D_KLF4_20T_CLR_5000Variable_M10_top50_genes_topics.txt", header = T)

DefaultAssay(scBC_SO) <- "RNA"
Model_T47D_KLF4 <- readRDS(file = "C:/Users/jainar/Documents/Arjun High School/High School Science Research/OHSU Internship - CEDAR/KLF4 Analysis of Chung et al dataset/Model_T47D_KLF4_20T_CLR_5000Variable_M10.rds")
Model_MCF7_KLF4 <- readRDS(file = "C:/Users/jainar/Documents/Arjun High School/High School Science Research/OHSU Internship - CEDAR/KLF4 Analysis of Chung et al dataset/Model_MCF7_KLF4_20T_CLR_5000Variable_M10.rds")


scBC_SO <- ImputeAndAddTopics(scBC_SO, Model_T47D_KLF4, TopicPrefix = "T47D_KLF4_Topics")
HeatmapTopic(Object = scBC_SO, topics = Embeddings(scBC_SO, "imputedLDA"), AnnoVector = scBC_SO@meta.data$orig.ident, AnnoName = "Timepoint")

