# Looking at each cell line of MCF10A data separately and finding DE genes across subclusters of each cell line

library(Seurat)
library(CytoTRACE)
library(tidyverse)
library(ggpubr)

mcf10A_seurat <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/mcf10A_seurat.rds")
DimPlot(mcf10A_seurat)

#________________________________________________________________________________________________________________________________________________________
# Subclusters in each cell line 

CTL_seurat <- subset(mcf10A_seurat, subset = Stage == "Control")
NeoT_seurat <- subset(mcf10A_seurat, subset = Stage == "NeoT")
AT1_seurat <- subset(mcf10A_seurat, subset = Stage == "AT1")
DCIS_seurat <- subset(mcf10A_seurat, subset = Stage == "DCIS")
CA1_seurat <- subset(mcf10A_seurat, subset = Stage == "CA1")

CTL_seurat <- NormalizeData(CTL_seurat)
NeoT_seurat <- NormalizeData(NeoT_seurat)
AT1_seurat <- NormalizeData(AT1_seurat)
DCIS_seurat <- NormalizeData(DCIS_seurat)
CA1_seurat <- NormalizeData(CA1_seurat)

CTL_seurat <- FindVariableFeatures(object = CTL_seurat, selection.method = "vst", nfeatures = 2000)
NeoT_seurat <- FindVariableFeatures(object = NeoT_seurat, selection.method = "vst", nfeatures = 2000)
AT1_seurat <- FindVariableFeatures(object = AT1_seurat, selection.method = "vst", nfeatures = 2000)
DCIS_seurat <- FindVariableFeatures(object = DCIS_seurat, selection.method = "vst", nfeatures = 2000)
CA1_seurat <- FindVariableFeatures(object = CA1_seurat, selection.method = "vst", nfeatures = 2000)

CTLall.genes <- rownames(x = CTL_seurat)
NeoTall.genes <- rownames(x = NeoT_seurat)
AT1all.genes <- rownames(x = AT1_seurat)
DCISall.genes <- rownames(x = DCIS_seurat)
CA1all.genes <- rownames(x = CA1_seurat)

CTL_seurat <- ScaleData(object = CTL_seurat, features = CTLall.genes)
NeoT_seurat <- ScaleData(object = NeoT_seurat, features = NeoTall.genes)
AT1_seurat <- ScaleData(object = AT1_seurat, features = AT1all.genes)
DCIS_seurat <- ScaleData(object = DCIS_seurat, features = DCISall.genes)
CA1_seurat <- ScaleData(object = CA1_seurat, features = CA1all.genes)

CTL_seurat <- RunPCA(object = CTL_seurat, features = VariableFeatures(object = CTL_seurat))
NeoT_seurat <- RunPCA(object = NeoT_seurat, features = VariableFeatures(object = NeoT_seurat))
AT1_seurat <- RunPCA(object = AT1_seurat, features = VariableFeatures(object = AT1_seurat))
DCIS_seurat <- RunPCA(object = DCIS_seurat, features = VariableFeatures(object = DCIS_seurat))
CA1_seurat <- RunPCA(object = CA1_seurat, features = VariableFeatures(object = CA1_seurat))

CTL_seurat <- FindNeighbors(object = CTL_seurat, dims = 1:30)
NeoT_seurat <- FindNeighbors(object = NeoT_seurat, dims = 1:30)
AT1_seurat <- FindNeighbors(object = AT1_seurat, dims = 1:30)
DCIS_seurat <- FindNeighbors(object = DCIS_seurat, dims = 1:30)
CA1_seurat <- FindNeighbors(object = CA1_seurat, dims = 1:30)

cc.genes <- readLines(con = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

CTL_seurat <- CellCycleScoring(CTL_seurat, s.features = s.genes, g2m.features = g2m.genes)
NeoT_seurat <- CellCycleScoring(NeoT_seurat, s.features = s.genes, g2m.features = g2m.genes)
AT1_seurat <- CellCycleScoring(AT1_seurat, s.features = s.genes, g2m.features = g2m.genes)
DCIS_seurat <- CellCycleScoring(DCIS_seurat, s.features = s.genes, g2m.features = g2m.genes)
CA1_seurat <- CellCycleScoring(CA1_seurat, s.features = s.genes, g2m.features = g2m.genes)

CTL_seurat <- ScaleData(CTL_seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(CTL_seurat))
NeoT_seurat <- ScaleData(NeoT_seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(NeoT_seurat))
AT1_seurat <- ScaleData(AT1_seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(AT1_seurat))
DCIS_seurat <- ScaleData(DCIS_seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(DCIS_seurat))
CA1_seurat <- ScaleData(CA1_seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(CA1_seurat))

CTL_seurat <- RunPCA(object = CTL_seurat, features = VariableFeatures(object = CTL_seurat))
NeoT_seurat <- RunPCA(object = NeoT_seurat, features = VariableFeatures(object = NeoT_seurat))
AT1_seurat <- RunPCA(object = AT1_seurat, features = VariableFeatures(object = AT1_seurat))
DCIS_seurat <- RunPCA(object = DCIS_seurat, features = VariableFeatures(object = DCIS_seurat))
CA1_seurat <- RunPCA(object = CA1_seurat, features = VariableFeatures(object = CA1_seurat))

CTL_seurat <- FindNeighbors(object = CTL_seurat, dims = 1:30)
NeoT_seurat <- FindNeighbors(object = NeoT_seurat, dims = 1:30)
AT1_seurat <- FindNeighbors(object = AT1_seurat, dims = 1:30)
DCIS_seurat <- FindNeighbors(object = DCIS_seurat, dims = 1:30)
CA1_seurat <- FindNeighbors(object = CA1_seurat, dims = 1:30)

CTL_seurat <- FindClusters(object = CTL_seurat, resolution = 0.5)
NeoT_seurat <- FindClusters(object = NeoT_seurat, resolution = 0.5)
AT1_seurat <- FindClusters(object = AT1_seurat, resolution = 0.5)
DCIS_seurat <- FindClusters(object = DCIS_seurat, resolution = 0.5)
CA1_seurat <- FindClusters(object = CA1_seurat, resolution = 0.5)

CTL_seurat <- RunUMAP(object = CTL_seurat, dims = 1:30, n.neighbors = 20, min.dist = 0.1)
NeoT_seurat <- RunUMAP(object = NeoT_seurat, dims = 1:30, n.neighbors = 20, min.dist = 0.1)
AT1_seurat <- RunUMAP(object = AT1_seurat, dims = 1:30, n.neighbors = 20, min.dist = 0.1)
DCIS_seurat <- RunUMAP(object = DCIS_seurat, dims = 1:30, n.neighbors = 20, min.dist = 0.1)
CA1_seurat <- RunUMAP(object = CA1_seurat, dims = 1:30, n.neighbors = 20, min.dist = 0.1)

saveRDS(CTL_seurat, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/CTL_seurat")
saveRDS(NeoT_seurat, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/NeoT_seurat")
saveRDS(AT1_seurat, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/AT1_seurat")
saveRDS(DCIS_seurat, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/DCIS_seurat")
saveRDS(CA1_seurat, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/CA1_seurat")

CTL_seurat <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/CTL_seurat")
NeoT_seurat <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/NeoT_seurat")
AT1_seurat <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/AT1_seurat")
DCIS_seurat <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/DCIS_seurat")
CA1_seurat <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/CA1_seurat")

DimPlot(object = CTL_seurat, reduction = "umap", group.by = "seurat_clusters")


# Finding DE genes across subclusters of each cell line 

Idents(CTL_seurat) <- "seurat_clusters"
CTL_seurat.markers <- FindAllMarkers(object = CTL_seurat, logfc.threshold = 0.2)
CTL_seurat_cluster0.markers <- CTL_seurat.markers[which(CTL_seurat.markers$cluster == 0), ]
CTL_seurat_cluster1.markers <- CTL_seurat.markers[which(CTL_seurat.markers$cluster == 1), ]
CTL_seurat_cluster2.markers <- CTL_seurat.markers[which(CTL_seurat.markers$cluster == 2), ]

CTL_seurat_cluster0.markers <- CTL_seurat_cluster0.markers[order(CTL_seurat_cluster0.markers$avg_logFC), ]
CTL_seurat_cluster1.markers <- CTL_seurat_cluster1.markers[order(CTL_seurat_cluster1.markers$avg_logFC), ]
CTL_seurat_cluster2.markers <- CTL_seurat_cluster2.markers[order(CTL_seurat_cluster2.markers$avg_logFC), ]

DE_genes_shared_across_subclusters_CTL <- intersect(CTL_seurat_cluster0.markers$gene, CTL_seurat_cluster1.markers$gene)

CTL_seurat.markers <- rbind(CTL_seurat_cluster0.markers, CTL_seurat_cluster1.markers, CTL_seurat_cluster2.markers)

Idents(NeoT_seurat) <- "seurat_clusters"
NeoT_seurat.markers <- FindAllMarkers(object = NeoT_seurat, logfc.threshold = 0.2)
NeoT_seurat_cluster0.markers <- NeoT_seurat.markers[which(NeoT_seurat.markers$cluster == 0), ]
NeoT_seurat_cluster1.markers <- NeoT_seurat.markers[which(NeoT_seurat.markers$cluster == 1), ]

NeoT_seurat_cluster0.markers <- NeoT_seurat_cluster0.markers[order(NeoT_seurat_cluster0.markers$avg_logFC), ]
NeoT_seurat_cluster1.markers <- NeoT_seurat_cluster1.markers[order(NeoT_seurat_cluster1.markers$avg_logFC), ]

NeoT_seurat.markers <- rbind(NeoT_seurat_cluster0.markers, NeoT_seurat_cluster1.markers)

Idents(AT1_seurat) <- "seurat_clusters"
AT1_seurat.markers <- FindAllMarkers(object = AT1_seurat, logfc.threshold = 0.2)
AT1_seurat_cluster0.markers <- AT1_seurat.markers[which(AT1_seurat.markers$cluster == 0), ]
AT1_seurat_cluster1.markers <- AT1_seurat.markers[which(AT1_seurat.markers$cluster == 1), ]
AT1_seurat_cluster2.markers <- AT1_seurat.markers[which(AT1_seurat.markers$cluster == 2), ]

AT1_seurat_cluster0.markers <- AT1_seurat_cluster0.markers[order(AT1_seurat_cluster0.markers$avg_logFC), ]
AT1_seurat_cluster1.markers <- AT1_seurat_cluster1.markers[order(AT1_seurat_cluster1.markers$avg_logFC), ]
AT1_seurat_cluster2.markers <- AT1_seurat_cluster2.markers[order(AT1_seurat_cluster2.markers$avg_logFC), ]

AT1_seurat.markers <- rbind(AT1_seurat_cluster0.markers, AT1_seurat_cluster1.markers, AT1_seurat_cluster2.markers)

Idents(DCIS_seurat) <- "seurat_clusters"
DCIS_seurat.markers <- FindAllMarkers(object = DCIS_seurat, logfc.threshold = 0.2)
DCIS_seurat_cluster0.markers <- DCIS_seurat.markers[which(DCIS_seurat.markers$cluster == 0), ]
DCIS_seurat_cluster1.markers <- DCIS_seurat.markers[which(DCIS_seurat.markers$cluster == 1), ]

DCIS_seurat_cluster0.markers <- DCIS_seurat_cluster0.markers[order(DCIS_seurat_cluster0.markers$avg_logFC), ]
DCIS_seurat_cluster1.markers <- DCIS_seurat_cluster1.markers[order(DCIS_seurat_cluster1.markers$avg_logFC), ]

DCIS_seurat.markers <- rbind(DCIS_seurat_cluster0.markers, DCIS_seurat_cluster1.markers)

Idents(CA1_seurat) <- "seurat_clusters"
CA1_seurat.markers <- FindAllMarkers(object = CA1_seurat, logfc.threshold = 0.2)
CA1_seurat_cluster0.markers <- CA1_seurat.markers[which(CA1_seurat.markers$cluster == 0), ]
CA1_seurat_cluster1.markers <- CA1_seurat.markers[which(CA1_seurat.markers$cluster == 1), ]
CA1_seurat_cluster2.markers <- CA1_seurat.markers[which(CA1_seurat.markers$cluster == 2), ]

CA1_seurat_cluster0.markers <- CA1_seurat_cluster0.markers[order(CA1_seurat_cluster0.markers$avg_logFC), ]
CA1_seurat_cluster1.markers <- CA1_seurat_cluster1.markers[order(CA1_seurat_cluster1.markers$avg_logFC), ]
CA1_seurat_cluster2.markers <- CA1_seurat_cluster2.markers[order(CA1_seurat_cluster2.markers$avg_logFC), ]

CA1_seurat.markers <- rbind(CA1_seurat_cluster0.markers, CA1_seurat_cluster1.markers, CA1_seurat_cluster2.markers)

library(xlsx)

saveRDS(CTL_seurat.markers, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/CTL_seurat.markers")
saveRDS(NeoT_seurat.markers, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/NeoT_seurat.markers")
saveRDS(AT1_seurat.markers, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/AT1_seurat.markers")
saveRDS(DCIS_seurat.markers, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/DCIS_seurat.markers")
saveRDS(CA1_seurat.markers, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/CA1_seurat.markers")

write.xlsx(CTL_seurat.markers, "c:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/CTL_seurat_markers.xlsx")
write.xlsx(NeoT_seurat.markers, "c:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/NeoT_seurat_markers.xlsx")
write.xlsx(AT1_seurat.markers, "c:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/AT1_seurat_markers.xlsx")
write.xlsx(DCIS_seurat.markers, "c:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/DCIS_seurat_markers.xlsx")
write.xlsx(CA1_seurat.markers, "c:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/CA1_seurat_markers.xlsx")

CTL_seurat.genes <-CTL_seurat.markers$gene
NeoT_seurat.genes <- NeoT_seurat.markers$gene
AT1_seurat.genes <- AT1_seurat.markers$gene
DCIS_seurat.genes <- DCIS_seurat.markers$gene
CA1_seurat.genes <- CA1_seurat.markers$gene

# SERPINB1, MT2A, KRT16, SAA1, S100P, KRT6A, MT1E

# 5 DE genes across subclusters shared across all cell lines (so small because DCIS has only 28 DE genes)
# "KRT6A"   "NCL"     "MT-ND4L" "MT-ND6"  "FST" 

DE_genes_shared_across_cell_lines2 <- intersect(intersect(intersect(intersect(CTL_seurat.genes, NeoT_seurat.genes), AT1_seurat.genes), CA1_seurat.genes), DCIS_seurat.genes)

DE_genes_shared_across_cell_lines <- intersect(intersect(intersect(CTL_seurat.genes, NeoT_seurat.genes), AT1_seurat.genes), CA1_seurat.genes)

# 20 DE genes across subclusters shared across all cell lines excluding DCIS since it only has 28 DE genes
# "KRT6A"    "KRT16"    "SLPI"     "HSP90AA1" "NCL"      "MT-ND4L"  "MT-ND6"  
# "S100A10"  "MT2A"     "MT1E"     "CLDN7"    "SERPINB1" "LGALS3"   "KLK10"   
# "S100P"    "FST"      "TACSTD2"  "LAMC2"    "FXYD3"    "SAA1"   

# Ran genes through EnrichR and found that 7 of the genes are enriched in the Oncostatin M Signaling Pathway
# SERPINB1, MT2A, KRT16, SAA1, S100P, KRT6A, MT1E

Oncostatin_M_pathway_genes <- c("SERPINB1", "MT2A", "KRT16", "SAA1", "S100P", "KRT6A", "MT1E")

