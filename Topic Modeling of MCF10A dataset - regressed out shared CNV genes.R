# Creating MCF10A object while regressing out genes in CNV regions that overlap with variable genes

library(Seurat)

mcf10A.data <- Read10X(data.dir = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/filtered_feature_bc_matrix/")

mcf10A_seurat <- CreateSeuratObject(counts = mcf10A.data$`Gene Expression`, min.cells = 5, min.features = 200)
barcode_ID <- mcf10A.data[["Antibody Capture"]][1:5,]
ADTassay <- (counts= barcode_ID)
missing_link <- setdiff(colnames(x = CreateAssayObject(counts = barcode_ID)), colnames(x = mcf10A_seurat))
ADTassay <- ADTassay[ , -which(colnames(ADTassay) %in% missing_link)]
mcf10A_seurat[["ADT"]] <- CreateAssayObject(ADTassay)

expression_ID <- mcf10A.data[["Antibody Capture"]][6:9,]
EXPassay <- (counts= expression_ID)
exp_missing_link <- setdiff(colnames(x = CreateAssayObject(counts = expression_ID)), colnames(x = mcf10A_seurat))
EXPassay <- EXPassay[ , -which(colnames(EXPassay) %in% exp_missing_link)]
mcf10A_seurat[["EXP"]] <- CreateAssayObject(EXPassay)

mcf10A_seurat <- NormalizeData(mcf10A_seurat, assay = "EXP", normalization.method = "CLR")
mcf10A_seurat <- ScaleData(mcf10A_seurat, assay = "EXP")

mcf10A_seurat[["percent.mt"]] <- PercentageFeatureSet(object = mcf10A_seurat, pattern = "^MT-")

VlnPlot(object = mcf10A_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

mcf10A_seurat <- subset(x = mcf10A_seurat, subset = nFeature_RNA > 300 & nFeature_RNA < 3300 & percent.mt < 15)

mcf10A_seurat <- NormalizeData(mcf10A_seurat)
mcf10A_seurat <- FindVariableFeatures(object = mcf10A_seurat, selection.method = "vst", nfeatures = 2000)
Mall.genes <- rownames(x = mcf10A_seurat)

mcf10A_seurat <- ScaleData(object = mcf10A_seurat, features = Mall.genes)

######

mcf10a_var_features <- VariableFeatures(mcf10A_seurat)

shared_var_features <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/InferCNV/shared_var_features")

features <- mcf10a_var_features[! mcf10a_var_features %in% shared_var_features]

mcf10A_seurat <- RunPCA(object = mcf10A_seurat, features = features)
#ElbowPlot(mcf10A_seurat, ndims = 50)

mcf10A_seurat <- NormalizeData(mcf10A_seurat, assay = "ADT", normalization.method = "CLR")
mcf10A_seurat <- ScaleData(mcf10A_seurat, assay = "ADT")
mcf10A_seurat <- HTODemux(mcf10A_seurat, assay = "ADT", positive.quantile = 0.99)

table(mcf10A_seurat$ADT_classification.global)

mcf10A_seurat@meta.data$Stage <- "HERE"
for(i in seq(1,length(mcf10A_seurat$ADT_classification))) {
  if(mcf10A_seurat$ADT_classification[i] == "CTL") {
    mcf10A_seurat@meta.data$Stage[i] <- "Control"
  }
  if(mcf10A_seurat$ADT_classification[i] == "NeoT") {
    mcf10A_seurat@meta.data$Stage[i] <- "NeoT"
  }
  if(mcf10A_seurat$ADT_classification[i] == "AT1") {
    mcf10A_seurat@meta.data$Stage[i] <- "AT1"
  }
  if(mcf10A_seurat$ADT_classification[i] == "DCIS") {
    mcf10A_seurat@meta.data$Stage[i] <- "DCIS"
  }
  if(mcf10A_seurat$ADT_classification[i] == "CA1.1") {
    mcf10A_seurat@meta.data$Stage[i] <- "CA1"
  }
}

mcf10A_seurat <- subset(mcf10A_seurat, subset = ADT_classification.global == "Singlet")

my_levels <- c("Control", "NeoT", "AT1", "DCIS", "CA1")

mcf10A_seurat@meta.data$Stage <- factor(x = mcf10A_seurat@meta.data$Stage, levels = my_levels)

Idents(mcf10A_seurat) <- "Stage"

cc.genes <- readLines(con = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/regev_lab_cell_cycle_genes.txt")

s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

mcf10A_seurat <- CellCycleScoring(mcf10A_seurat, s.features = s.genes, g2m.features = g2m.genes)
mcf10A_seurat <- FindNeighbors(object = mcf10A_seurat, dims = 1:30)
mcf10A_seurat <- FindClusters(object = mcf10A_seurat, resolution = 0.3)
mcf10A_seurat <- RunUMAP(object = mcf10A_seurat, dims = 1:30, n.neighbors = 20, min.dist = 0.2)

saveRDS(mcf10A_seurat, file = "~/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/mcf10A_seurat_cnv.rds")

###############################################################################################################################
# Building LDA model on MCF10a data

library(TITAN)
library(Seurat)
library(tidyverse)
library(knitr)

mcf10A_seurat_cnv <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/mcf10A_seurat_cnv.rds")

memory.limit()
memory.limit(size = 24000)

LDA_model_cnv <- runLDA(mcf10A_seurat_cnv, ntopics = 30, normalizationMethod = "CLR")
saveRDS(LDA_model_cnv, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/LDA_model_cnv")

##########################################################################################################################
# Adding topic info to Seurat object 

## Store topic info within Seurat object 
mcf10A_seurat_cnv <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/mcf10A_seurat_cnv.rds")
LDA_model_cnv <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/LDA_model_cnv")

mcf10A_seurat_topicModeling_cnv <- addTopicsToSeuratObject(model = LDA_model_cnv, Object = mcf10A_seurat_cnv)
saveRDS(mcf10A_seurat_topicModeling_cnv, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/mcf10A_seurat_topicModeling_cnv")

############################################################################################################################
# Analyzing LDA model

## GeneScores() outputs a matrix with a table of the scores for each gene for each topic
## Sorting by topic, one can see the disribution of genes and which genes contribute most to each topic. Below is an example of showing only the top topic 1 genes and subsequent scores of said for the first 10 topics.

mcf10A_seurat_topicModeling_cnv <- readRDS("~/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/mcf10A_seurat_topicModeling_cnv")
LDA_model_cnv <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/LDA_model_cnv")

GeneDistribution <- GeneScores(LDA_model_cnv)

## head and sort the matrix by genes most contributing to Topic 4
head(GeneDistribution[order(GeneDistribution[,"Topic_4"], decreasing = T),],n = 10)

## TopTopicGenes() output a matrix with top n scoring genes for each topic  
## These genes can be used to connect each topic to a certain gene network or ontology

TopicGenes <- TopTopicGenes(LDA_model_cnv, ngenes = 200)
head(TopicGenes[,1:10], 10)

## GetTopics() outputs the contribution each cell belongs to a topic, calculated from the gene score
## Outputs a scaled cell-topic matrix

LDA_topics_cnv <- GetTopics(LDA_model_cnv, mcf10A_seurat_topicModeling_cnv)
kable(head(LDA_topics_cnv[,1:10]))

saveRDS(LDA_topics_cnv, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/LDA_topics_cnv")

## Heatmap to visualize which cells expression which topics

library(Seurat)
library(TITAN)

mcf10A_seurat_topicModeling_cnv <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/mcf10A_seurat_topicModeling_cnv")

LDA_topics_cnv <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/LDA_topics_cnv")

LDA_model_cnv <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/LDA_model_cnv")


HeatmapTopic(Object = mcf10A_seurat_topicModeling_cnv,
             topics =  LDA_topics_cnv,
             AnnoVector = mcf10A_seurat_topicModeling_cnv@meta.data$Stage,
             AnnoName = "Cell Lines")

## Visualize the expression of the topics through Seurat's FeaturePlot function

FeaturePlot(mcf10A_seurat_topicModeling_cnv, pt.size = 0.01, features  = "Topic_4", min.cutoff = "q1")

FeaturePlot(mcf10A_seurat_topicModeling_cnv, pt.size = 0.01, features  = "Topic_26", min.cutoff = 'q1')

FeaturePlot(mcf10A_seurat_topicModeling_cnv, pt.size = 0.01, features  = "Topic_15", min.cutoff = 'q1')

CombinePlots(list(topic14, topic15, topic18, topic27, topic28), ncol = 2)
