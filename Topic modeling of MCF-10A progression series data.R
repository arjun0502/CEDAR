# Background Info

## TITAN (Topic Inference of Transciptionally Associated Networks)
  ## R package that runs Topic Modeling on scRNA-seq data
  ## Topic modeling has been used in text mining to link words and documents to topics
  ## Topics can be thought of as themes, such as categories in news article inference (science, politics, etc)
  ## TITAN applies topic modeling to scRNA-seq to find latent transcriptional topics that link genes to topics and cells to topics
  ## TITAN uses a nlp alghorithm called Latent Dirchlet Allocation (LDA) that maximizes these two probability distributions that show:
      # 1) the probability a gene is assigned to a topic 
      # 2) the probability a cell is assigned to a topic

###############################################################################################################################
# Building LDA model on MCF10a data

library(TITAN)
library(Seurat)
library(tidyverse)
library(knitr)

mcf10A_seurat <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/mcf10A_seurat.rds")

memory.limit()

LDA_model <- runLDA(mcf10A_seurat, ntopics = 30, normalizationMethod = "CLR")
saveRDS(LDA_model, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/LDA_model")

##########################################################################################################################
# Adding topic info to Seurat object and Analyzing Model

## Store topic info within Seurat object 
mcf10A_seurat <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/mcf10A_seurat.rds")
LDA_model <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/LDA_model")

mcf10A_seurat_topicModeling <- addTopicsToSeuratObject(model = LDA_model, Object = mcf10A_seurat)
saveRDS(mcf10A_seurat_topicModeling, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/mcf10A_seurat_topicModeling")

## GeneScores() outputs a matrix with a table of the scores for each gene for each topic
## Sorting by topic, one can see the disribution of genes and which genes contribute most to each topic. Below is an example of showing only the top topic 1 genes and subsequent scores of said for the first 10 topics.

GeneDistribution <- GeneScores(LDA_model)

## head and sort the matrix by genes most contributing to Topic 14
head(GeneDistribution[order(GeneDistribution[,"Topic_14"], decreasing = T),],n = 10)

## TopTopicGenes() output a matrix with top n scoring genes for each topic  
## These genes can be used to connect each topic to a certain gene network or ontology

TopicGenes <- TopTopicGenes(LDA_model, ngenes = 200)
head(TopicGenes[,1:10], 10)

## GetTopics() outputs the contribution each cell belongs to a topic, calculated from the gene score
## Outputs a scaled cell-topic matrix

LDA_topics <- GetTopics(LDA_model, mcf10A_seurat_topicModeling)
kable(head(LDA_topics[,1:10]))

saveRDS(LDA_topics, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/LDA_topics")

## Heatmap to visualize which cells expression which topics

library(Seurat)
library(TITAN)

mcf10A_seurat_topicModeling <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/mcf10A_seurat_topicModeling")

LDA_topics <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/LDA_topics")

LDA_model <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Topic Modeling/LDA_model")


HeatmapTopic(Object = mcf10A_seurat_topicModeling,
             topics =  LDA_topics,
             AnnoVector = mcf10A_seurat_topicModeling@meta.data$Stage,
             AnnoName = "Cell Lines")

## Visualize the expression of the topics through Seurat's FeaturePlot function

FeaturePlot(mcf10A_seurat_topicModeling, pt.size = 0.01, features  = "Topic_14", min.cutoff = "q1")

FeaturePlot(mcf10A_seurat_topicModeling, pt.size = 0.01, features  = "Topic_15", min.cutoff = 'q1')

FeaturePlot(mcf10A_seurat_topicModeling, pt.size = 0.01, features  = "Topic_18", min.cutoff = 'q1')

FeaturePlot(mcf10A_seurat_topicModeling,  pt.size = 0.01, features  = "Topic_27", min.cutoff = 'q1')

FeaturePlot(mcf10A_seurat_topicModeling,  pt.size = 0.01, features  = "Topic_28", min.cutoff = 'q1')

CombinePlots(list(topic14, topic15, topic18, topic27, topic28), ncol = 2)

Oncostatin_genes <- read.table("~/Gap Year/CEDAR Student Worker Position/Analysis of MCF-10A Progression Series/Oncostatin M Signaling Pathway Genes.txt")
Oncostatin_genes <- Oncostatin_genes$V1

DoHeatmap(mcf10A_seurat, features = Oncostatin_genes)

kable(head(TopicGenes[,27], 100))

cat( paste( head(TopicGenes[,27], 200), collapse='\n' ) )

