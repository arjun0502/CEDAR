# Scoring basal-ness and luminal-ness in MCF7 cells
# 4 categories: normal basal (NB), normal luminal (NL), cancer basal (CB), cancer luminal (CL)
# Arjun Jain

# Load Seurat library
library(Seurat)

# Read in MCF7 and healthy objects
MCF7 <- readRDS("C:/Users/jainar/Documents/cc_reduced_time.rds")
Idents(object = MCF7) <- 'RNA_snn_res.0.3'
healthy <- readRDS("C:/Users/jainar/Documents/NDintegrated_data_res0.5")

# Get markers for four categories: normal basal (NB), normal luminal (NL), cancer basal (CB), cancer luminal (CL)

basal_markers <- FindMarkers(healthy, ident.1 = c(1, 4, 5, 6, 9, 11), ident.2 = c(0, 2, 3, 7, 8, 10), only.pos = T)
luminal_markers <- FindMarkers(healthy, ident.1 = c(0, 2, 3, 7, 8, 10), ident.2 = c(1, 4, 5, 6, 9, 11), only.pos = T)
basal_cancer_markers <- list(c("UBE2C", "PTTG1", "MYBL2", "BIRC5", "CCNB1", "TYMS", "MELK", "CEP55", "KNTC2", "UBE2T", "ANLN", "ORC6L", "KIF2C", "EXO1", "CDCA1", "CENPF", "CCNE1", "MK167", "MYC", "MIA", "FOXC1", "ACTR3B", "PHGCH", "CDH3"))
luminal_cancer_markers <- list(c("TMEM45B", "BAG1", "PGR", "MAPT", "GPR160", "CXXC5", "ESR1", "NAT1", "FOXA1", "BLVRA"))

# Add meta data columns for each category containing expression of markers for respective category
top50_basal <- list(rownames(basal_markers)[1:50])
top50_luminal <- list(rownames(luminal_markers)[1:50])
MCF7 <- AddModuleScore(MCF7, features = top50_basal, name = "NB_scores")
MCF7 <- AddModuleScore(MCF7, features = top50_luminal, name = "NL_scores")
MCF7 <- AddModuleScore(MCF7, features = basal_cancer_markers, name = "CB_scores")
MCF7 <- AddModuleScore(MCF7, features = luminal_cancer_markers, name = "CL_scores")

#Scale module scores to range of scale data in MCF7 data
MCF7@meta.data$NB_scores_scaled <- scale(MCF7@meta.data$NB_scores1, center = T, scale=T)
MCF7@meta.data$NL_scores_scaled <- scale(MCF7@meta.data$NL_scores1, center = T, scale=T)
MCF7@meta.data$CB_scores_scaled <- scale(MCF7@meta.data$CB_scores1, center = T, scale=T)
MCF7@meta.data$CL_scores_scaled <- scale(MCF7@meta.data$CL_scores1, center = T, scale=T)

saveRDS(MCF7, "C:/Users/jainar/Documents/MCF7_basal_luminal_scoring")
#___________________________________________________________________________
# CORRELATION between module scores and gene expression 

correlation_module_scores_gene_exp(MCF7, "NB_scores_scaled", "spearman", top50_basal)
correlation_module_scores_gene_exp(MCF7, "NL_scores_scaled", "spearman", top50_luminal)
correlation_module_scores_gene_exp(MCF7, "CB_scores_scaled", "spearman", basal_cancer_markers)
correlation_module_scores_gene_exp(MCF7, "CL_scores_scaled", "spearman", luminal_cancer_markers)

correlation_module_scores_gene_exp(MCF7, "NB_scores_scaled", "pearson", top50_basal)
correlation_module_scores_gene_exp(MCF7, "NL_scores_scaled", "pearson", top50_luminal)
correlation_module_scores_gene_exp(MCF7, "CB_scores_scaled", "pearson", basal_cancer_markers)
correlation_module_scores_gene_exp(MCF7, "CL_scores_scaled", "pearson", luminal_cancer_markers)
