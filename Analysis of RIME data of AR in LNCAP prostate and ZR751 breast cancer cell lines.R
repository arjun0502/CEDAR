library(proDA)

# proDA package 
## R package that implements a powerful probabilistic dropout model to identify differentially abundant proteins
## Specifically designed for label-free mass spectrometry data and in particular how to handle the many many missing values (good for RIME data)

#========================================================================================================================

# 1. Load the proteinGroups.txt MaxQuant output table

## MaxQuant - most popular tools for handling raw MS data which produces multiple files
## Important file that contains the protein intensities is called proteinGroups.txt (same as data Hisham gave to me)
## proteinGroups.txt - table with detailed information about the identification and quantification process for each protein group 

maxquant_protein_table <- read.delim(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/proteinGroups-16.txt" ,stringsAsFactors = FALSE)

# 2. Extract the intensity columns and create the abundance matrix

## Table contains a lot of information (359 columns!!), but only interested in the columns which contain the measured intensities

## Use a regular expression (regex) to select the intensity columns (specifically LFQ)

intensity_colnames <- grep("^Intensity", colnames(maxquant_protein_table), value=TRUE)

intensity_colnames <- intensity_colnames[2:22]

intensity_colnames_LFQ <- grep("^LFQ",colnames(maxquant_protein_table), value = TRUE)

head(intensity_colnames)

## Create the intensity matrix
abundance_matrix <- as.matrix(maxquant_protein_table[, intensity_colnames])
abundance_matrix_LFQ <- as.matrix(maxquant_protein_table[, intensity_colnames_LFQ])

## Adapt column and row maxquant_protein_table
rownames(abundance_matrix) <- maxquant_protein_table$Gene.names
rownames(abundance_matrix_LFQ) <- maxquant_protein_table$Gene.names

colnames(abundance_matrix) <- sub("^Intensity.", "", intensity_colnames)
colnames(abundance_matrix_LFQ) <- sub("^LFQ.intensity.", "", intensity_colnames_LFQ)

# 3. Replace the zeros with NAs and take the log2() of the data

## MaxQuant codes missing values as 0
## Misleading because the actual abundance probably was not 0, but just some value too small to be detected by MS
## So, replace all 0 with NA
## Secondly,raw intensity values are such that variance and mean are dependent
## This is undesirable, because change of x units can be large shift if mean is small or irrelevant if mean is large
## To make mean and variance independent, log2 transform the intensities
## Now a change of x units is as significant for highly abundant proteins, as it is for low abundant ones.

abundance_matrix[abundance_matrix == 0] <- NA
abundance_matrix <- log2(abundance_matrix)

abundance_matrix_LFQ[abundance_matrix_LFQ == 0] <- NA
abundance_matrix_LFQ <- log2(abundance_matrix_LFQ)


# 4. Normalize the data using median_normalization()

## Normalize the data to remove potential sample specific effects
## But this step is challenging, because the missing values cannot easily be corrected for
## Thus, a first helpful plot is to look how many missing values are in each sample
# 
# barplot(colSums(is.na(abundance_matrix)),
#          ylab = "# missing values", xlab = "Samples")
# 
# barplot(colSums(is.na(abundance_matrix_LFQ)),
#         ylab = "# missing values", xlab = "Samples")

## Number of missing values differs substantially between samples (between 30% and 90%) in this dataset
## If we take a look at the intensity distribution for each sample, we see that they differ substantially as well
# 
# boxplot(abundance_matrix,
#         ylab = "Intensity Distribution",
#         xlab = "Samples")


# boxplot(abundance_matrix_LFQ,
#         ylab = "Intensity Distribution",
#         xlab = "Samples")

## Intensity distribution shifted upwards for samples which also have a large number of missing values (for example the last one)
## This agrees with our idea that small values are more likely to be missing
## Also demonstrates why normalization methods like quantile normalization, which distort the data until all the distributions are equal, are problematic

## Apply "conservative" median normalization
## Ignores missing values and transforms values so that median difference between sample and average across all other samples is zero

normalized_abundance_matrix <- median_normalization(abundance_matrix)
saveRDS(normalized_abundance_matrix, file = "~/Gap Year/CEDAR Student Worker Position/RIME Project/normalized_abundance_matrix")
normalized_abundance_matrix_LFQ <- median_normalization(abundance_matrix_LFQ)
saveRDS(normalized_abundance_matrix_LFQ, file = "~/Gap Year/CEDAR Student Worker Position/RIME Project/normalized_abundance_matrix_LFQ")

# 5. Inspect sample structure with a heatmap of the distance matrix (dist_approx())

## Important tool to identify sample swaps and outliers in dataset is to look at the sample distance matrix
## It shows the distances of samples A to B, A to C, B to C and so on

## Use dist_approx() function 
## Takes either a fitted model (ie. the output from proDA()) or a simple matrix (for which it internally calls proDA())....
## ......and estimates the expected distance (shown in scale - blue to red) without imputing the missing values
## Also reports the associated uncertainty with every estimate (in every tile, 95% confidence interval)
## Estimates for samples with many missing values will be uncertain (wuide confidence interval), allowing the data analyst to discount them.


da <- dist_approx(normalized_abundance_matrix)
da_LFQ <- dist_approx(normalized_abundance_matrix_LFQ)
 
plot_mat <- as.matrix(da$mean)

# Remove diagonal elements, so that the colorscale is not distorted

plot_mat[diag(21) == 1] <- NA

## 95% conf interval is approx `sd * 1.96`

uncertainty <- matrix(paste0(" ± ",round(as.matrix(da$sd * 1.96), 1)), nrow=21)

library(pheatmap)
pheatmap::pheatmap(plot_mat, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers= uncertainty, number_color = "black")


# 6. Create data frame containing additional info on each sample, in particular to which condition that sample belongs

# The best way to create this data.frame depends on the column naming scheme
sample_info_df <- data.frame(name = colnames(normalized_abundance_matrix), stringsAsFactors = FALSE)
sample_info_df$condition <- substr(sample_info_df$name, 1, nchar(sample_info_df$name)  - 2)
sample_info_df$replicate <- as.numeric(substr(sample_info_df$name, nchar(sample_info_df$name), 20))
sample_info_df

sample_info_df_LNCAP <- sample_info_df[c(1:3, 10:12, 16:18), ]
sample_info_df_ZR751 <- sample_info_df[c(4:9, 13:15, 19:21), ]

#!!!!
## For example, condition CG1407 would be equivalent of AR_DHT_LNCAP and CG1407.01 would be equivalent of AR_DHT_LNCAP_1 
## Each intensity column in proteinGroups.xlsx excel sheet considered a sample (3 samples or replicates per condition)



# 7. Fit the linear probabilistic dropout model to the normalized data

## Call proDA() function to fit the model
## Specify design using the formula notation, referencing the condition column in the sample_info_df data.frame just made
## Negative control for us would be IgG pull down 
# This way automatically all coefficients measure how much each condition differs from the negative control.

normalized_abundance_matrix_LNCAP <- normalized_abundance_matrix[, c(1:3, 10:12, 16:18)]
normalized_abundance_matrix_LFQ_LNCAP <- normalized_abundance_matrix_LFQ[, c(1:3, 10:12, 16:18)]
normalized_abundance_matrix_ZR751 <- normalized_abundance_matrix[, c(4:9, 13:15, 19:21)]
normalized_abundance_matrix_LFQ_ZR751 <- normalized_abundance_matrix_LFQ[, c(4:9, 13:15, 19:21)]


fit_LNCAP <- proDA(normalized_abundance_matrix_LNCAP, design = ~ condition, col_data = sample_info_df_LNCAP, reference_level = "IgG_LNCAP")
saveRDS(fit_LNCAP, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_LNCAP")

fit_LNCAP_LFQ <- proDA(normalized_abundance_matrix_LFQ_LNCAP, design = ~ condition, col_data = sample_info_df_LNCAP, reference_level = "IgG_LNCAP")
saveRDS(fit_LNCAP_LFQ, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_LNCAP_LFQ")

fit_ZR751 <- proDA(normalized_abundance_matrix_ZR751, design = ~ condition, col_data = sample_info_df_ZR751, reference_level = "IgG_ZR751")
saveRDS(fit_ZR751, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_ZR751")

fit_ZR751_LFQ <- proDA(normalized_abundance_matrix_LFQ_ZR751, design = ~ condition, col_data = sample_info_df_ZR751, reference_level = "IgG_ZR751")
saveRDS(fit_ZR751_LFQ, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_ZR751_LFQ")

fit_noref <- proDA(normalized_abundance_matrix, design = ~ condition + 0, col_data = sample_info_df, reference_level = NULL)
saveRDS(fit_noref, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_noref")

fit_noref_LFQ <- proDA(normalized_abundance_matrix_LFQ, design = ~ condition + 0, col_data = sample_info_df, reference_level = NULL)
saveRDS(fit_noref_LFQ, file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_noref_LFQ")


# 8. Identify Differential Abundance 

## Use a Wald test to identify in which proteins a coefficient is significantly different from zero
## The test_diff() function takes first the fit object produced by proDA() and a contrast argument
## This can either be a string or an expression if we want to test more complex combinations
## For example conditionCG1407 - (conditionCG6017 + conditionCG5880) / 2 would test for....
##.....difference between CG1407 and the average of CG6017 and CG5880.
## Alternatively test_diff() also supports likelihood ratio F-tests
## In that case instead of the contrast argument specify the reduced_model argument.

## Test which proteins differ between conditions
## IgG is the default contrast, because it was specified as the `reference_level`

library(proDA)

fit_LNCAP <- readRDS("C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_LNCAP")
fit_LNCAP_LFQ <- readRDS("C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_LNCAP_LFQ")
fit_ZR751 <- readRDS("C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_ZR751")
fit_ZR751_LFQ <- readRDS("C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_ZR751_LFQ")

fit_noref <- readRDS("C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_noref")
fit_noref_LFQ <- readRDS("C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_noref_LFQ")


# test_res_ARvsIgG_LNCAP <- test_diff(fit_LNCAP, contrast = "conditionAR_VEH_LNCAP")
# write.csv(test_res_ARvsIgG_LNCAP, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_ARvsIgG_LNCAP.csv", row.names = F)
# 
# test_res_ARvsIgG_ZR751 <- test_diff(fit_ZR751, contrast = "conditionAR_VEH_ZR751")
# write.csv(test_res_ARvsIgG_ZR751, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_ARvsIgG_ZR751.csv", row.names = F)
# 
# test_res_ARvsIgG_LNCAP_LFQ <- test_diff(fit_LNCAP_LFQ, contrast = "conditionAR_VEH_LNCAP")
# write.csv(test_res_ARvsIgG_LNCAP_LFQ, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_ARvsIgG_LNCAP_LFQ.csv", row.names = F)
# 
# test_res_ARvsIgG_ZR751_LFQ <- test_diff(fit_ZR751_LFQ, contrast = "conditionAR_VEH_ZR751")
#write.csv(test_res_ARvsIgG_ZR751_LFQ, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_ARvsIgG_ZR751_LFQ.csv", row.names = F)

result_names(fit_LNCAP)

test_res_VEHvsDHT_LNCAP <- test_diff(fit_LNCAP, contrast = conditionAR_DHT_LNCAP - conditionAR_VEH_LNCAP)
write.csv(test_res_VEHvsDHT_LNCAP, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_LNCAP.csv", row.names = F)

test_res_VEHvsDHT_LNCAP_LFQ <- test_diff(fit_LNCAP_LFQ, contrast = conditionAR_DHT_LNCAP - conditionAR_VEH_LNCAP)
write.csv(test_res_VEHvsDHT_LNCAP_LFQ, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_LNCAP_LFQ.csv", row.names = F)


result_names(fit_ZR751)

test_res_VEHvsDHT_ZR751 <- test_diff(fit_ZR751, contrast = conditionAR_DHT_ZR751 - conditionAR_VEH_ZR751)
write.csv(test_res_VEHvsDHT_ZR751, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_ZR751.csv", row.names = F)

test_res_VEHvsDHT_ZR751_LFQ <- test_diff(fit_ZR751_LFQ, contrast = conditionAR_DHT_ZR751 - conditionAR_VEH_ZR751)
write.csv(test_res_VEHvsDHT_ZR751_LFQ, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_ZR751_LFQ.csv", row.names = F)


result_names(fit_noref)

test_res_LNCAP_VEH_vs._ZR751_VEH <- test_diff(fit_noref, contrast = conditionAR_VEH_LNCAP - conditionAR_VEH_ZR751)
write.csv(test_res_LNCAP_VEH_vs._ZR751_VEH, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_VEH_vs._ZR751_VEH.csv", row.names = F)

test_res_LNCAP_DHT_vs._ZR751_DHT <- test_diff(fit_noref, contrast = conditionAR_DHT_LNCAP - conditionAR_DHT_ZR751)
write.csv(test_res_LNCAP_DHT_vs._ZR751_DHT, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_DHT_vs._ZR751_DHT.csv", row.names = F)

test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ <- test_diff(fit_noref_LFQ, contrast = conditionAR_VEH_LNCAP - conditionAR_VEH_ZR751)
write.csv(test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ.csv", row.names = F)

test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ <- test_diff(fit_noref_LFQ, contrast = conditionAR_DHT_LNCAP - conditionAR_DHT_ZR751)
write.csv(test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ.csv", row.names = F)



#############################################################################################################################
##Generating Volcano Plots
library(xlsx)

test_res_VEHvsDHT_LNCAP <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_LNCAP.csv", header = T)
test_res_VEHvsDHT_LNCAP_LFQ <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_LNCAP_LFQ.csv", header = T)
test_res_VEHvsDHT_ZR751 <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_ZR751.csv", header = T)
test_res_VEHvsDHT_ZR751_LFQ <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_ZR751_LFQ.csv", header = T)

test_res_LNCAP_VEH_vs._ZR751_VEH <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_VEH_vs._ZR751_VEH.csv")
test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ.csv")
test_res_LNCAP_DHT_vs._ZR751_DHT <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_DHT_vs._ZR751_DHT.csv")
test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ.csv")

TF_Chromatin_genes <- read.xlsx(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/List of TF and Chromatin modifiers.xlsx", sheetName = "Combined list", header = F)
TF_Chromatin_genes <- TF_Chromatin_genes$X1


matching_genes <- function(df, key_genes) {
  
  df$is_key_gene <- NULL
  
  for(i in c(1:nrow(df))) {
    gene_list <- strsplit(as.character(df$name[i]), ";")[[1]]
    boolean <- gene_list %in% key_genes
    if(any(boolean) == TRUE) {
      df$is_key_gene[i] = "yes"
    } else { 
      df$is_key_gene[i] = "no"
    }
  }
  return(df)
}

test_res_VEHvsDHT_LNCAP <- matching_genes(df = test_res_VEHvsDHT_LNCAP, key_genes = TF_Chromatin_genes)
test_res_VEHvsDHT_LNCAP_LFQ <- matching_genes(df = test_res_VEHvsDHT_LNCAP_LFQ, key_genes = TF_Chromatin_genes)
test_res_VEHvsDHT_ZR751 <- matching_genes(df = test_res_VEHvsDHT_ZR751, key_genes = TF_Chromatin_genes)
test_res_VEHvsDHT_ZR751_LFQ <- matching_genes(df = test_res_VEHvsDHT_ZR751_LFQ, key_genes = TF_Chromatin_genes)

test_res_LNCAP_VEH_vs._ZR751_VEH <- matching_genes(df = test_res_LNCAP_VEH_vs._ZR751_VEH, key_genes = TF_Chromatin_genes)
test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ <- matching_genes(df = test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ, key_genes = TF_Chromatin_genes)
test_res_LNCAP_DHT_vs._ZR751_DHT <- matching_genes(df = test_res_LNCAP_DHT_vs._ZR751_DHT, key_genes = TF_Chromatin_genes)
test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ <- matching_genes(df = test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ, key_genes = TF_Chromatin_genes)

#####

write.csv(test_res_VEHvsDHT_LNCAP, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_LNCAP.csv", row.names = F)
write.csv(test_res_VEHvsDHT_LNCAP_LFQ, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_LNCAP_LFQ.csv", row.names = F)
write.csv(test_res_VEHvsDHT_ZR751, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_ZR751.csv", row.names = F)
write.csv(test_res_VEHvsDHT_ZR751_LFQ, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_ZR751_LFQ.csv", row.names = F)

write.csv(test_res_LNCAP_VEH_vs._ZR751_VEH, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_VEH_vs._ZR751_VEH.csv", row.names = F)
write.csv(test_res_LNCAP_DHT_vs._ZR751_DHT, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_DHT_vs._ZR751_DHT.csv", row.names = F)
write.csv(test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ.csv", row.names = F)
write.csv(test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ.csv", row.names = F)




test_res_VEHvsDHT_LNCAP$label[test_res_VEHvsDHT_LNCAP$is_key_gene != "no"] <- test_res_VEHvsDHT_LNCAP$name[test_res_VEHvsDHT_LNCAP$is_key_gene != "no"]
test_res_VEHvsDHT_LNCAP_LFQ$label[test_res_VEHvsDHT_LNCAP_LFQ$is_key_gene != "no"] <- test_res_VEHvsDHT_LNCAP_LFQ$name[test_res_VEHvsDHT_LNCAP_LFQ$is_key_gene != "no"]
test_res_VEHvsDHT_ZR751$label[test_res_VEHvsDHT_ZR751$is_key_gene != "no"] <- test_res_VEHvsDHT_ZR751$name[test_res_VEHvsDHT_ZR751$is_key_gene != "no"]
test_res_VEHvsDHT_ZR751_LFQ$label[test_res_VEHvsDHT_ZR751_LFQ$is_key_gene != "no"] <- test_res_VEHvsDHT_ZR751_LFQ$name[test_res_VEHvsDHT_ZR751_LFQ$is_key_gene != "no"]

test_res_LNCAP_VEH_vs._ZR751_VEH$label[test_res_LNCAP_VEH_vs._ZR751_VEH$is_key_gene != "no"] <- test_res_LNCAP_VEH_vs._ZR751_VEH$name[test_res_LNCAP_VEH_vs._ZR751_VEH$is_key_gene != "no"]
test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$label[test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$is_key_gene != "no"] <- test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$name[test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$is_key_gene != "no"]
test_res_LNCAP_DHT_vs._ZR751_DHT$label[test_res_LNCAP_DHT_vs._ZR751_DHT$is_key_gene != "no"] <- test_res_LNCAP_DHT_vs._ZR751_DHT$name[test_res_LNCAP_DHT_vs._ZR751_DHT$is_key_gene != "no"]
test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$label[test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$is_key_gene != "no"] <- test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$name[test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$is_key_gene != "no"]


library(ggplot2)
library(ggrepel)

ggplot(data= test_res_VEHvsDHT_LNCAP, aes(x= diff, y= -log10(pval), label = label)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_VEHvsDHT_LNCAP_LFQ, aes(x= diff, y= -log10(pval), label = label)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_VEHvsDHT_ZR751, aes(x= diff, y= -log10(pval), label = label)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_VEHvsDHT_ZR751_LFQ, aes(x= diff, y= -log10(pval), label = label)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_LNCAP_VEH_vs._ZR751_VEH, aes(x= diff, y= -log10(pval), label = label)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ, aes(x= diff, y= -log10(pval), label = label)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_LNCAP_DHT_vs._ZR751_DHT, aes(x= diff, y= -log10(pval), label = label)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ, aes(x= diff, y= -log10(pval), label = label)) + geom_point() + theme_minimal() + geom_text_repel()

#____________________________________________________________________________________________________________________________________________________________________
#Plotting volcano plots with just TF genes 

test_res_VEHvsDHT_LNCAP <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_LNCAP.csv", header = T)
test_res_VEHvsDHT_LNCAP_LFQ <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_LNCAP_LFQ.csv", header = T)
test_res_VEHvsDHT_ZR751 <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_ZR751.csv", header = T)
test_res_VEHvsDHT_ZR751_LFQ <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_ZR751_LFQ.csv", header = T)

test_res_LNCAP_VEH_vs._ZR751_VEH <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_VEH_vs._ZR751_VEH.csv")
test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ.csv")
test_res_LNCAP_DHT_vs._ZR751_DHT <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_DHT_vs._ZR751_DHT.csv")
test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ.csv")

test_res_VEHvsDHT_LNCAP <- test_res_VEHvsDHT_LNCAP[which(test_res_VEHvsDHT_LNCAP$is_key_gene == "yes"), ]
test_res_VEHvsDHT_LNCAP_LFQ <- test_res_VEHvsDHT_LNCAP_LFQ[which(test_res_VEHvsDHT_LNCAP_LFQ$is_key_gene == "yes"), ]
test_res_VEHvsDHT_ZR751 <- test_res_VEHvsDHT_ZR751[which(test_res_VEHvsDHT_ZR751$is_key_gene == "yes"), ]
test_res_VEHvsDHT_ZR751_LFQ <- test_res_VEHvsDHT_ZR751_LFQ[which(test_res_VEHvsDHT_ZR751_LFQ$is_key_gene == "yes"), ]

test_res_LNCAP_VEH_vs._ZR751_VEH <- test_res_LNCAP_VEH_vs._ZR751_VEH[which(test_res_LNCAP_VEH_vs._ZR751_VEH$is_key_gene == "yes"), ]
test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ <- test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ[which(test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$is_key_gene == "yes"), ]
test_res_LNCAP_DHT_vs._ZR751_DHT <- test_res_LNCAP_DHT_vs._ZR751_DHT[which(test_res_LNCAP_DHT_vs._ZR751_DHT$is_key_gene == "yes"), ]
test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ <- test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ[which(test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$is_key_gene == "yes"), ]

library(ggplot2)
library(ggrepel)

ggplot(data= test_res_VEHvsDHT_LNCAP, aes(x= diff, y= -log10(pval), label = test_res_VEHvsDHT_LNCAP$name)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_VEHvsDHT_LNCAP_LFQ, aes(x= diff, y= -log10(pval), label = test_res_VEHvsDHT_LNCAP_LFQ$name)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_VEHvsDHT_ZR751, aes(x= diff, y= -log10(pval), label = test_res_VEHvsDHT_ZR751$name)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_VEHvsDHT_ZR751_LFQ, aes(x= diff, y= -log10(pval), label = test_res_VEHvsDHT_ZR751_LFQ$name)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_LNCAP_VEH_vs._ZR751_VEH, aes(x= diff, y= -log10(pval), label = test_res_LNCAP_VEH_vs._ZR751_VEH$name)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ, aes(x= diff, y= -log10(pval), label = test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$name)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_LNCAP_DHT_vs._ZR751_DHT, aes(x= diff, y= -log10(pval), label = test_res_LNCAP_DHT_vs._ZR751_DHT$name)) + geom_point() + theme_minimal() + geom_text_repel()

ggplot(data= test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ, aes(x= diff, y= -log10(pval), label = test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$name)) + geom_point() + theme_minimal() + geom_text_repel()

#____________________________________________________________________________________________________________________________________________________________________
#Plotting volcano plots with TF genes labeled and top genes labeled 

test_res_VEHvsDHT_LNCAP <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_LNCAP.csv", header = T)
test_res_VEHvsDHT_LNCAP_LFQ <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_LNCAP_LFQ.csv", header = T)
test_res_VEHvsDHT_ZR751 <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_ZR751.csv", header = T)
test_res_VEHvsDHT_ZR751_LFQ <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_VEHvsDHT_ZR751_LFQ.csv", header = T)

test_res_LNCAP_VEH_vs._ZR751_VEH <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_VEH_vs._ZR751_VEH.csv")
test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ.csv")
test_res_LNCAP_DHT_vs._ZR751_DHT <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_DHT_vs._ZR751_DHT.csv")
test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ <- read.csv(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ.csv")

#

test_res_VEHvsDHT_LNCAP$diffexpressed <- "NO"
test_res_VEHvsDHT_LNCAP$diffexpressed[test_res_VEHvsDHT_LNCAP$diff > 2 & test_res_VEHvsDHT_LNCAP$pval < 0.05] <- "UP"
test_res_VEHvsDHT_LNCAP$diffexpressed[test_res_VEHvsDHT_LNCAP$diff < -2 & test_res_VEHvsDHT_LNCAP$pval < 0.05] <- "DOWN"

test_res_VEHvsDHT_LNCAP$label <- NA
test_res_VEHvsDHT_LNCAP$label[test_res_VEHvsDHT_LNCAP$diffexpressed != "NO"] <- test_res_VEHvsDHT_LNCAP$name[test_res_VEHvsDHT_LNCAP$diffexpressed != "NO"]
test_res_VEHvsDHT_LNCAP$label[test_res_VEHvsDHT_LNCAP$is_key_gene != "no"] <- test_res_VEHvsDHT_LNCAP$name[test_res_VEHvsDHT_LNCAP$is_key_gene != "no"]

#

test_res_VEHvsDHT_LNCAP_LFQ$diffexpressed <- "NO"
test_res_VEHvsDHT_LNCAP_LFQ$diffexpressed[test_res_VEHvsDHT_LNCAP_LFQ$diff > 2 & test_res_VEHvsDHT_LNCAP_LFQ$pval < 0.05] <- "UP"
test_res_VEHvsDHT_LNCAP_LFQ$diffexpressed[test_res_VEHvsDHT_LNCAP_LFQ$diff < -2 & test_res_VEHvsDHT_LNCAP_LFQ$pval < 0.05] <- "DOWN"

test_res_VEHvsDHT_LNCAP_LFQ$label <- NA
test_res_VEHvsDHT_LNCAP_LFQ$label[test_res_VEHvsDHT_LNCAP_LFQ$diffexpressed != "NO"] <- test_res_VEHvsDHT_LNCAP_LFQ$name[test_res_VEHvsDHT_LNCAP_LFQ$diffexpressed != "NO"]
test_res_VEHvsDHT_LNCAP_LFQ$label[test_res_VEHvsDHT_LNCAP_LFQ$is_key_gene != "no"] <- test_res_VEHvsDHT_LNCAP_LFQ$name[test_res_VEHvsDHT_LNCAP_LFQ$is_key_gene != "no"]

#

test_res_VEHvsDHT_ZR751$diffexpressed <- "NO"
test_res_VEHvsDHT_ZR751$diffexpressed[test_res_VEHvsDHT_ZR751$diff > 2 & test_res_VEHvsDHT_ZR751$pval < 0.05] <- "UP"
test_res_VEHvsDHT_ZR751$diffexpressed[test_res_VEHvsDHT_ZR751$diff < -2 & test_res_VEHvsDHT_ZR751$pval < 0.05] <- "DOWN"

test_res_VEHvsDHT_ZR751$label <- NA
test_res_VEHvsDHT_ZR751$label[test_res_VEHvsDHT_ZR751$diffexpressed != "NO"] <- test_res_VEHvsDHT_ZR751$name[test_res_VEHvsDHT_ZR751$diffexpressed != "NO"]
test_res_VEHvsDHT_ZR751$label[test_res_VEHvsDHT_ZR751$is_key_gene != "no"] <- test_res_VEHvsDHT_ZR751$name[test_res_VEHvsDHT_ZR751$is_key_gene != "no"]

#

test_res_VEHvsDHT_ZR751_LFQ$diffexpressed <- "NO"
test_res_VEHvsDHT_ZR751_LFQ$diffexpressed[test_res_VEHvsDHT_ZR751_LFQ$diff > 2 & test_res_VEHvsDHT_ZR751_LFQ$pval < 0.05] <- "UP"
test_res_VEHvsDHT_ZR751_LFQ$diffexpressed[test_res_VEHvsDHT_ZR751_LFQ$diff < -2 & test_res_VEHvsDHT_ZR751_LFQ$pval < 0.05] <- "DOWN"

test_res_VEHvsDHT_ZR751_LFQ$label <- NA
test_res_VEHvsDHT_ZR751_LFQ$label[test_res_VEHvsDHT_ZR751_LFQ$diffexpressed != "NO"] <- test_res_VEHvsDHT_ZR751_LFQ$name[test_res_VEHvsDHT_ZR751_LFQ$diffexpressed != "NO"]
test_res_VEHvsDHT_ZR751_LFQ$label[test_res_VEHvsDHT_ZR751_LFQ$is_key_gene != "no"] <- test_res_VEHvsDHT_ZR751_LFQ$name[test_res_VEHvsDHT_ZR751_LFQ$is_key_gene != "no"]

#

test_res_LNCAP_VEH_vs._ZR751_VEH$diffexpressed <- "NO"
test_res_LNCAP_VEH_vs._ZR751_VEH$diffexpressed[test_res_LNCAP_VEH_vs._ZR751_VEH$diff > 5 & test_res_LNCAP_VEH_vs._ZR751_VEH$pval < 0.05] <- "UP"
test_res_LNCAP_VEH_vs._ZR751_VEH$diffexpressed[test_res_LNCAP_VEH_vs._ZR751_VEH$diff < -5 & test_res_LNCAP_VEH_vs._ZR751_VEH$pval < 0.05] <- "DOWN"

test_res_LNCAP_VEH_vs._ZR751_VEH$label <- NA
test_res_LNCAP_VEH_vs._ZR751_VEH$label[test_res_LNCAP_VEH_vs._ZR751_VEH$diffexpressed != "NO"] <- test_res_LNCAP_VEH_vs._ZR751_VEH$name[test_res_LNCAP_VEH_vs._ZR751_VEH$diffexpressed != "NO"]
test_res_LNCAP_VEH_vs._ZR751_VEH$label[test_res_LNCAP_VEH_vs._ZR751_VEH$is_key_gene != "no"] <- test_res_LNCAP_VEH_vs._ZR751_VEH$name[test_res_LNCAP_VEH_vs._ZR751_VEH$is_key_gene != "no"]

#

test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$diffexpressed <- "NO"
test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$diffexpressed[test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$diff > 4 & test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$pval < 0.05] <- "UP"
test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$diffexpressed[test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$diff < -4 & test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$pval < 0.05] <- "DOWN"

test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$label <- NA
test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$label[test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$diffexpressed != "NO"] <- test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$name[test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$diffexpressed != "NO"]
test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$label[test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$is_key_gene != "no"] <- test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$name[test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ$is_key_gene != "no"]

#

test_res_LNCAP_DHT_vs._ZR751_DHT$diffexpressed <- "NO"
test_res_LNCAP_DHT_vs._ZR751_DHT$diffexpressed[test_res_LNCAP_DHT_vs._ZR751_DHT$diff > 5 & test_res_LNCAP_DHT_vs._ZR751_DHT$pval < 0.05] <- "UP"
test_res_LNCAP_DHT_vs._ZR751_DHT$diffexpressed[test_res_LNCAP_DHT_vs._ZR751_DHT$diff < -5 & test_res_LNCAP_DHT_vs._ZR751_DHT$pval < 0.05] <- "DOWN"

test_res_LNCAP_DHT_vs._ZR751_DHT$label <- NA
test_res_LNCAP_DHT_vs._ZR751_DHT$label[test_res_LNCAP_DHT_vs._ZR751_DHT$diffexpressed != "NO"] <- test_res_LNCAP_DHT_vs._ZR751_DHT$name[test_res_LNCAP_DHT_vs._ZR751_DHT$diffexpressed != "NO"]
test_res_LNCAP_DHT_vs._ZR751_DHT$label[test_res_LNCAP_DHT_vs._ZR751_DHT$is_key_gene != "no"] <- test_res_LNCAP_DHT_vs._ZR751_DHT$name[test_res_LNCAP_DHT_vs._ZR751_DHT$is_key_gene != "no"]

#

test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$diffexpressed <- "NO"
test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$diffexpressed[test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$diff > 3 & test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$pval < 0.05] <- "UP"
test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$diffexpressed[test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$diff < -4 & test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$pval < 0.05] <- "DOWN"

test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$label <- NA
test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$label[test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$diffexpressed != "NO"] <- test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$name[test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$diffexpressed != "NO"]
test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$label[test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$is_key_gene != "no"] <- test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$name[test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ$is_key_gene != "no"]

#

library(ggplot2)
library(ggrepel)

ggplot(data= test_res_VEHvsDHT_LNCAP, aes(x= diff, y= -log10(pval), label = label, col=diffexpressed)) + geom_point() + theme_minimal() + geom_text_repel() + scale_color_manual(values=c("blue", "black", "red")) + geom_vline(xintercept=c(-2, 2), col="red") + geom_hline(yintercept= -log10(0.05), col="red")

ggplot(data= test_res_VEHvsDHT_LNCAP_LFQ, aes(x= diff, y= -log10(pval), label = label, col=diffexpressed)) + geom_point() + theme_minimal() + geom_text_repel() + scale_color_manual(values=c("blue", "black", "red")) + geom_vline(xintercept=c(-2, 2), col="red") + geom_hline(yintercept= -log10(0.05), col="red")

ggplot(data= test_res_VEHvsDHT_ZR751, aes(x= diff, y= -log10(pval), label = label, col=diffexpressed)) + geom_point() + theme_minimal() + geom_text_repel() + scale_color_manual(values=c("blue", "black", "red")) + geom_vline(xintercept=c(-2, 2), col="red") + geom_hline(yintercept= -log10(0.05), col="red")

ggplot(data= test_res_VEHvsDHT_ZR751_LFQ, aes(x= diff, y= -log10(pval), label = label, col=diffexpressed)) + geom_point() + theme_minimal() + geom_text_repel() + scale_color_manual(values=c("blue", "black", "red")) + geom_vline(xintercept=c(-2, 2), col="red") + geom_hline(yintercept= -log10(0.05), col="red")

ggplot(data= test_res_LNCAP_VEH_vs._ZR751_VEH, aes(x= diff, y= -log10(pval), label = label, col=diffexpressed)) + geom_point() + theme_minimal() + geom_text_repel() + scale_color_manual(values=c("blue", "black", "red")) + geom_vline(xintercept=c(-5, 5), col="red") + geom_hline(yintercept= -log10(0.05), col="red")

ggplot(data= test_res_LNCAP_VEH_vs._ZR751_VEH_LFQ, aes(x= diff, y= -log10(pval), label = label,col=diffexpressed)) + geom_point() + theme_minimal() + geom_text_repel() + scale_color_manual(values=c("blue", "black", "red")) + geom_vline(xintercept=c(-4, 4), col="red") + geom_hline(yintercept= -log10(0.05), col="red")

ggplot(data= test_res_LNCAP_DHT_vs._ZR751_DHT, aes(x= diff, y= -log10(pval), label = label, col=diffexpressed)) + geom_point() + theme_minimal() + geom_text_repel() + scale_color_manual(values=c("blue", "black", "red")) + geom_vline(xintercept=c(-5, 5), col="red") + geom_hline(yintercept= -log10(0.05), col="red")

ggplot(data= test_res_LNCAP_DHT_vs._ZR751_DHT_LFQ, aes(x= diff, y= -log10(pval), label = label, col=diffexpressed)) + geom_point() + theme_minimal() + geom_text_repel() + scale_color_manual(values=c("blue", "black", "red")) + geom_vline(xintercept=c(-4, 3), col="red") + geom_hline(yintercept= -log10(0.05), col="red")

#_____________________________________________________________________________________________________________________
# Comparing results to Zwart AR RIME Paper

library(proDA)

fit_LNCAP_LFQ <- readRDS("C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_LNCAP_LFQ")

test_res_LNCAP_DHT_vs._IgG_LFQ <- test_diff(fit_LNCAP_LFQ, contrast = "conditionAR_DHT_LNCAP")

test_res_LNCAP_DHT_vs._IgG_LFQ <- test_res_LNCAP_DHT_vs._IgG_LFQ[order(test_res_LNCAP_DHT_vs._IgG_LFQ$pval), ]

Zwart_proteins <- read.xlsx(file = "~/Gap Year/CEDAR Student Worker Position/RIME Project/Zwart RIME Paper/41388_2018_BFonc2017330_MOESM364_ESM.xlsx", sheetName = "GeneList", header = F)
Zwart_proteins <- Zwart_proteins$X1

overlap <- intersect(Zwart_proteins, test_res_LNCAP_DHT_vs._IgG_LFQ$name[1:100])
cat(overlap, sep = "\n")
not_overlap_Zwart <- Zwart_proteins[!Zwart_proteins %in% overlap]
not_overlap_Arjun <- test_res_LNCAP_DHT_vs._IgG_LFQ$name[1:100][!test_res_LNCAP_DHT_vs._IgG_LFQ$name[1:100] %in% overlap]

results_table <- cbind(overlap, not_overlap_Zwart, not_overlap_Arjun)
results_table[, 2][54:80] <- "N/A"
results_table[, 1][21:80] <- "N/A"

colnames(results_table) <- c("Overlap", "Zwart RIME Paper", "Arjun Analysis")

#write.xlsx(results_table, file = "~/Gap Year/CEDAR Student Worker Position/RIME Project/results_table.xlsx", row.names = F, col.names = T)

##
library(proDA)
library(ggplot2)
library(ggrepel)

fit_LNCAP_LFQ <- readRDS("C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/RIME Project/fit_LNCAP_LFQ")

test_res_LNCAP_DHT_vs._IgG_LFQ <- test_diff(fit_LNCAP_LFQ, contrast = "conditionAR_DHT_LNCAP")

test_res_LNCAP_DHT_vs._IgG_LFQ$diffexpressed <- "NO"
test_res_LNCAP_DHT_vs._IgG_LFQ$diffexpressed[test_res_LNCAP_DHT_vs._IgG_LFQ$diff > 1.5 & test_res_LNCAP_DHT_vs._IgG_LFQ$pval < 0.05] <- "Enriched over IgG"

test_res_LNCAP_DHT_vs._IgG_LFQ$label <- NA
test_res_LNCAP_DHT_vs._IgG_LFQ$label[test_res_LNCAP_DHT_vs._IgG_LFQ$diffexpressed != "NO"] <- test_res_LNCAP_DHT_vs._IgG_LFQ$name[test_res_LNCAP_DHT_vs._IgG_LFQ$diffexpressed != "NO"]

ggplot(data= test_res_LNCAP_DHT_vs._IgG_LFQ, aes(x= diff, y= -log10(pval), label = label, col=diffexpressed)) + geom_point() + theme_minimal() + geom_text_repel() + scale_color_manual(values=c("blue", "black"))

#######################################################################################################
# Heatmap of RIME results

library(proDA)
library(pheatmap)
library(xlsx)

## Proteins restricted to TFs and chromatin modifiers 

normalized_abundance_matrix <- readRDS(file = "~/Gap Year/CEDAR Student Worker Position/RIME Project/normalized_abundance_matrix")

TFs_chromatin_modifiers <- read.xlsx(file = "~/Gap Year/CEDAR Student Worker Position/RIME Project/List of TF and Chromatin modifiers.xlsx", sheetName = "Combined list", header = F)

TFs_chromatin_modifiers <- TFs_chromatin_modifiers$X1

matching_genes_matrix <- function(matrix, key_genes) {
  library(dplyr)
  gene_match <- vector()
  
  for(i in c(1:nrow(matrix))) {
    gene_list <- strsplit(as.character(rownames(matrix)[i]), ";")[[1]]
    boolean <- gene_list %in% key_genes
    if(any(boolean) == TRUE) {
      gene_match[i] = TRUE
    } else { 
      gene_match[i] = FALSE
    }
  }
  new_matrix <- cbind(matrix, gene_match)
  return(new_matrix)
}

plot_mat <- matching_genes_matrix(normalized_abundance_matrix, TFs_chromatin_modifiers)
plot_mat <- plot_mat[plot_mat[, "gene_match"] == 1, ]
plot_mat <- plot_mat[, colnames(plot_mat) != "gene_match"]

### pheatmap(plot_mat)
### Too many missing values so had to remove rows with more than 15 missing values (>75% of all samples)

plot_mat_subsetted <- plot_mat[!rowSums(is.na(plot_mat)) > 15,]
pheatmap(plot_mat_subsetted)

## All proteins (2114), but had to remove rows with more than 10 missing values so down to 1180 proteins
### Heatmap too messy since so many rows

normalized_abundance_matrix_subsetted <- normalized_abundance_matrix[!rowSums(is.na(normalized_abundance_matrix)) > 10,]
pheatmap(normalized_abundance_matrix_subsetted)


