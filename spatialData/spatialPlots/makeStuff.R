library(Matrix)
library(singlecellmethods)
library(dplyr)
library(tidyverse)
library(pals)
library(ggplot2)
library(caret)
library(pROC)
library(ggrepel)
library(glmnet)
library(ggnewscale)

set.seed(1)

load("../cleaneddata.RData")
cannot <- read.csv("../extraAnalyses/processed_data/cell positions wrt gloms.csv") # loads "cannot"

# Normalize data
myeloidCells <- names(clust)[clust %in% c("macrophage", "monocyte", "Mesangial.cell", "mDC", "mast", "neturophil", "pDC")]
cannot <- cannot[cannot$id %in% myeloidCells, ]

countsNorm <- raw[, colnames(raw) %in% myeloidCells] %>%
  singlecellmethods::normalizeData(method = "log")

cvfit <- readRDS("danaherSAMclassifier.rds")

# Calculate gene sparsity
# Read in AMP data
ampNorm <- readRDS("../../finalObjects_v2/myeloid_tissue_norm.Rds")
ampMeta <- readRDS("../../finalObjects_v2/myeloid_tissue_meta.Rds")
ampMeta <- ampMeta[ampMeta$dataset == "scRNAseq" & ampMeta$sample != "AMPSLEkid_cells_1138", ]

# Assign meta clustering
ampMetaSub <- ampMeta %>%
  mutate(cell_group = case_when(
    grepl("Monocyte", final_annotation) ~ "MONO",
    final_annotation %in% c(
      "M5. GPNMBhigh NUPR1high Macrophage",
      "M7. SPP1high FABP5high Macrophage",
      "M9. MERTKhigh FABP5high Macrophage",
      "M11. GPMNBhigh NUPR1low Macrophage"
    ) ~ "SAMs",
    grepl("Macrophage", final_annotation) ~ "MAC",
    grepl("^DC", final_annotation) ~ "DC",
    grepl("MAST", final_annotation) ~ "MAST",
    grepl("^M15", final_annotation) ~ "PROLIF",
    TRUE ~ NA_character_
  )) %>%
  mutate(SAMorNot = ifelse(cell_group == "SAMs", "SAM", "nonSAMs"))

gene_sparsity <- as.data.frame(rowSums(ampNorm == 0) / ncol(ampNorm))
colnames(gene_sparsity) <- "Sparsity"
gene_sparsity$gene <- rownames(gene_sparsity)

# Plot
coefs <- coef(cvfit, s = "lambda.min")
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df$gene <- rownames(coefs_df)
top_genes <- coefs_df[order(abs(coefs_df$lambda.min), decreasing = TRUE), ]

mergedDf <- merge(top_genes, gene_sparsity, by = "gene")

# Order genes by beta (lambda.min)
mergedDf <- mergedDf %>%
  arrange(lambda.min) %>%
  mutate(gene = factor(gene, levels = gene))  # preserve order for plotting

#Take non sparse genes
nonSparseGenes <- gene_sparsity %>%
  filter(Sparsity < 0.75) %>%
  pull(gene)

nonSparseGenes <- nonSparseGenes[nonSparseGenes %in% rownames(countsNorm)]

countsNorm <- countsNorm[nonSparseGenes, ]

# Predict on Spatial data
# Predict class labels
trainIndex <- which(ampMetaSub$sample %in% sample(unique(ampMetaSub$sample), round(0.9 * length(unique(ampMetaSub$sample)))))
y_test  <- factor(ampMetaSub$SAMorNot)[-trainIndex]

X_new <- t(countsNorm)

pred_probs <- predict(cvfit, newx = X_new, s = "lambda.min", type = "response")
pred_labels <- case_when(
  pred_probs[,1] > 0.8 ~ "SAM",
  pred_probs[,1] < 0.2 ~ "nonSAMs",
  TRUE ~ "ambiguous"
)

saveRDS(pred_labels, "classifiedLabels.Rds")
saveRDS(colnames(countsNorm), "cellNames.Rds")