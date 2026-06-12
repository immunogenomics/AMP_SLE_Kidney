library(SeuratObject)
library(Matrix)
library(singlecellmethods)
library(dplyr)
library(tidyverse)
library(harmony)
library(pals)
library(ggplot2)
library(caret)
library(pROC)
library(ggrepel)
library(glmnet)

set.seed(1)

load("../cleaneddata.RData")

# Normalize data
myeloidCells <- names(clust)[clust %in% c("macrophage", "monocyte", "Mesangial.cell", "mDC", "mast", "neturophil", "pDC")]

countsNorm <- raw[, colnames(raw) %in% myeloidCells] %>%
  singlecellmethods::normalizeData(method = "log")

# Read in AMP data
ampNorm <- readRDS("../../../finalObjects_v4/myeloid_tissue_norm.Rds")
ampMeta <- readRDS("../../../finalObjects_v4/myeloid_tissue_meta.Rds")
ampMeta <- ampMeta[ampMeta$dataset == "scRNAseq", ]

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

# Calculate gene sparsity
gene_sparsity <- as.data.frame(rowSums(ampNorm == 0) / ncol(ampNorm))
colnames(gene_sparsity) <- "Sparsity"
gene_sparsity$gene <- rownames(gene_sparsity)

#Take non sparse genes
nonSparseGenes <- gene_sparsity %>%
  filter(Sparsity < 0.75) %>%
  pull(gene)

nonSparseGenes <- nonSparseGenes[nonSparseGenes %in% rownames(countsNorm)]

ampNormSub <- ampNorm[nonSparseGenes, colnames(ampNorm) %in% ampMetaSub$cell]
countsNorm <- countsNorm[nonSparseGenes, ]

set.seed(1)
samples <- unique(ampMetaSub$sample)
k <- 10

# Randomly assign samples to folds
sample_folds <- cut(seq_along(samples), breaks = k, labels = FALSE)
sample_folds <- sample(sample_folds)  # shuffle

fold_list <- lapply(1:k, function(i) {
  test_samples <- samples[sample_folds == i]
  test_cells <- which(ampMetaSub$sample %in% test_samples)
  return(test_cells)
})

results <- vector("numeric", k)

roc_list <- list()
auc_list <- numeric(k)

for(i in 1:k) {
  print(i)
  test_idx <- fold_list[[i]]
  train_idx <- setdiff(seq_len(ncol(ampNormSub)), test_idx)
  
  X_train <- t(ampNormSub)[train_idx, ]
  X_test  <- t(ampNormSub)[test_idx, ]
  y_train <- factor(ampMetaSub$SAMorNot)[train_idx]
  y_test  <- factor(ampMetaSub$SAMorNot)[test_idx]
  
  # Fit LASSO glmnet with internal CV for lambda
  cvfit <- cv.glmnet(
    x = X_train,
    y = y_train,
    family = "binomial",
    alpha = 1,
    type.measure = "class"
  )
  
  # Predict probabilities for "SAM"
  probs <- predict(cvfit, newx = X_test, s = "lambda.min", type = "response")
  
  roc_obj <- roc(as.numeric(y_test == "SAM"), probs)
  roc_list[[i]] <- roc_obj
  auc_list[i] <- auc(roc_obj)
}

# Step 2: Compute mean ROC
# pROC has a convenient "roc" function for this, using all test data together
all_true <- unlist(lapply(roc_list, function(r) r$labels))
all_probs <- unlist(lapply(roc_list, function(r) r$predictor))
mean_roc <- roc(all_true, all_probs)

# Step 3: Plot
# Save the ROC plot as a PDF
pdf("10fold_ROC_mean.pdf", width = 4.5, height = 4.5)

plot(0,0,type="n", xlim=c(1,0), ylim=c(0,1),
     xlab="Specificity", ylab="Sensitivity",
     main="10-fold ROC Curves with Mean ROC")
abline(a=1, b=-1, lty=2, col="gray")

colors <- rainbow(length(roc_list))
for(i in seq_along(roc_list)){
  lines(roc_list[[i]]$specificities,      # flip FPR -> Specificity
        roc_list[[i]]$sensitivities,
        col=colors[i], lwd=1.5)
}

# Mean ROC
lines(mean_roc$specificities, mean_roc$sensitivities,
      col="black", lwd=3)

legend("bottomright",
       legend=c(paste0("Fold ", 1:length(roc_list)), "Mean ROC"),
       col=c(colors,"black"),
       lwd=c(rep(1.5,length(roc_list)),3))

dev.off()  # Close the PDF device