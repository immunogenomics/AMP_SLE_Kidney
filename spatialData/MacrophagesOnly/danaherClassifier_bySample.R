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

load("cleaneddata.RData")
cannot <- read.csv("extraAnalyses/processed_data/cell positions wrt gloms.csv") # loads "cannot"
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

print(mergedDf[mergedDf$gene %in% c("TREM2", "CD9", "GPNMB", "FABP5", "CD63", "SPP1"), ])

#Take non sparse genes
nonSparseGenes <- gene_sparsity %>%
  filter(Sparsity < 0.75) %>%
  pull(gene)

nonSparseGenes <- nonSparseGenes[nonSparseGenes %in% rownames(countsNorm)]

countsNorm <- countsNorm[nonSparseGenes, ]

X_new <- t(countsNorm)

pred_probs <- predict(cvfit, newx = X_new, s = "lambda.min", type = "response")
pred_labels <- case_when(
  pred_probs[,1] > 0.8 ~ "SAM",
  pred_probs[,1] < 0.2 ~ "nonSAMs",
  TRUE ~ "ambiguous"
)

cannot$label <- pred_labels
cannot$prob_SAMs <- pred_probs[,1]
mergedCannot <- merge(cannot, annot, by.x = "id", by.y = "cell_ID")
mergedCannot$Disease <- ifelse(grepl("SLE", mergedCannot$tissuename), "SLE", "Control")

# Compute percentage for each combination of label, predicted_label, and location
mergedCannotBars <- mergedCannot %>%
  group_by(label, tissuename, position.vs.glom) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(label, tissuename) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ungroup()

mergedCannotBars <- mergedCannotBars %>%
  mutate(position.vs.glom = factor(position.vs.glom, levels = c("inside glomerulous", "bordering glomerulous", "tubulointerstitium")))  # bottom → top in stack

dir.create("bySample")
pdf("bySample/Spatial_distribution_barplot_allCells.pdf", width = 20, height = 12)
  ggplot(mergedCannotBars, aes(x = label, y = percent, fill = position.vs.glom)) +
  geom_bar(stat = "identity", position = "stack") +  # stacked by location
  geom_text(aes(label = paste0(count, " (", round(percent, 1), "%)")), 
          position = position_stack(vjust = 0.5), size = 3) +
  facet_wrap(~ tissuename, nrow = 2) +  # separates bars by predicted_label (MAC/SAMs)
  scale_fill_manual(values = c("inside glomerulous" = "#e01c1cff", "bordering glomerulous" = "#f8c12bff", "tubulointerstitium" = "#72d917ff")) + 
  labs(
    x = "Macrophage Type",
    y = "Percentage of cells",
    fill = "Location",
    title = "Localization of Macrophages and SAMs"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis text
dev.off()

mergedCannotBars <- mergedCannotBars[mergedCannotBars$label != "ambiguous", ]

pdf("bySample/Spatial_distribution_barplot_highConfidence.pdf", width = 20, height = 12)
  ggplot(mergedCannotBars, aes(x = label, y = percent, fill = position.vs.glom)) +
  geom_bar(stat = "identity", position = "stack") +  # stacked by location
  geom_text(aes(label = paste0(count, " (", round(percent, 1), "%)")), 
          position = position_stack(vjust = 0.5), size = 3) +
  facet_wrap(~ tissuename, nrow = 2) +  # separates bars by predicted_label (MAC/SAMs)
  scale_fill_manual(values = c("inside glomerulous" = "#e01c1cff", "bordering glomerulous" = "#f8c12bff", "tubulointerstitium" = "#72d917ff")) + 
  labs(
    x = "Macrophage Type",
    y = "Percentage of cells",
    fill = "Location",
    title = "Localization of Macrophages and SAMs"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis text
dev.off()

pdf("bySample/Spatial_distribution_barplot_highConfidence_byCounts.pdf", width = 20, height = 12)
  ggplot(mergedCannotBars, aes(x = label, y = count, fill = position.vs.glom)) +
  geom_bar(stat = "identity", position = "stack") +  # stacked by location
  geom_text(aes(label = paste0(count, " (", round(percent, 1), "%)")), 
          position = position_stack(vjust = 0.5), size = 3) +
  facet_wrap(~ tissuename, nrow = 2, scales = "free_y") +  # separates bars by predicted_label (MAC/SAMs)
  scale_fill_manual(values = c("inside glomerulous" = "#e01c1cff", "bordering glomerulous" = "#f8c12bff", "tubulointerstitium" = "#72d917ff")) + 
  labs(
    x = "Macrophage Type",
    y = "Percentage of cells",
    fill = "Location",
    title = "Localization of Macrophages and SAMs"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis text
dev.off()

# See if more or less confident predictions are in glomeruli
# Convert is_in_glom to numeric (0/1)
mergedCannot <- mergedCannot %>%
  mutate(
    is_in_glom_numeric = ifelse(position.vs.glom == "inside glomerulous", 1, 0)
  )

mergedCannot$prob_SAMs <- pred_probs[,1]
pdf("PredictionConfidence_vs_GlomerularLocation_allGenes.pdf", width = 8, height = 6)
    # Make a plot
    ggplot(mergedCannot, aes(x = prob_SAMs, y = is_in_glom_numeric, color = Disease)) +
    geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) + # logistic regression line
    labs(
        x = "Predicted probability of being SAMs",
        y = "Is in glomeruli (0/1)",
        color = "True label",
        title = "Relationship between prediction confidence and glomerular location"
    ) +
    theme_minimal()
dev.off()

saveRDS(obsPredict, file = "obs_with_predictions_nonSparse.Rds")

run_or <- function(freq_table) {
    output_df <- data.frame()
    for (i in seq(1, nrow(freq_table))) {
        pred_level <-  freq_table[i, "cell_type_pred_knn"]
        annot_level <- freq_table[i, "annotation"]
        iter_res <- or(pred_level, annot_level, freq_table)
        output_df <- output_df %>% rbind(iter_res)
    }
    colnames(output_df) <- c("cell_type_pred_knn", "annotation", "OR", "pvalue", "log_OR")
    return(output_df)
}
or <- function(pred_level, annot_level, freq_table) {
    inAandB <- freq_table %>% filter(annotation == annot_level, 
        cell_type_pred_knn == pred_level) %>% pull(Freq) %>% 
        sum()
    inAnotB <- freq_table %>% filter(annotation == annot_level, 
        cell_type_pred_knn != pred_level) %>% pull(Freq) %>% 
        sum()
    inBnotA <- freq_table %>% filter(annotation != annot_level, 
        cell_type_pred_knn == pred_level) %>% pull(Freq) %>% 
        sum()
    notBnotA <- freq_table %>% filter(annotation != annot_level, 
        cell_type_pred_knn != pred_level) %>% pull(Freq) %>% 
        sum()
    contin <- as.matrix(data.frame(c(inAandB, inAnotB), c(inBnotA, notBnotA)))

    fisher_res <- fisher.test(contin)
    chisq_res <- chisq.test(contin, simulate.p.value = TRUE, B = 100000)
    res_df <- t(data.frame(c(pred_level, 
                           annot_level, 
                           fisher_res$estimate, 
                           chisq_res$p.value, 
                           log(fisher_res$estimate))))
    return(res_df)
} 

freq_table <- as.data.frame.matrix(table(obsPredict$predicted_label, obsPredict$celltype_l1))

freq_table_long <- freq_table %>%
  pivot_longer(
    cols = everything(),           # Pivot all columns
    names_to = "annotation",       # Name for the new column with the original column names
    values_to = "Freq"             # Name for the new column with the original values
  ) %>%
  mutate(cell_type_pred_knn = rep(rownames(freq_table), each = ncol(freq_table))) %>%  # Add row names as a new column
  select(cell_type_pred_knn, annotation, Freq)  # Rearrange columns to have row_name, column_name, and value
freq_table_long <- as.data.frame(freq_table_long)

or_stats <- run_or(freq_table_long)

options(warn=-1)
p1 <- ggplot(or_stats, 
       aes(y = cell_type_pred_knn, 
           x = annotation, 
           fill = as.numeric(log_OR))) + 
   geom_tile(height = 0.8, width = 0.8) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                          na.value = "grey",
                          name = "log(Odds Ratio)", 
                         limits = c(-5, 5),
                         oob = scales::squish) +  
    theme_classic(base_size = 40) + 
    theme(legend.text = element_text(size = 30),
          legend.title = element_text(hjust = -0.5),
          axis.text.x = element_text(angle = 90)) +
    labs(x = "Original Cell States", y = "Transferred Cell States") + 
    theme(legend.key.size = unit(1.5, "cm"),
          legend.position = "right")
legend <- cowplot::get_legend(p1)
p1 <- p1 + theme(legend.position = "none")
legend <- cowplot::plot_grid(legend)
cowplot::save_plot("predLabels.png",
       p1,
       base_height = 20,
       base_width = 20)
cowplot::save_plot("predLabels_legend.png",
       legend,
       base_height = 5,
       base_width = 5)
