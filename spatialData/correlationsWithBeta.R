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
library(wesanderson)

set.seed(1)

cvfit <- readRDS("danaherSAMclassifier.rds")
coefs <- coef(cvfit, s = "lambda.min")
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df$gene <- rownames(coefs_df)

# Load DE gene correlates
deGenes <- read.csv("all DE results/DE immune cells - glom vs tubulointerstitium - macrophage.csv")
mergedDf <- merge(coefs_df, deGenes, by.x = "gene", by.y = "target")

# Calculate correlation
cor_test <- cor.test(mergedDf$lambda.min[mergedDf$fdr < 0.05], mergedDf$log2foldchange[mergedDf$fdr < 0.05])
cor_val <- signif(cor_test$estimate, 3)
p_val <- signif(cor_test$p.value, 3)

# Prepare subtitle
subtitle_text <- paste0("Pearson r = ", cor_val, ", p = ", p_val)

# Plot correlations
pdf("BetaDECorrelations.pdf", width = 8, height = 6)
ggplot(mergedDf[mergedDf$fdr < 0.05, ], aes(x = log2foldchange, y = lambda.min)) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_point(alpha = 0.5) +
  geom_text_repel(
    aes(label = gene),
    size = 3,
    max.overlaps = 20  # adjust as needed
  ) +
  labs(title = "Correlation between DE logFC and Classifier Coefficients",
       x = "Log Fold Change (Glom vs Tubulointerstitium)",
       y = "SAM Beta",
       subtitle = subtitle_text) +
  theme_minimal()
dev.off()

# Create SAM score
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

# Get only relevant, significant genes
deGenes_sig <- deGenes %>%
  filter(fdr < 0.05) %>%
  select(target, log2foldchange, fdr)

# Macrophages Only
ampMetaMacs <- ampMetaSub %>%
  filter(grepl("Macrophage", final_annotation))

deGenes_sig <- deGenes_sig[deGenes_sig$target %in% rownames(ampNorm), ]

ampNormMacs <- ampNorm[rownames(ampNorm) %in% deGenes_sig$target, colnames(ampNorm) %in% ampMetaMacs$cell] %>% singlecellmethods::scaleData()
ampNormMacs <- ampNorm[deGenes_sig$target, ampMetaMacs$cell] # Reorer the rows

# Calculate SAM Score
# Compute SAM score: weighted sum across rows for each column (cell)
sam_score <- as.numeric(t(ampNormMacs) %*% deGenes_sig$log2foldchange)
names(sam_score) <- colnames(ampNormMacs)

# Attach scores to metadata (matched by cell name)
ampMetaMacs$InfilScore <- sam_score

saveRDS(ampMetaMacs, file = "MacrophagesOnly/macrophagesOnly_InfilScore.Rds")
# Visualize SAM scores
dir.create("MacrophagesOnly/")
pdf("MacrophagesOnly/SAMScores_Macrophages_broad.pdf", width = 8, height = 6)
    ggplot(ampMetaMacs, aes(x = cell_group, y = InfilScore, fill = SAMorNot)) +
        geom_violin() +
        geom_boxplot(width = 0.2) +
        labs(title = "SAM Scores in Macrophage Subtypes",
            x = "Cell Group",
            y = "InfilScore") +
        scale_fill_manual(values = c("nonSAMs" = "#0073c2", "SAM" = "#EFC00F")) +
        theme_minimal() + 
        geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) + 
        theme(legend.position = "none") # Remove legend
dev.off()

ampMetaMacs$shortCluster <- sub("\\..*", "", ampMetaMacs$final_annotation)

pdf("MacrophagesOnly/SAMScores_Macrophages_byCluster.pdf", width = 8, height = 6)
    ggplot(ampMetaMacs, aes(x = reorder(shortCluster, InfilScore, FUN = median), y = InfilScore, fill = SAMorNot)) +
        geom_violin() +
        geom_boxplot(width = 0.2) +
        labs(title = "SAM Scores in Macrophage Subtypes",
            x = "Cell Group",
            y = "InfilScore") +
        scale_fill_manual(values = c("nonSAMs" = "#0073c2", "SAM" = "#EFC00F")) +
        theme_minimal() +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) + 
        theme(legend.position = "none")  # Remove legend
dev.off()

# UMAP Plot
pdf("MacrophagesOnly/SAMScores_Macrophages_UMAP.pdf", width = 8, height = 6)
    ggplot() +
        ggrastr::rasterise(geom_point(data = ampMeta, aes(x = hUMAP1, y = hUMAP2), color = "lightgrey", size = 0.1), dpi = 300) +
        ggrastr::rasterise(geom_point(data = ampMetaMacs %>% arrange(abs(InfilScore)), aes(x = hUMAP1, y = hUMAP2, color = InfilScore), size = 0.1), dpi = 300) + 
        scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous")) + 
        labs(title = "UMAP of Macrophages Colored by InfilScore",
            x = "UMAP 1",
            y = "UMAP 2") +
        theme_classic()
dev.off()

# Take only SAMs. Correlate with clinical data later
ampMetaSAMs <- ampMetaMacs %>%
  filter(SAMorNot == "SAM")

clinicalData <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects_v2/clinicalData.rds")

samsPerSample <- ampMetaSAMs %>%
  group_by(sample) %>%
  summarise(medianInfilScore = median(InfilScore, na.rm = TRUE))

samsPerSample <- merge(samsPerSample, clinicalData, by = "sample") %>% filter(!is.na(Final_Chronicity))

# Correlations with Chronicity
correlations = cor.test(samsPerSample$medianInfilScore, samsPerSample$Final_Chronicity)
pdf("MacrophagesOnly/InfilScore_vs_Chronicity.pdf", width = 6, height = 6)
  ggplot(samsPerSample, aes(x = medianInfilScore, y = Final_Chronicity)) +
      geom_point() +
      geom_smooth(method = "lm", color = "plum3", se = FALSE) +
      labs(title = "Median InfilScore vs Chronicity Score",
          x = "Median InfilScore",
          y = "Chronicity") +
      theme_minimal() + 
      labs(title = "Sample InfilScore vs Chronicity", subtitle = paste("R = ", signif(correlations$estimate, 3), ", p = ", signif(correlations$p.value, 3)))
dev.off()

# Correlations with Activity
correlations = cor.test(samsPerSample$medianInfilScore, samsPerSample$Final_Activity)
pdf("MacrophagesOnly/InfilScore_vs_Activity.pdf", width = 6, height = 6)
  ggplot(samsPerSample, aes(x = medianInfilScore, y = Final_Activity)) +
      geom_point() +
      geom_smooth(method = "lm", color = "plum3", se = FALSE) +
      labs(title = "Median InfilScore vs Activity",
          x = "Median InfilScore",
          y = "Activity") +
      theme_minimal() + 
      labs(title = "Sample InfilScore vs Activity", subtitle = paste("R = ", signif(correlations$estimate, 3), ", p = ", signif(correlations$p.value, 3)))
dev.off()

# Correlations with ISN Class
pdf("MacrophagesOnly/InfilScore_vs_ISN.pdf", width = 6, height = 6)
  ggplot(samsPerSample, aes(x = Final_ISN, y = medianInfilScore)) +
      geom_boxplot(fill = "plum3") +
      labs(title = "Median InfilScore vs ISN Class",
          x = "ISN Class",
          y = "Median InfilScore") +
      theme_minimal()
dev.off()

# All Myeloid Cells
ampScaledAll <- ampNorm[rownames(ampNorm) %in% deGenes_sig$target, colnames(ampNorm) %in% ampMetaSub$cell] %>% singlecellmethods::scaleData()
ampScaledAll <- ampNorm[deGenes_sig$target, ampMetaSub$cell] # Reorer the rows

sam_score <- as.numeric(t(ampScaledAll) %*% deGenes_sig$log2foldchange)
names(sam_score) <- colnames(ampScaledAll)

ampMetaSub$InfilScore <- sam_score


# Visualize SAM scores
dir.create("AllMyeloid/")
pdf("AllMyeloid/SAMScores_Macrophages_broad.pdf", width = 10, height = 6)
    ggplot(ampMetaSub, aes(x = cell_group, y = InfilScore, fill = SAMorNot)) +
        geom_violin() +
        geom_boxplot(width = 0.2) +
        labs(title = "InfilScore in Macrophage Subtypes",
            x = "Cell Group",
            y = "InfilScore") +
        scale_fill_manual(values = c("nonSAMs" = "#0073c2", "SAM" = "#EFC00F")) +
        theme(legend.position = "none") + # Remove legend
        geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) + 
        theme_minimal()
dev.off()

ampMetaSub$shortCluster <- sub("\\..*", "", ampMetaSub$final_annotation)

pdf("AllMyeloid/SAMScores_Macrophages_byCluster.pdf", width = 16, height = 6)
    ggplot(ampMetaSub, aes(x = reorder(shortCluster, InfilScore, FUN = median), y = InfilScore, fill = cell_group)) +
        geom_violin() +
        geom_boxplot(width = 0.2) +
        labs(title = "InfilScore in Macrophage Subtypes",
            x = "Cell Group",
            y = "InfilScore") +
        theme_minimal() + 
        geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) + 
        theme(legend.position = "none")
dev.off()

# UMAP Plot
pdf("AllMyeloid/SAMScores_Macrophages_UMAP.pdf", width = 8, height = 6)
    ggplot() +
        ggrastr::rasterise(geom_point(data = ampMetaSub %>% arrange(abs(InfilScore)), aes(x = hUMAP1, y = hUMAP2, color = InfilScore), size = 0.1), dpi = 300) + 
        scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous")) + 
        labs(title = "UMAP of Macrophages Colored by InfilScore",
            x = "UMAP 1",
            y = "UMAP 2") +
        theme_classic()
dev.off()

saveRDS(ampMetaSub, file = "AllMyeloid/allMyeloid_InfilScore.Rds")