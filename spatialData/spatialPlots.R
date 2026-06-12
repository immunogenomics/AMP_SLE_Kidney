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

# See histogram of probabilities
# Extract probabilities as a vector
pred_probs_vec <- as.numeric(pred_probs[, 1])

# relod cannot
cannot <- read.csv("extraAnalyses/processed_data/cell positions wrt gloms.csv") # loads "cannot"

macs <- data.frame(cell = colnames(countsNorm), label = pred_labels)

viz_df <- data.frame(viz)
viz_df$id <- annot$cell_ID
fullSlides <- merge(cannot[, c("id", "position.vs.glom")], annot[, c("cell_ID", "tissuename", "fov")], by.x = "id", by.y = "cell_ID")
fullSlides <- merge(fullSlides, viz_df, by = "id")

samMarkers <- data.frame(data.matrix(t(countsNorm[c("CD63", "CD9", "FABP5", "GPNMB", "SPP1"), ])))

myeloidExpr <- merge(fullSlides, samMarkers, by.x = "id", by.y = "row.names")
myeloidExpr <- merge(myeloidExpr, macs, by.x = "id", by.y = "cell")
dir.create("spatialPlots", showWarnings = FALSE)
for (fov in unique(fullSlides$fov)) {
  print(fov)
  tissue <- unique(fullSlides$tissuename[fullSlides$fov == fov])
  dir.create(paste0("spatialPlots/", tissue), showWarnings = FALSE)
  for (marker in c("CD63", "CD9", "FABP5", "GPNMB", "SPP1")) {
    pdf(file = paste0("spatialPlots/", tissue, "/spatialPlots_", marker, "_", fov, ".pdf"),
            width = 8, height = 6)
        
        p <- ggplot() +
          # First layer: categorical color for position
          ggrastr::rasterise(
            geom_point(
              data = fullSlides[fullSlides$fov == fov, ],
              aes(x = sdimx, y = sdimy, color = position.vs.glom),
              alpha = 0.7, size = 3),
            dpi = 300
          ) +
          scale_color_manual(values = c(
            "tubulointerstitium" = "lightgrey",
            "bordering glomerulous" = "steelblue1",
            "inside glomerulous" = "tomato"
          )) +
          
          # Tell ggplot you want a *new* color scale for the next geom
          ggnewscale::new_scale_color() +
          
          # Second layer: continuous color for marker expression
          ggrastr::rasterise(
            geom_point(
              data = myeloidExpr[myeloidExpr$fov == fov, ] %>% arrange(.data[[marker]]),
              aes(x = sdimx, y = sdimy, color = .data[[marker]], shape = label),
              alpha = 0.9, size = 3
            ),
            dpi = 300
          ) +
          scale_color_viridis_c() +
          scale_shape_manual(values = c(
            nonSAMs = 15,    # square
            ambiguous = 16,  # circle
            SAM = 17         # triangle
          )) + 
          labs(title = paste("Tissue:", tissue)) +
          theme_minimal()
        
        print(p)
        dev.off()
  }
}

# Plot of called SAMs and nonSAMs in Glomeruli
samsOrNot <- merge(fullSlides, macs, by.x = "id", by.y = "cell")
