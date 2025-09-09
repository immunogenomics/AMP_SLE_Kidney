suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(pheatmap)
})

# Load the data
meta <- readRDS("../reorderedClusters/myeloid_tissue_meta_reOrdered.Rds")
norm <- readRDS("../finalObjects/myeloid_tissue_norm.Rds")

# Genes 
genes <- c("TREM2", "CD9", "GPNMB", "SPP1", "FABP5")

subset <- norm[genes, ]

subset <- as.matrix(t(subset))
subset_df <- data.frame(subset, cluster = meta$final_annotation) 
subset_df <- aggregate(. ~ cluster, data = subset_df, FUN = mean)

macrophages <- subset(subset_df, grepl("Macrophage", cluster))

plot_m <- as.matrix(macrophages[, 2:6])

rownames(plot_m) <- macrophages$cluster
colnames(plot_m) <- genes

pdf("SAM_heatmap.pdf", width = 9, height = 6)
    print(pheatmap(plot_m, 
            cluster_rows = TRUE, 
            cluster_cols = FALSE, 
            show_rownames = TRUE, 
            show_colnames = TRUE, 
            scale = "column",
            fontsize_row = 10, 
            fontsize_col = 10,
            main = "Scaled Expression of Selected Genes in Macrophages",
            color = colorRampPalette(c("blue", "white", "red"))(100),
            border_color = NA))
dev.off()