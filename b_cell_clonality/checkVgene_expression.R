suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
})

norm <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects_v3/b_tissue_norm.Rds")
meta <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects_v3/b_tissue_meta.Rds")

vgenes <- rownames(norm)[grep("^IG[HKL]V", rownames(norm))]

lightChain <- vgenes[grep("^IG[LK]V", vgenes)]  # 178 genes
heavyChain <- vgenes[grep("^IGHV", vgenes)]  # 149 genes

lExpr <- norm[lightChain, ]

# Check to see if any genes aren't expressed
nonZero_light <- lightChain[rowSums(lExpr) > 0]  # 112 genes
write.csv(as.data.frame(nonZero_light), "nonZero_lightChainGenes.csv")

# Histogram of total normalized light chain expression per cell
total_lExpr_perCell <- colSums(lExpr[nonZero_light, ])

plot_df <- data.frame(expression = total_lExpr_perCell, 
    annotation = meta$final_annotation, 
    broad = ifelse(grepl("^P", meta$final_annotation), "Plasma", ifelse(grepl("^B", meta$final_annotation), "B Cell", meta$final_annotation)))

pdf("lightChain_expression_histogram.pdf", width = 6, height = 4)
    ggplot(data.frame(plot_df), aes(x = expression)) +
        geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black") +
        labs(title = "Total Normalized Light Chain Expression per Cell",
            x = "Total Normalized Expression",
            y = "Number of Cells") +
        theme_minimal()
dev.off()

pdf("lightChain_expression_histogram_broad.pdf", width = 12, height = 4)
    ggplot(data.frame(plot_df), aes(x = expression)) +
        geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black") +
        labs(title = "Total Normalized Light Chain Expression per Cell (by Broad)",
            x = "Total Normalized Expression",
            y = "Number of Cells") +
        facet_wrap(~ broad, scales = "free_y") +
        theme_minimal()
dev.off()

pdf("lightChain_expression_histogram_annotation.pdf", width = 6, height = 4)
    ggplot(data.frame(plot_df), aes(x = expression)) +
        geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black") +
        labs(title = "Total Normalized Light Chain Expression per Cell (by Annotation)",
            x = "Total Normalized Expression",
            y = "Number of Cells") +
        facet_wrap(~ annotation, scales = "free_y") +
        theme_minimal()
dev.off()

# Heatmap of light chain gene expression across all cells
beatmap <- as.data.frame(as.matrix(t(lExpr[nonZero_light, colSums(lExpr[nonZero_light, ]) > 0])))

beatmapPlasma <- beatmap[row.names(beatmap) %in% meta$cell[grep("^P", meta$final_annotation)], ]
beatmapBcells <- beatmap[row.names(beatmap) %in% meta$cell[grep("^B", meta$final_annotation)], ]

# Function to Z-score rows (cells)
zscore_rows <- function(mat) {
  t(apply(mat, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
}

# Z-score across cells
beatmapPlasma_z <- zscore_rows(beatmapPlasma)
beatmapPlasma_z <- beatmapPlasma_z[, order(colnames(beatmapPlasma_z))]

beatmapBcells_z  <- zscore_rows(beatmapBcells)
beatmapBcells_z <- beatmapBcells_z[, order(colnames(beatmapBcells_z))]

clinData <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects_v2/clinicalData.rds")

cellChronicities <- meta %>%
  left_join(
    clinData %>% select(sample, Final_Chronicity),
    by = "sample"
  )

# Create a lookup: cell -> Final_Chronicity
cell_chronicity_lookup <- cellChronicities %>%
  select(cell = cell, Final_Chronicity) %>%
  filter(!is.na(Final_Chronicity))

# --- Plasma cells ---
# Keep only cells that are in meta
common_cells_plasma <- intersect(rownames(beatmapPlasma_z), cell_chronicity_lookup$cell)

# Order by Final_Chronicity
ordered_cells_plasma <- cell_chronicity_lookup %>%
  filter(cell %in% common_cells_plasma) %>%
  arrange(Final_Chronicity) %>%
  pull(cell)

# Subset and reorder
beatmapPlasma_z <- beatmapPlasma_z[ordered_cells_plasma, ]

# --- B cells ---
common_cells_b <- intersect(rownames(beatmapBcells_z), cell_chronicity_lookup$cell)

ordered_cells_b <- cell_chronicity_lookup %>%
  filter(cell %in% common_cells_b) %>%
  arrange(Final_Chronicity) %>%
  pull(cell)

beatmapBcells_z <- beatmapBcells_z[ordered_cells_b, ]

# Make a data frame: row names are cells (rows of your heatmap)
# --- Plasma cells ---
chronicity_annotation_plasma <- data.frame(
  Chronicity = cell_chronicity_lookup$Final_Chronicity
)
rownames(chronicity_annotation_plasma) <- cell_chronicity_lookup$cell
# Keep only rows present in the heatmap
chronicity_annotation_plasma <- chronicity_annotation_plasma[rownames(beatmapPlasma_z), , drop = FALSE]

# --- B cells ---
chronicity_annotation_b <- data.frame(
  Chronicity = cell_chronicity_lookup$Final_Chronicity
)
rownames(chronicity_annotation_b) <- cell_chronicity_lookup$cell
# Keep only rows present in the B cell heatmap
chronicity_annotation_b <- chronicity_annotation_b[rownames(beatmapBcells_z), , drop = FALSE]

# Heatmap of plasma cells
pdf("lightChain_expression_heatmap_plasmaCells.pdf", width = 14, height = 8)
    # Plasma cells heatmap
    pheatmap::pheatmap(
        beatmapPlasma_z,
        scale = "none",
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
        treeheight_row = 0, 
        treeheight_col = 0,
        show_rownames = FALSE,
        annotation_row = chronicity_annotation_plasma,
        annotation_colors = list(Chronicity = colorRampPalette(c("yellow", "red"))(100)),
        use_raster = TRUE,
        fontsize = 14,          # global text size
        fontsize_row = 14,      # gene labels
        fontsize_col = 14,       # sample labels (if shown)
        main = paste("Light Chain Gene Expression in", nrow(beatmapPlasma_z),"Plasma Cells; n =", length(unique(cellChronicities$sample[cellChronicities$cell %in% rownames(beatmapPlasma_z)])), "samples")
    )
dev.off()

# Heatmap of b cells
pdf("lightChain_expression_heatmap_bCells.pdf", width = 14, height = 8)
    pheatmap::pheatmap(
        beatmapBcells_z,
        scale = "none",
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
        treeheight_row = 0, 
        treeheight_col = 0,
        show_rownames = FALSE,
        show_colnames = FALSE,
        annotation_row = chronicity_annotation_b,
        use_raster = TRUE,
        fontsize = 14,          # global text size
        fontsize_row = 14,      # gene labels
        fontsize_col = 14,       # sample labels (if shown)
        annotation_colors = list(Chronicity = colorRampPalette(c("yellow", "red"))(100)),
        main = paste("Light Chain Gene Expression in", nrow(beatmapBcells_z),"Plasma Cells; n =", length(unique(cellChronicities$sample[cellChronicities$cell %in% rownames(beatmapBcells_z)])), "samples")
    )
dev.off()

# Heatmap of plasma cells
pdf("lightChain_expression_heatmap_plasmaCells_fullCluster.pdf", width = 14, height = 8)
    # Plasma cells heatmap
    pheatmap::pheatmap(
        beatmapPlasma_z,
        scale = "none",
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
        treeheight_row = 0, 
        treeheight_col = 0,
        show_rownames = FALSE,
        annotation_row = chronicity_annotation_plasma,
        annotation_colors = list(Chronicity = colorRampPalette(c("yellow", "red"))(100)),
        main = paste("Light Chain Gene Expression in", nrow(beatmapPlasma_z),"Plasma Cells; n =", length(unique(cellChronicities$sample[cellChronicities$cell %in% rownames(beatmapPlasma_z)])), "samples")
    )
dev.off()

# Heatmap of b cells
pdf("lightChain_expression_heatmap_bCells_fullCluster.pdf", width = 14, height = 8)
    pheatmap::pheatmap(
        beatmapBcells_z,
        scale = "none",
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
        treeheight_row = 0, 
        treeheight_col = 0,
        show_rownames = FALSE,
        annotation_row = chronicity_annotation_b,
        annotation_colors = list(Chronicity = colorRampPalette(c("yellow", "red"))(100)),
        main = paste("Light Chain Gene Expression in", nrow(beatmapBcells_z),"Plasma Cells; n =", length(unique(cellChronicities$sample[cellChronicities$cell %in% rownames(beatmapBcells_z)])), "samples")
    )
dev.off()

pseudobulk <- function(counts, groups) {
    colnames(counts) = groups
    
    # Convert group labels into a sparse indicator matrix
    group_levels <- unique(groups)
    group_matrix <- sparseMatrix(
      i = match(groups, group_levels), 
      j = seq_len(length(groups)), 
    )

    # Perform sparse matrix multiplication to sum within groups
    pseudo_raw <- group_matrix %*% (counts %>% t)

    # Assign row names to match groups
    rownames(pseudo_raw) <- group_levels

    pseudo_raw = pseudo_raw %>% t
    
    return(pseudo_raw)
    
}

plasmaMeta <- meta %>% filter(meta$cell %in% rownames(beatmapPlasma))
pb <- pseudobulk(t(as.matrix(beatmapPlasma[, c("IGKV3-11", "IGKV3D-11")])), plasmaMeta$sample)


clinDataFiltered <- clinData %>% filter(!(is.na(Final_Chronicity)))
# Keep only samples present in both pseudobulk and clinical data
common_samples <- intersect(colnames(pb), clinDataFiltered$sample)

pb_sub <- pb[, common_samples, drop = FALSE]

clin_sub <- clinDataFiltered %>%
  dplyr::filter(sample %in% common_samples) %>%
  dplyr::arrange(Final_Chronicity)

# Order pseudobulk columns by increasing chronicity
ordered_samples <- clin_sub$sample
pb_ordered <- pb_sub[, ordered_samples, drop = FALSE]

# Column annotation for chronicity
s <- data.frame(
  Final_Chronicity = clin_sub$Final_Chronicity
)

annotation_col <- data.frame(
  Final_Chronicity = clin_sub$Final_Chronicity
)
rownames(annotation_col) <- clin_sub$sample

# Save heatmap to PDF
pdf("plasma_cells_pb.pdf", width = 10, height = 3)

pheatmap::pheatmap(
  pb_ordered,
  scale = "row",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_cols = FALSE,   # preserve chronicity ordering
  cluster_rows = FALSE,
  annotation_col = annotation_col,
  annotation_colors = list(Final_Chronicity = colorRampPalette(c("yellow", "red"))(100)),
  show_colnames = FALSE,
  main = "Plasma Cell Pseudobulk IGKV Expression Ordered by Chronicity"
)

dev.off()
