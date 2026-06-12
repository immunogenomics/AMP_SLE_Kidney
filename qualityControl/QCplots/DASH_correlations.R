library(Seurat)
library(singlecellmethods)

path_to_no_dash3 <- "/data/srlab2/jmears/SLE_SOP/Nuc-seq-Live-Comparison-2/191212_SL-NVQ_0129_BHMLFTDMXX/cellranger-3.1.0/GRCh38-1/LN-2837-No-Dash/outs/raw_feature_bc_matrix/"
noDash <- Read10X(data.dir = path_to_no_dash3)

# Count genes detected per cell
nUMIsPerCell <- colSums(noDash)

# Keep only cells with > 500 genes
noDash_filtered <- noDash[, nUMIsPerCell > 500] %>% normalizeData(method = "log")

# Identify mitochondrial genes
mito_genes <- grep("^MT-", rownames(noDash_filtered), value = TRUE)

# Specify genes
target_genes <- c("MT-ND5", "MT-ND6")

# Subset counts for target and mito genes
ND5_counts <- Matrix::colSums(noDash_filtered["MT-ND5", , drop = FALSE])
ND6_counts <- Matrix::colSums(noDash_filtered["MT-ND6", , drop = FALSE])

# All mitochondrial genes except ND5 and ND6
other_mito_genes <- setdiff(mito_genes, c("MT-ND5", "MT-ND6"))
other_mito_counts <- Matrix::colSums(noDash_filtered[other_mito_genes, , drop = FALSE])

# Total counts per cell
total_counts <- Matrix::colSums(noDash_filtered)

# Build dataframe with percentages
mito_df <- data.frame(
  cell = colnames(noDash_filtered),
  percent_ND5 = ND5_counts / total_counts * 100,
  percent_ND6 = ND6_counts / total_counts * 100,
  both = (ND5_counts + ND6_counts) / total_counts * 100,
  percent_other_mito = other_mito_counts / total_counts * 100
)

library(ggplot2)

# Function to make scatter plot and save
plot_vs_other <- function(gene_name, df, file_name) {
  
  # Compute slope of line
  fit <- lm(as.formula(paste("percent_other_mito ~", gene_name)), data = df)
  r_val <- cor(df[[gene_name]], df$percent_other_mito)

  # Create plot
  p <- ggplot(df, aes_string(x = gene_name, y = "percent_other_mito")) +
    ggrastr::rasterise(geom_point(alpha = 0.6, color = "steelblue"), dpi = 300) +
    geom_smooth(method = "lm", color = "darkred", se = FALSE) +
    labs(
      x = paste0(gene_name, " (%)"),
      y = "% of other Mitochondrial Genes",
      title = paste0(gene_name, " vs Other MT genes"),
      subtitle = paste0("Slope = ", round(coef(fit)[2], 3), "; R = ", round(r_val, 3))
    ) +
    theme_minimal(base_size = 16)
  
  # Save to PDF
  pdf(file_name, width = 4, height = 4)
  print(p)
  dev.off()
}

# Plot MT-ND4 vs Other_MT
plot_vs_other("percent_ND5", mito_df, "MT_ND5_vs_OtherMT.pdf")

# Plot MT-ND5 vs Other_MT
plot_vs_other("percent_ND6", mito_df, "MT_ND6_vs_OtherMT.pdf")

# Plot both vs Other_MT
plot_vs_other("both", mito_df, "Both_ND5_ND6_vs_OtherMT.pdf")