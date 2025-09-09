# Use with the sc environment
library(SeuratObject)
library(singlecellmethods)
library(dplyr)
library(purrr)

set.seed(0)

args <- commandArgs(trailingOnly=TRUE)

# Tissues
tissue_norm <- readRDS(as.character(args[1]))
tissue_meta <- readRDS(as.character(args[2]))

# PBMCs
pbmc_norm <- readRDS(as.character(args[3]))[['RNA']]@data
pbmc_meta <- readRDS(as.character(args[4]))

name <- as.character(args[6])

# Sample to unified visit
clinical <- readRDS('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/clinical_data_05042023.rds')
sampleToDonor <- data.frame(sample = clinical$sample, unified_visit = clinical$AMP.ID)

# Add in Unified_Visit column
tissue_meta <- merge(tissue_meta, sampleToDonor, by = "sample")

tissue_meta <- tissue_meta[tissue_meta$unified_visit %in% pbmc_meta$Unified_Visit, ]
pbmc_meta <- pbmc_meta[pbmc_meta$Unified_Visit %in% tissue_meta$unified_visit, ]

# Subset to take only shared cells
pbmc_norm <- pbmc_norm[, colnames(pbmc_norm) %in% pbmc_meta$Cell]
tissue_norm <- tissue_norm[, colnames(tissue_norm) %in% tissue_meta$cell]

# Calculate the weights of the combined references
y1 <- factor(c(tissue_meta$sample, pbmc_meta$Unified_Visit))
weights1 <- as.numeric((1 / prop.table(table(y1)))[y1]) / nlevels(y1)

names <- c(tissue_meta$cell, pbmc_meta$Cell)

# Merge the two matrices
combined <- RowMergeSparseMatrices(tissue_norm, pbmc_norm)
combined <- combined[, names]  # Make sure that the order of cells are the same as the weights

# Find the variable genes
pbmc_genes <- readRDS(as.character(args[3]))@assays$RNA@var.features  # 3,000
tissue_genes <- readRDS(as.character(args[5]))$vargenes$symbol  # 5,484

genes <- unique(c(pbmc_genes, tissue_genes))  # 7,653

# Scale Data and PCs
combined_pcs <- singlecellmethods::weighted_pca(combined, weights=weights1, genes_use=genes)
saveRDS(combined_pcs, paste0("PCs/", name, "_PCs.Rds"))

# Save meta files for posterity
combined_meta <- data.frame(
    barcode = c(pbmc_meta$Cell, tissue_meta$cell),
    sample = c(paste0(pbmc_meta$Unified_Visit, "_PBMC"), paste0(tissue_meta$sample, "_Tissue")),
    annotation = c(pbmc_meta$annotation, tissue_meta$final_annotation),
    batch = c(paste0(pbmc_meta$Batch, "_PBMC"), paste0(tissue_meta$processing.batch, "_Tissue")),
    origin = c(rep("PBMC", length(pbmc_meta$Cell)), rep("Tissue", length(tissue_meta$cell))),
    unified_visit = c(pbmc_meta$Unified_Visit, tissue_meta$unified_visit)
) %>% arrange(match(barcode, names))

write.csv(combined_meta, paste0("PCs/", name, "_combinedMeta.csv"), row.names = FALSE)