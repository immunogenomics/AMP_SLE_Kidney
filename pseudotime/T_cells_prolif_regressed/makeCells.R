library(tidyverse)
library(ggplot2)
library(Seurat)
library(Matrix)

# Test with normalized counts
tissue <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/data/qcd/T_NK_clusterQCd_rawcounts_03062023.rds")
blood <- readRDS("/data/srlab1/Yu/tissue_blood_eQTL/AMPII_SLE/AMPII_blood/data/pheno/AllSample/GEX_raw.rds")

# Test with hPCs
meta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/TNK_combinedMeta.csv")

normalize_channel <- function(x) {
  x %>%
    gsub("channel_", "channel", .) %>%
    gsub("_", "-", .) %>%
    gsub("channel([0-9]+)-", "channel\\1-", .)
}

meta <- meta %>%
  mutate(
    normalized_channels = if_else(
      origin == "PBMC",
      normalize_channel(barcode),
      barcode
    )
  )

tissue <- tissue[, colnames(tissue) %in% meta$normalized_channels]
blood <- blood[, colnames(blood) %in% meta$normalized_channels]

joint <- RowMergeSparseMatrices(tissue, blood)
joint <- joint[, match(meta$normalized_channels, colnames(joint))]

writeMM(joint, file = "jointT.mtx")

write.table(rownames(joint),
            file = "jointT_rows.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(colnames(joint),
            file = "jointT_cols.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)