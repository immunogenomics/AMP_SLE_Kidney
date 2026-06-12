library(tidyverse)
library(ggplot2)
library(Seurat)
library(Matrix)

# Load hPCs and metadata
meta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/TNK_combinedMeta.csv")
hPCs <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/TNK/harmony_0/TNK_combined_hPCs.Rds")

meta <- meta %>%
  mutate(
    normalized_channels = if_else(
      origin == "PBMC",
      normalize_channel(barcode),
      barcode
    )
  )

usages <- read.csv("tcat_usage.csv")

cellCycle <- usages[, colnames(usages) %in% c("gene_symbol", "CellCycle.G2M", "CellCycle.S", "CellCycle.Late.S")]

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


rownames(hPCs) <- meta$normalized_channels

meta <- meta %>% filter(normalized_channels %in% cellCycle$gene_symbol)
meta <- meta[match(cellCycle$gene_symbol, meta$normalized_channels), ]
hPCs <- hPCs[rownames(hPCs) %in% meta$normalized_channels, ]
hPCs <- hPCs[match(meta$normalized_channels, rownames(hPCs)), ]

X <- cbind(Intercept = 1, as.matrix(cellCycle[, c("CellCycle.G2M", "CellCycle.S", "CellCycle.Late.S")]))

# Fit all PCs at once
fit <- lm.fit(x = X, y = as.matrix(hPCs))

# Get the residuals
hPCs_corrected <- residuals(fit)

saveRDS(hPCs_corrected, "corrected_Tcells_hPCs.Rds")
saveRDS(meta, "corrected_Tcells_meta.Rds")
