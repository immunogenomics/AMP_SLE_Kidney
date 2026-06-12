suppressPackageStartupMessages({
    library(destiny)
    library(tidyverse)
    library(ggplot2)
})

# Test with normalized counts

# Test with hPCs
hPCs <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/B/harmony_0/B_combined_hPCs.Rds")
umap <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/B/harmony_0/B_umap.Rds")
meta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/B_combinedMeta.csv")

meta <- meta[startsWith(meta$annotation, "B") | startsWith(meta$annotation, "bl-B"), ]

hPCs <- as.matrix(hPCs[rownames(hPCs) %in% meta$barcode ,])

dm <- DiffusionMap(hPCs)
saveRDS(dm, paste0("B_cells_diffusionMap.Rds"))
saveRDS(dm@eigenvectors, "B_cells_diffusionMap_eigenvectors.Rds")