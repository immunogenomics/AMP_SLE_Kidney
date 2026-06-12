library(destiny)
library(tidyverse)
library(ggplot2)

# Test with normalized counts

# Test with hPCs
hPCs <- readRDS("corrected_Tcells_hPCs.Rds")
meta <- readRDS("corrected_Tcells_meta.Rds")

tnk_meta <- meta[startsWith(meta$annotation, "T") | startsWith(meta$annotation, "bl-T"), ]

hPCs <- as.matrix(hPCs[rownames(hPCs) %in% tnk_meta$normalized_channels ,])

dm <- DiffusionMap(hPCs)
saveRDS(dm, paste0("T_cells_diffusionMap.Rds"))
saveRDS(dm@eigenvectors, paste0("T_cells_diffusionMap_eigenvectors.Rds"))