library(destiny)
library(tidyverse)
library(ggplot2)

# Test with hPCs
hPCs <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/TNK/harmony_0/TNK_combined_hPCs.Rds")
meta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/TNK_combinedMeta.csv")

# Removing proliferating T cells
tnk_meta <- meta[startsWith(meta$annotation, "T") | startsWith(meta$annotation, "bl-T"), ]
tnk_meta <- tnk_meta[!(tnk_meta$annotationNumber %in% c("T4", "bl-T15", "bl-T18", "bl-T19")), ]

hPCs <- as.matrix(hPCs[rownames(hPCs) %in% tnk_meta$barcode ,])

dm <- DiffusionMap(hPCs)

saveRDS(dm, paste0("T_cells_diffusionMap.Rds"))
saveRDS(dm@eigenvectors, paste0("T_cells_diffusionMap_eigenvectors.Rds"))

write.csv(tnk_meta, "noProlif_meta.csv")