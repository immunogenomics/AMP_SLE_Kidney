library(destiny)
library(tidyverse)
library(ggplot2)

# Test with hPCs
hPCs <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/TNK/harmony_0/TNK_combined_hPCs.Rds")
meta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/TNK_combinedMeta.csv")

# Take CD8 cells only
cd8_idx <- grepl("CD8", meta$annotation, ignore.case = TRUE) & 
    meta$unified_visit != "200-1138-V0" & 
    !(meta$annotationNumber %in% c("T4", "bl-T15", "bl-T18", "bl-T19"))

meta_8 <- meta[cd8_idx, ]
hPCs_8 <- hPCs[cd8_idx, ]  # subset by the SAME row indices

dir.create("CD8_inputs/")
saveRDS(hPCs_8, "CD8_inputs/CD8_T_cells_hPCs.Rds")
write.csv(meta_8, "CD8_inputs/CD8_T_cells_meta.csv")

# Take CD4 cells only
cd4_idx <- grepl("CD4", meta$annotation, ignore.case = TRUE) & 
    meta$unified_visit != "200-1138-V0" & 
    !(meta$annotationNumber %in% c("T4", "bl-T15", "bl-T18", "bl-T19"))

meta_4 <- meta[cd4_idx, ]
hPCs_4 <- hPCs[cd4_idx, ]

dir.create("CD4_inputs/")
saveRDS(hPCs_4, "CD4_inputs/CD4_T_cells_hPCs.Rds")
write.csv(meta_4, "CD4_inputs/CD4_T_cells_meta.csv")

dir.create("pt_outputs/")
dm8 <- DiffusionMap(hPCs_8)

saveRDS(dm8, paste0("pt_outputs/CD8_diffusionMap.Rds"))
saveRDS(dm8@eigenvectors, paste0("CD8_diffusionMap_eigenvectors.Rds"))

dm4 <- DiffusionMap(hPCs_4)
saveRDS(dm4, paste0("pt_outputs/CD4_diffusionMap.Rds"))
saveRDS(dm4@eigenvectors, paste0("pt_outputs/CD4_diffusionMap_eigenvectors.Rds"))