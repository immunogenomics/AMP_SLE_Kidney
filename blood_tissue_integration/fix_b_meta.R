suppressPackageStartupMessages({
    library(tidyverse)
})

# Fix B Cells
# Load the data
meta <- read.csv("PCs/B_combinedMeta.csv")
newMeta <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/reorderedClusters/b_tissue_meta_reOrdered.Rds")

# Replace tissue clusters
meta <- meta %>%
  left_join(newMeta %>% select(cell, final_annotation), by = c("barcode" = "cell")) %>%
  mutate(annotation = if_else(!is.na(final_annotation), final_annotation, annotation)) %>%
  select(-final_annotation)

# Replace PBMC clusters
pbmcMeta <- read.csv("../newClusterConversion/B_PBMC.csv")

meta <- meta %>%
  left_join(pbmcMeta, by = c("annotation" = "oldAnnotation")) %>%
  mutate(annotation = if_else(!is.na(newAnnotation), newAnnotation, annotation)) %>%
  select(-newAnnotation)

meta$annotationNumber <- sub("\\..*", "", meta$annotation)
write.csv(meta, "PCs/B_combinedMeta.csv", row.names = FALSE)

# Fix T Cells
meta <- read.csv("PCs/TNK_combinedMeta.csv")
tissueMeta <- read.csv("../newClusterConversion/TNK_Tissue.csv")
pbmcMeta <- read.csv("../newClusterConversion/TNK_PBMC.csv")

annotation_map <- rbind(tissueMeta, pbmcMeta)

# Join and update annotation values
meta <- meta %>%
  left_join(annotation_map, by = c("annotation" = "oldAnnotation")) %>%
  mutate(annotation = if_else(!is.na(newAnnotation), newAnnotation, annotation)) %>%
  select(-newAnnotation)

meta$annotationNumber <- sub("\\..*", "", meta$annotation)
write.csv(meta, "PCs/TNK_combinedMeta.csv", row.names = FALSE)

# Fix Myeloid Cells
meta <- read.csv("PCs/Myeloid_combinedMeta.csv")
tissueMeta <- read.csv("../newClusterConversion/Myeloid_Tissue.csv")
pbmcMeta <- read.csv("../newClusterConversion/Myeloid_PBMC.csv")

annotation_map <- rbind(tissueMeta, pbmcMeta)

# Join and update annotation values
meta <- meta %>%
  left_join(annotation_map, by = c("annotation" = "oldAnnotation")) %>%
  mutate(annotation = if_else(!is.na(newAnnotation), newAnnotation, annotation)) %>%
  select(-newAnnotation)

meta$annotationNumber <- sub("\\..*", "", meta$annotation)
write.csv(meta, "PCs/Myeloid_combinedMeta.csv", row.names = FALSE)