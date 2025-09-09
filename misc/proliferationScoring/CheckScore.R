suppressPackageStartupMessages({
  library(Matrix)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(singlecellmethods)
})

meta <- readRDS("../finalObjects_v2/myeloid_tissue_meta.Rds") %>% select(cell, sample, dataset, final_annotation) %>% filter(str_detect(final_annotation, regex("Macrophage|Proliferating")))
norm <- readRDS("../finalObjects_v2/myeloid_tissue_norm.Rds")

clinData <- readRDS("../finalObjects/clinicalData.rds") %>% select(sample, Final_Chronicity, Final_Activity)
clinData$Final_Chronicity <- as.vector(as.character(clinData$Final_Chronicity))

clinData_clean <- data.frame(sample = clinData$sample, Chronicity = as.character(clinData$Final_Chronicity), Activity = clinData$Final_Activity)

meta <- merge(meta, clinData_clean, by = "sample", all.x = TRUE)
meta <- meta %>% filter(!is.na(Chronicity))
meta$Chronicity <- as.numeric(meta$Chronicity)

# Cell Cycle Genes
cc_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
cc_genes <- cc_genes[cc_genes != "E2F8"]
norm <- norm[rownames(norm) %in% cc_genes, colnames(norm) %in% meta$cell]

# Calculate Proliferation Scores
scaled <- scaleData(norm)
proliferation <- colMeans(scaled)
proliferation <- proliferation[names(proliferation) %in% meta$cell]

profScores <- as.data.frame(proliferation) %>%
  rownames_to_column("cell") %>%
  left_join(meta, by = "cell") %>%
  select(cell, sample, Chronicity, Activity, final_annotation, proliferation)

profScores$Chronicity <- factor(as.character(profScores$Chronicity), levels = as.character(0:10))
# Plot with Chronicity
p <- ggplot(data = profScores, aes(x = final_annotation, y = proliferation, group = final_annotation, fill = final_annotation)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  labs(title = "Proliferation Scores by Cell Type", x = "Cell Type", y = "Proliferation Score")

ggsave("profScores_by_cellType_byCellType.pdf", plot = p, width = 8, height = 6)

# Plot with Activity
p <- ggplot(data = profScores, aes(x = Activity, y = proliferation, group = Activity)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Proliferation Scores by Activity", x = "Activity", y = "Proliferation Score") +
  theme(legend.position = "none") + scale_x_continuous(breaks = seq(0, max(profScores$Activity), by = 1))

ggsave("SAMs_proliferation_scores_by_activity.pdf", plot = p, width = 8, height = 6)
