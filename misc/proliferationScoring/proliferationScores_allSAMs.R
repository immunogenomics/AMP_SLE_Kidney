suppressPackageStartupMessages({
  library(Matrix)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(singlecellmethods)
})

SAMs <- c("M5. GPNMBhigh NUPR1high Macrophage", "M7. SPP1high FABP5high Macrophage", "M9. MERTKhigh FABP5high Macrophage", "M11. GPMNBhigh NUPR1low Macrophage")
meta <- readRDS("../finalObjects_v2/myeloid_tissue_meta.Rds") %>% select(cell, sample, dataset, final_annotation) %>% filter(final_annotation %in% SAMs)
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
  select(cell, sample, Chronicity, Activity, proliferation)

# Plot with Chronicity
p <- ggplot(data = profScores, aes(x = Chronicity, y = proliferation, group = Chronicity)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Proliferation Scores by Chronicity", x = "Chronicity", y = "Proliferation Score") +
  theme(legend.position = "none") + scale_x_continuous(breaks = seq(0, 10, by = 1))

ggsave("SAMs_proliferation_scores_by_chronicity.pdf", plot = p, width = 8, height = 6)

# Plot with Activity
p <- ggplot(data = profScores, aes(x = Activity, y = proliferation, group = Activity)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Proliferation Scores by Activity", x = "Activity", y = "Proliferation Score") +
  theme(legend.position = "none") + scale_x_continuous(breaks = seq(0, max(profScores$Activity), by = 1))

ggsave("SAMs_proliferation_scores_by_activity.pdf", plot = p, width = 8, height = 6)
