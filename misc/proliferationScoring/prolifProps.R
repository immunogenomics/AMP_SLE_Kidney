suppressPackageStartupMessages({
  library(Matrix)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(singlecellmethods)
})

meta <- readRDS("../finalObjects_v2/myeloid_tissue_meta.Rds") %>% select(cell, sample, dataset, final_annotation)

# Define the annotation of interest
target_annotation <- "M15. CENPF+ MKI67+ Proliferating"

# Calculate the percentage per sample
proportions <- meta %>%
  group_by(sample) %>%
  summarise(
    total = n(),
    count = sum(final_annotation == target_annotation),
    percentage = count / total
)

clinData <- readRDS("../finalObjects/clinicalData.rds") %>% select(sample, Final_Chronicity, Final_Activity)
clinData$Final_Chronicity <- as.vector(as.character(clinData$Final_Chronicity))
clinData_clean <- data.frame(sample = clinData$sample, Chronicity = as.character(clinData$Final_Chronicity), Activity = clinData$Final_Activity)

# Take only sc data
proportions <- merge(proportions, clinData_clean, by = "sample")
proportions <- proportions %>% filter(!is.na(Chronicity) & !is.na(Activity))  # Filter out rows with NA values; 141 rows remain
proportions$Chronicity <- as.integer(proportions$Chronicity)

# Plot with Chronicity
p <- ggplot(data = proportions, aes(x = Chronicity, y = percentage, group = Chronicity)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Proportion of Proliferating Macrophages by Chronicity", x = "Chronicity", y = "Proliferating Proportion") +
  theme(legend.position = "none") + scale_x_continuous(breaks = seq(0, 10, by = 1))

ggsave("proliferation_proportion_by_chronicity.pdf", plot = p, width = 8, height = 6)

# Plot with Activity
p <- ggplot(data = proportions, aes(x = Activity, y = percentage, group = Activity)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Proportion of Proliferating Macrophages by Activity", x = "Activity", y = "Proliferating Proportion") +
  theme(legend.position = "none") + scale_x_continuous(breaks = seq(0, max(proportions$Activity), by = 1))

ggsave("proliferation_proportion_by_activity.pdf", plot = p, width = 8, height = 6)