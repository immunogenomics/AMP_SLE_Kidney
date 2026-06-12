# Read in single cell data
suppressPackageStartupMessages({
    library(tidyverse)
    library(Matrix)
    library(ggplot2)
    library(scales)
    library(ggrepel)
})

path <- "/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/raw_metadata/"

# Get all files ending with "_raw_meta.rds"
files <- list.files(
  path = path,
  pattern = "_raw_meta\\.rds$",
  full.names = TRUE
)

# Read and combine
meta_list_cells <- lapply(files, function(f) {
  df <- readRDS(f)
  df$source_file <- basename(f)
  df
})

super_meta_cells <- bind_rows(meta_list_cells) %>% filter(!(source_file %in% c("AMPSLEkid_cells_1258_raw_meta.rds", "AMPSLEkid_cells_0135_raw_meta.rds", "AMPSLEkid_cells_1138_raw_meta.rds")))

# Define a sequence of thresholds
thresholds <- seq(0, max(super_meta_cells$percent.mito.sub, na.rm = TRUE), by = 0.01)

goodCells <- super_meta_cells %>% filter(nGene > 500, nUMI > 1000)

senstivity <- data.frame(
  twoPercent = goodCells %>% filter(percent.mito.sub <= 0.02) %>% nrow() / nrow(goodCells),
  twoFivePercent = goodCells %>% filter(percent.mito.sub <= 0.025) %>% nrow() / nrow(goodCells),
  threePercent = goodCells %>% filter(percent.mito.sub <= 0.03) %>% nrow() / nrow(goodCells),
  threeFivePercent = goodCells %>% filter(percent.mito.sub <= 0.035) %>% nrow() / nrow(goodCells),
  fourPercent = goodCells %>% filter(percent.mito.sub <= 0.04) %>% nrow() / nrow(goodCells)
)

# Convert to long format for ggplot
sensitivity_long <- senstivity %>%
  pivot_longer(cols = everything(), names_to = "Threshold", values_to = "Proportion")

sensitivity_long$Threshold <- factor(sensitivity_long$Threshold, levels = c("twoPercent", "twoFivePercent", "threePercent", "threeFivePercent", "fourPercent"),
                                      labels = c("2%", "2.5%", "3%", "3.5%", "4%"))

pdf("mito_sensitivity.pdf", width = 6, height = 4)
  # Plot
ggplot(sensitivity_long, aes(x = Threshold, y = Proportion)) +
  geom_col(fill = "steelblue", width = 0.6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Mitochondrial Percent Threshold",
    y = "Proportion of Cells Retained",
    title = "Sensitivity of Cell Retention to Mitochondrial Cutoff"
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")
dev.off()
