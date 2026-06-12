library(tidyverse)

path <- "/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/raw_metadata/"

# Get all files ending with "_raw_meta.rds"
files <- list.files(
  path = path,
  pattern = "_raw_meta\\.rds$",
  full.names = TRUE
)

# Read and combine
meta_list <- lapply(files, function(f) {
  df <- readRDS(f)
  df$source_file <- basename(f)
  df
})

super_meta <- bind_rows(meta_list) %>% filter(!(source_file %in% c("AMPSLEkid_cells_1258_raw_meta.rds", "AMPSLEkid_cells_0135_raw_meta.rds", "AMPSLEkid_cells_1138_raw_meta.rds")))

# Prepare a long‑format dataframe
hist_df <- super_meta %>%
  select(nGene, nUMI, percent.mito.sub) %>%
  pivot_longer(cols = everything(),
               names_to = "metric",
               values_to = "value") %>%
  mutate(
    plot_value = ifelse(metric %in% c("nGene", "nUMI"),
                        log10(value),
                        value)
  )

# Define cutoff values
cutoffs <- data.frame(
  metric = c("nGene", "nUMI", "percent.mito.sub"),
  cutoff = c(500, 1000, 0.03)
) %>%
  mutate(
    plot_cutoff = ifelse(metric %in% c("nGene", "nUMI"),
                         log10(cutoff),
                         cutoff)
)


pdf("QC_metrics_violins.pdf", width = 12, height = 6)

  ggplot(hist_df, aes(x = "", y = plot_value)) +
    geom_violin(fill = "#ffffb6", color = "black", trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "#ffffb6", outlier.shape = NA) +
    geom_hline(data = cutoffs,
              aes(yintercept = plot_cutoff),
              color = "red",
              linetype = "dashed",
              linewidth = 1) +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    theme_classic(base_size = 15) +
    labs(x = NULL, y = "Value (log10 for nGene/nUMI)",
        title = "QC violin plots with cutoffs") +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

dev.off()

# Single-nuc
path <- "/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/raw_metadata/"

# Get all files ending with "_raw_meta.rds"
files <- list.files(
  path = path,
  pattern = "_raw_meta\\.rds$",
  full.names = TRUE
)

# Read and combine
meta_list <- lapply(files, function(f) {
  df <- readRDS(f)
  df$source_file <- basename(f)
  df
})

samples <- read.delim("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/list_of_nuc_samples.txt",
                            sep = "\n", header = F)
samples <- as.character(samples[,1])

# Read and combine
meta_list <- lapply(paste0("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/raw_metadata/", samples, "_raw_meta.rds"), function(f) {
  print(f)
  df <- readRDS(f)
  df <- df %>% filter(nUMI > 10)
  str(df)
  df$source_file <- basename(f)
  df
})

gc()
super_meta_nuc <- bind_rows(meta_list) %>% filter(!(source_file %in% c("AMPSLEkid_nuc_1258_raw_meta.rds", "AMPSLEkid_nuc_0135_raw_meta.rds", "AMPSLEkid_nuc_1138_raw_meta.rds")))

saveRDS(super_meta_nuc, "super_meta_nuc.rds")

# Prepare a long-format dataframe
hist_df <- super_meta_nuc %>%
  select(nGene, nUMI, percent.mito) %>%
  pivot_longer(cols = everything(),
               names_to = "metric",
               values_to = "value")

# Define multiple cutoff values per metric


pdf("QC_metrics_histograms_sn.pdf", width = 10, height = 8)
    # Plot histograms with vertical cutoff lines
    ggplot(hist_df, aes(x = value)) +
    geom_histogram(bins = 100, fill = "#bbb9d8", color = "black") +
    geom_vline(data = cutoffs, aes(xintercept = cutoff), 
                color = "red", linetype = "dashed", size = 1) +
    facet_wrap(~metric, scales = "free", ncol = 1) +
    theme_classic(base_size = 15) +
    labs(x = "Value", y = "Number of Cells",
        title = "QC Histograms with Cutoffs")
dev.off()

# Prepare a long‑format dataframe
hist_df <- super_meta_nuc %>%
  select(nGene, nUMI, percent.mito) %>%
  pivot_longer(cols = everything(),
               names_to = "metric",
               values_to = "value") %>%
  mutate(
    plot_value = ifelse(metric %in% c("nGene", "nUMI"),
                        log10(value),
                        value)
  )

# Define cutoff values
cutoffs <- data.frame(
  metric = c("nGene", "nGene", "nUMI", "nUMI", "percent.mito"),
  cutoff = c(500, 7500, 1000, 40000, 0.01)
) %>%
  mutate(
    plot_cutoff = ifelse(metric %in% c("nGene", "nUMI"),
                         log10(cutoff),
                         cutoff)
)


pdf("QC_metrics_violins_sn.pdf", width = 12, height = 6)

  ggplot(hist_df, aes(x = "", y = plot_value)) +
    geom_violin(fill = "#bbb9d8", color = "black", trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "#bbb9d8", outlier.shape = NA) +
    geom_hline(data = cutoffs,
              aes(yintercept = plot_cutoff),
              color = "red",
              linetype = "dashed",
              linewidth = 1) +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    theme_classic(base_size = 15) +
    labs(x = NULL, y = "Value (log10 for nGene/nUMI)",
        title = "QC violin plots with cutoffs") +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

dev.off()
