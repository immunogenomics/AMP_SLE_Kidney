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
super_meta$percent.mito <- super_meta$percent.mito.sub

super_meta_nuc <- readRDS("super_meta_nuc.rds")

superUltraMeta <- bind_rows(
  super_meta %>% select(cells, nGene, nUMI, percent.mito) %>% mutate(type = "cell"),
  super_meta_nuc %>% mutate(type = "nuclei")
)

superUltraMeta_filtered <- bind_rows(
  super_meta %>% filter(nUMI > 1000, nGene > 500, percent.mito.sub < 0.03) %>% select(cells, nGene, nUMI, percent.mito) %>% mutate(type = "cell"),
  super_meta_nuc %>% filter(nUMI > 1000, nUMI < 40000,nGene > 500, nGene < 7500, percent.mito < 0.01) %>% select(cells, nGene, nUMI, percent.mito) %>% mutate(type = "nuclei")
)

superDuperUltraMeta <- bind_rows(
  superUltraMeta %>% mutate(status = "pre-QC"), superUltraMeta_filtered %>% mutate(status = "post-QC")
)

# Define the function
plot_violin_metric <- function(df, metric, y_label = NULL, title = NULL) {
  
  metric <- enquo(metric)  # tidy evaluation
  
  # Set default y-axis label
  if (is.null(y_label)) y_label <- quo_name(metric)
  
  # Factor for plotting order: cell pre-QC, cell post-QC, nuclei pre-QC, nuclei post-QC
  df <- df %>%
    mutate(
      type_status = factor(paste(type, ifelse(type %in% c("cell", "nuclei") & status != "pre-QC", "post-QC", "pre-QC")), 
                           levels = c("cell pre-QC", "cell post-QC", "nuclei pre-QC", "nuclei post-QC"))
    )
  
  # Colors
  fill_colors <- c("cell" = "#ffffb6", "nuclei" = "#bbb9d8")
  
  pdf(paste0("QC_metrics_violins_", quo_name(metric), ".pdf"), width = 8, height = 4)
    print(ggplot(df, aes(x = type_status, y = !!metric, fill = type)) +
      geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
      geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.5, outlier.shape = NA) +
      scale_fill_manual(values = fill_colors) +
      labs(x = "", y = y_label, title = title) +
      theme_minimal(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      ) + scale_y_log10())
  dev.off()
}

# Example usage:
# Violin plot of nGene
plot_violin_metric(superDuperUltraMeta, nGene, y_label = "Number of Genes", title = "nGene Distribution pre- vs post-QC")

# Violin plot of nUMI
plot_violin_metric(superDuperUltraMeta, nUMI, y_label = "Number of UMIs", title = "nUMI Distribution pre- vs post-QC")

# Violin plot of percent.mito
plot_violin_metric(superDuperUltraMeta, percent.mito, y_label = "Percent Mitochondrial Reads", title = "Percent Mito Distribution pre- vs post-QC")