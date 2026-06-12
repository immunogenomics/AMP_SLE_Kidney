library(RANN)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pals)
library(ggthemes)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)
name <- as.character(args[1])
theta <- as.numeric(args[2])

numNeigh = 100

# Read in PCs
meta <- read.csv(paste0("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/", name, "_combinedMeta.csv"), header = TRUE)

# Find nearest neighbors of all
neighbors <- readRDS(paste0(name, "_100.Rds"))

tissueOfOrigin <- apply(neighbors$nn.idx, 1:2, function(x, row) meta$origin[x])
tissueOfOrigin <- cbind(meta$annotation, tissueOfOrigin)

weightFactor <- (length(meta$origin[meta$origin == "Tissue"]) / 
        length(meta$origin[meta$origin == "PBMC"]))

results <- data.frame(matrix(nrow = nrow(tissueOfOrigin), ncol = 0))

for (i in seq(20, 100, by = 10)) {
    numNeigh = i
    TSM <- apply(tissueOfOrigin[, 1:(numNeigh + 1)], 1, function(x){
        sum(x[2:length(x)] == "Tissue") / (sum(x[2:length(x)] == "Tissue") + weightFactor * sum(x[2:length(x)] == "PBMC"))
    })
    results <- cbind(results, TSM)
}

colnames(results) <- paste0("k_", seq(20, 100, by = 10))
results$annotation <- tissueOfOrigin[, 1]

# Group annotations by similar names and calculate median TSM with error bars
df_long <- results %>%
  pivot_longer(
    cols = starts_with("k_"),
    names_to = "k",
    values_to = "value"
  )

df_long <- df_long %>%
  mutate(
    k = str_remove(k, "^k_"),
    k = as.numeric(k)
  )

df_long <- df_long %>%
  mutate(
    annotation = str_remove(annotation, "\\..*$")
  )

df_summary <- df_long %>%
  group_by(annotation, k) %>%
  summarise(
    median = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

df_summary <- df_summary %>%
  mutate(
    is_bl = str_starts(annotation, regex("bl-", ignore_case = TRUE)),
    num   = as.numeric(str_extract(annotation, "\\d+$"))
  ) %>%
  arrange(is_bl, num) %>%
  mutate(
    annotation = factor(annotation, levels = unique(annotation))
  ) %>%
  select(-is_bl, -num)

palette <- as.vector(c(
    as.vector(rev(polychrome(26)))[1:length(unique(meta$annotation[meta$origin == "Tissue"]))],
    as.vector(ggthemes_data$tableau$'color-palettes'$regular$`Classic 20`[1])[1:length(unique(meta$annotationp[meta$origin == "PBMC"]))]$value
))

palette_named <- setNames(
  palette[seq_along(levels(df_summary$annotation))],
  levels(df_summary$annotation)
)

p <- ggplot(df_summary, aes(x = k, y = median, color = annotation)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "black") +
    theme_bw(base_size = 20) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_blank()
    ) +
    labs(
        x = "k",
        y = "Median TSM",
        color = "Cell State"
    ) + scale_color_manual(values = palette) 

# Extract legend
legend <- get_legend(p)

# Save legend as its own image
ggsave(
  filename = paste0(name, "_robustPlot_legend.pdf"),
  plot = legend,
  width = 4,
  height = 6
)

p_nolegend <- p + theme(legend.position = "none")

pdf(paste0(name, "_robustPlot.pdf"), width = 8, height = 6)
  print(p_nolegend)
dev.off()