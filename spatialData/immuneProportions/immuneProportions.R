library(SeuratObject)
library(Matrix)
library(singlecellmethods)
library(dplyr)
library(tidyverse)
library(harmony)
library(pals)
library(ggplot2)
library(caret)
library(pROC)
library(ggrepel)
library(glmnet)
library(ggpubr)

set.seed(1)

load("../cleaneddata.RData")
cannot <- read.csv("../extraAnalyses/processed_data/cell positions wrt gloms.csv") # loads "cannot"

immuneCells <- c("macrophage", "monocyte", "Mesangial.cell", "mDC", "mast", "neturophil", "pDC",
                 "B-cell", "plasmablast", 
                 "Treg", "T CD8 memory", "T CD8 naive", "T CD4 memory", "T CD4 naive")

immuneSpatial <- names(clust)[clust %in% immuneCells]

cannotImmune <- cannot[cannot$id %in% immuneSpatial, ]

clustImmune <- clust[immuneSpatial]

cannotImmune$celltype <- clustImmune

cannotImmune$broad <- case_when(
  cannotImmune$celltype %in% c("macrophage", "monocyte", "Mesangial.cell", "mDC", "mast", "neturophil", "pDC") ~ "Myeloid",
  cannotImmune$celltype %in% c("B-cell", "plasmablast") ~ "Bcell",
  cannotImmune$celltype %in% c("Treg", "T CD8 memory", "T CD8 naive", "T CD4 memory", "T CD4 naive") ~ "TNK"
)

caseSamples <- annot %>% 
  filter(grepl("SLE", tissuename)) %>% 
  pull(cell_ID)

cannotImmune <- filter(cannotImmune, id %in% caseSamples)

# Make sure position.vs.glom is a factor
cannotImmune <- cannotImmune %>%
  mutate(position.vs.glom = factor(position.vs.glom, levels = c("inside glomerulous", "bordering glomerulous", "tubulointerstitium")))

# Plot stacked proportion bar chart
pdf("broad_position_stacked_proportion.pdf", width = 8, height = 6)
ggplot(cannotImmune, aes(x = broad, fill = position.vs.glom)) +
  geom_bar(position = "fill") +   # <-- normalize bars to 1
  scale_y_continuous(labels = scales::percent_format()) +  # show y-axis in %
  labs(
    x = "Broad Cell Category",
    y = "Proportion of Cells",
    fill = "Position vs Glomerulus",
    title = "Cell Position Distribution by Broad Category"
  ) +
  theme_minimal(base_size = 14) + 
  scale_fill_manual(values = c(
    "tubulointerstitium" = "steelblue", 
    "bordering glomerulous" = "gold", 
    "inside glomerulous" = "red"
  ))
dev.off()

# Compute counts for labels
count_df <- cannotImmune %>%
  group_by(broad, position.vs.glom) %>%
  summarise(n = n(), .groups = "drop")

pdf("broad_position_stacked_counts_facet.pdf", width = 10, height = 6)
ggplot(count_df, aes(x = "", y = n, fill = position.vs.glom)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 4) +
  facet_wrap(~broad, scales = "free_y") +
  labs(
    x = "",
    y = "Number of Cells",
    fill = "Position vs Glomerulus",
    title = "Cell Position Counts by Broad Category"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c(
    "tubulointerstitium" = "steelblue", 
    "bordering glomerulous" = "gold", 
    "inside glomerulous" = "red"
  ))
dev.off()

# --------------------------
# 3. Overall stacked bar ignoring broad
# --------------------------
overall_df <- cannotImmune %>%
  group_by(position.vs.glom) %>%
  summarise(n = n(), .groups = "drop")

pdf("overall_position_stacked_counts.pdf", width = 6, height = 6)
ggplot(overall_df, aes(x = "", y = n, fill = position.vs.glom)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 4) +
  labs(
    x = "",
    y = "Number of Cells",
    fill = "Position vs Glomerulus",
    title = "Overall Cell Position Distribution"
  ) +
  theme_minimal(base_size = 14) + 
  scale_fill_manual(values = c(
    "tubulointerstitium" = "steelblue", 
    "bordering glomerulous" = "gold", 
    "inside glomerulous" = "red"
  ))
dev.off()