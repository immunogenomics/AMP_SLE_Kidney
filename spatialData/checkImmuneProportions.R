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
library(ggnewscale)

set.seed(1)

load("cleaneddata.RData")
cannot <- read.csv("extraAnalyses/processed_data/cell positions wrt gloms.csv") # loads "cannot"

cellProportions <- data.frame(cells = names(clust), cluster = clust, tissue = annot$tissuename)

cellProportions$broad <- case_when(
  cellProportions$cluster %in% c("monocyte", "mast", "macrophage", "neutrophil", "mDC", "pDC") ~ "Myeloid",
  cellProportions$cluster %in% c("B-cell", "plasmablast") ~ "Bcell",
  cellProportions$cluster %in% c("T CD4 naive", "T CD8 memory", "T CD8 naive", "T CD4 memory", "NK", "Treg") ~ "TNK",
  TRUE ~ "Other"
)

# Bar chart of immune proportions by tissue
immuneProps <- cellProportions %>% filter(broad %in% c("Myeloid", "Bcell", "TNK")) %>%
  group_by(tissue, broad) %>%
  summarise(count = n()) %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count))

pdf("immune_proportions_by_tissue.pdf", width = 6, height = 4)
    ggplot(immuneProps, aes(x = tissue, y = proportion, fill = broad)) +
        geom_bar(stat = "identity", position = "fill") +
        theme_minimal() +
        labs(title = "Immune Cell Proportions by Tissue", x = "Tissue", y = "Proportion") +
        scale_fill_brewer(palette = "Set2", name = "Immune Cell Type") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off() 


# Bar chart of immune proportions by tissue
cellProps <- cellProportions %>%
  group_by(tissue, broad) %>%
  summarise(count = n()) %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count))

cellProps <- cellProps[cellProps$broad != "Other", ]

pdf("immune_proportions_by_tissue_total.pdf", width = 6, height = 4)
    ggplot(cellProps, aes(x = tissue, y = proportion, fill = broad)) +
        geom_bar(stat = "identity", position = "stack") +
        theme_minimal() +
        labs(title = "Immune Cell Proportions by Tissue", x = "Tissue", y = "Proportion") +
        scale_fill_brewer(palette = "Set2", name = "Immune Cell Type") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off() 

saveRDS(cellProps, "Danaher_immuneProps.Rds")
