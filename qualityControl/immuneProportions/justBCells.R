library(tidyverse)
library(ggplot2)
library(scales)

load("../../SpatialRevisions/DanaherEtAl/cleaneddata.RData")
cannot <- read.csv("../../SpatialRevisions/DanaherEtAl/extraAnalyses/processed_data/cell positions wrt gloms.csv") # loads "cannot"

cellProportions <- data.frame(cells = names(clust), cluster = clust, tissue = annot$tissuename)

cellProportions$broad <- case_when(
  cellProportions$cluster %in% c("B-cell", "plasmablast") ~ "B Cells",
  TRUE ~ "nonImmune"
)

cellProportions$Type <- case_when(
  str_detect(cellProportions$tissue, "SLE") ~ "LN",
  TRUE ~ "Control"
)

# Bar chart of immune proportions by tissue
immuneProps <- cellProportions %>%
  group_by(tissue, broad, Type) %>%
  summarise(count = n()) %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count))

# Combine AMP dataset
myeloid <- readRDS("../../finalObjects_v3/myeloid_tissue_meta.Rds")  %>% select(cell, sample, dataset, Type) %>% mutate(broad = "immune")
b <- readRDS("../../finalObjects_v3/b_tissue_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "B Cells")
tnk <- readRDS("../../finalObjects_v3/t_nk_tissue_meta.Rds")  %>% select(cell, dataset, sample, Type) %>% mutate(broad = "immune")
glom <- readRDS("../../finalObjects_v3/glom_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonImmune")
loh <- readRDS("../../finalObjects_v3/loh_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonImmune")
pt <- readRDS("../../finalObjects_v3/pt_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonImmune")
intl <- readRDS("../../finalObjects_v3/intl_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonImmune")
dn <- readRDS("../../finalObjects_v3/dn_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonImmune")

amp <- rbind(myeloid, b, tnk, glom, loh, pt, intl, dn)

broad_levels <- amp %>%
  distinct(broad) %>%
  pull(broad)

ampProps <- amp %>%
  group_by(sample, broad, dataset, Type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample, dataset, Type) %>%
  complete(
    broad = broad_levels,
    fill = list(count = 0)
  ) %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count))

# Compare immune proportions for diseased smaples
allDanaher <- immuneProps %>% select(tissue, broad, Type, proportion) %>% mutate(dataset = "Spatial")
allAmp <- ampProps %>% mutate(tissue = sample) %>%
    ungroup() %>% 
    select(tissue, broad, Type, proportion, dataset)

allSamples <- rbind(allDanaher, allAmp)

allSamplesImmune <- allSamples %>% filter(broad == "B Cells")
allSamplesImmune$dataset <- factor(allSamplesImmune$dataset, levels = c("Spatial", "scRNAseq", "snRNAseq"))

pdf("BCellProp.pdf", width = 8, height = 6)
    print(ggplot(allSamplesImmune, aes(x = dataset, y = proportion, fill = dataset)) +
        geom_boxplot() +
        ggpubr::stat_compare_means(
          comparisons = list(c("Spatial", "scRNAseq"),
            c("Spatial", "snRNAseq")),
          label = "p.signif"
        ) +
        scale_fill_brewer(palette = "Set3", name = "Dataset") +
        scale_y_continuous(
          trans = 'pseudo_log',
        ) +
        facet_wrap(~ Type, scale = "free") +
        theme_minimal(base_size = 18) +
        labs(title = "B Cell Proportions by Dataset", x = "Dataset", y = "B Cell Proportion")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 18),       # facet labels
    axis.title = element_text(size = 18),       # axis labels
    legend.position = "none"
  )) 
dev.off()

pdf("BCellProp_noFacet.pdf", width = 8, height = 6)
    print(ggplot(allSamplesImmune, aes(x = dataset, y = proportion, fill = dataset)) +
        geom_boxplot() +
        # Danaher vs scRNAseq (exists everywhere)
        ggpubr::stat_compare_means(
          method = "t.test",
          comparisons = list(c("Danaher", "scRNAseq")),
          label = "p.signif"
        ) +
        # Danaher vs snRNAseq (ONLY where it exists)
        ggpubr::stat_compare_means(
          method = "t.test",
          data = subset(allSamplesImmune, Type != "Controls"),
          comparisons = list(c("Danaher", "snRNAseq")),
          label = "p.signif",
          label.y = 0.04
        ) +
        scale_fill_brewer(palette = "Set3", name = "Dataset") +
        scale_y_log10(breaks=c(0.001,.01,.1,1)) +        
        # facet_wrap(~ Type, scale = "fixed") +
        theme_minimal(base_size = 18) +
        labs(title = "B Cell Proportions by Dataset", x = "Dataset", y = "B Cell Proportion")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 18),       # facet labels
    axis.title = element_text(size = 18),       # axis labels
    legend.position = "none"
  )) 
dev.off()

## Immune Cells Only
load("../../SpatialRevisions/DanaherEtAl/cleaneddata.RData")
cannot <- read.csv("../../SpatialRevisions/DanaherEtAl/extraAnalyses/processed_data/cell positions wrt gloms.csv") # loads "cannot"

cellProportions <- data.frame(cells = names(clust), cluster = clust, tissue = annot$tissuename)

cellProportions$broad <- case_when(
  cellProportions$cluster %in% c("monocyte", "mast", "macrophage", "neutrophil", "mDC", "pDC") ~ "immune",
  cellProportions$cluster %in% c("B-cell", "plasmablast") ~ "B Cells",
  cellProportions$cluster %in% c("T CD4 naive", "T CD8 memory", "T CD8 naive", "T CD4 memory", "NK", "Treg") ~ "immune",
  TRUE ~ "nonImmune"
)

cellProportions$Type <- case_when(
  str_detect(cellProportions$tissue, "SLE") ~ "LN",
  TRUE ~ "Control"
)

cellProportions <- cellProportions[ cellProportions$broad != "nonImmune", ]

# Bar chart of immune proportions by tissue
immuneProps <- cellProportions %>%
  group_by(tissue, broad, Type) %>%
  summarise(count = n()) %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count))
  
amp <- rbind(myeloid, b, tnk)

ampProps <- amp %>%
  group_by(sample, broad, dataset, Type) %>%
  summarise(count = n()) %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count))

# Compare immune proportions for diseased smaples
allDanaher <- immuneProps %>% select(tissue, broad, Type, proportion) %>% mutate(dataset = "Danaher")
allAmp <- ampProps %>% mutate(tissue = sample) %>%
    ungroup() %>% 
    select(tissue, broad, Type, proportion, dataset)

allSamples <- rbind(allDanaher, allAmp)

allSamplesImmune <- allSamples %>% filter(broad == "B Cells")

pdf("BCellProp_immuneOnly.pdf", width = 8, height = 6)
    print(ggplot(allSamplesImmune, aes(x = dataset, y = proportion, fill = dataset)) +
        geom_boxplot() +
        # Danaher vs scRNAseq (exists everywhere)
        ggpubr::stat_compare_means(
          method = "t.test",
          comparisons = list(c("Danaher", "scRNAseq")),
          label = "p.signif",
        ) +
        # Danaher vs snRNAseq (ONLY where it exists)
        ggpubr::stat_compare_means(
          method = "t.test",
          data = subset(allSamplesImmune, Type != "Controls"),
          comparisons = list(c("Danaher", "snRNAseq")),
          label = "p.signif",
          label.y = 0.08
        ) +
        scale_fill_brewer(palette = "Set3", name = "Dataset") +
        scale_y_log10(breaks=c(0.001,.01,.1,1)) +        
        facet_wrap(~ Type, scale = "fixed") +
        theme_minimal(base_size = 18) +
        labs(title = "B Cell Proportions by Dataset", x = "Dataset", y = "B Cell Proportion")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 18),       # facet labels
    axis.title = element_text(size = 18),       # axis labels
    legend.position = "none"
  )) 
dev.off()

pdf("BCellProp_immuneOnly_noFacet.pdf", width = 8, height = 6)
    print(ggplot(allSamplesImmune, aes(x = dataset, y = proportion, fill = dataset)) +
        geom_boxplot() +
        scale_fill_brewer(palette = "Set3", name = "Dataset") +
        scale_y_log10(breaks=c(0,.01,.1,1)) +        
        theme_minimal(base_size = 18) +
        labs(title = "B Cell Proportions by Dataset", x = "Dataset", y = "B Cell Proportion")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 18),       # facet labels
    axis.title = element_text(size = 18),       # axis labels
    legend.position = "none"
  )) 
dev.off()