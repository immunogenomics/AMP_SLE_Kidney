library(tidyverse)
library(ggplot2)

load("../../SpatialRevisions/DanaherEtAl/cleaneddata.RData")
cannot <- read.csv("../../SpatialRevisions/DanaherEtAl/extraAnalyses/processed_data/cell positions wrt gloms.csv") # loads "cannot"

cellProportions <- data.frame(cells = names(clust), cluster = clust, tissue = annot$tissuename)

cellProportions$broad <- case_when(
  cellProportions$cluster %in% c("monocyte", "mast", "macrophage", "neutrophil", "mDC", "pDC") ~ "immune",
  cellProportions$cluster %in% c("B-cell", "plasmablast") ~ "immune",
  cellProportions$cluster %in% c("T CD4 naive", "T CD8 memory", "T CD8 naive", "T CD4 memory", "NK", "Treg") ~ "immune",
  TRUE ~ "nonImmune"
)

cellProportions$cellType <- case_when(
  cellProportions$cluster %in% c("monocyte", "mast", "macrophage", "neutrophil", "mDC", "pDC") ~ "Myeloid",
  cellProportions$cluster %in% c("B-cell", "plasmablast") ~ "B Cells",
  cellProportions$cluster %in% c("T CD4 naive", "T CD8 memory", "T CD8 naive", "T CD4 memory", "NK", "Treg") ~ "T and NK",
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
myeloid <- readRDS("../../finalObjects_v3/myeloid_tissue_meta.Rds")  %>% select(cell, sample, dataset, Type) %>% mutate(broad = "immune", cellType = "Myeloid")
b <- readRDS("../../finalObjects_v3/b_tissue_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "immune", cellType = "B Cells")
tnk <- readRDS("../../finalObjects_v3/t_nk_tissue_meta.Rds")  %>% select(cell, dataset, sample, Type) %>% mutate(broad = "immune", cellType = "T and NK")
glom <- readRDS("../../finalObjects_v3/glom_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonImmune", cellType = "nonImmune")
loh <- readRDS("../../finalObjects_v3/loh_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonImmune", cellType = "nonImmune")
pt <- readRDS("../../finalObjects_v3/pt_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonImmune", cellType = "nonImmune")
intl <- readRDS("../../finalObjects_v3/intl_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonImmune", cellType = "nonImmune")
dn <- readRDS("../../finalObjects_v3/dn_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonImmune", cellType = "nonImmune")

amp <- rbind(myeloid, b, tnk, glom, loh, pt, intl, dn)

ampProps <- amp %>%
  group_by(sample, broad, dataset, Type) %>%
  summarise(count = n()) %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count))

# Compare immune proportions for diseased smaples
allDanaher <- immuneProps %>% select(tissue, broad, Type, proportion) %>% mutate(dataset = "Spatial")
allAmp <- ampProps %>% mutate(tissue = sample) %>%
    ungroup() %>% 
    select(tissue, broad, Type, proportion, dataset)

allSamples <- rbind(allDanaher, allAmp)

allSamplesImmune <- allSamples %>% filter(broad == "immune")
allSamplesImmune$dataset <- factor(allSamplesImmune$dataset, levels = c("Spatial", "scRNAseq", "snRNAseq"))

# pdf("immuneVsNonImmune.pdf", width = 8, height = 6)
#     print(ggplot(allSamplesImmune, aes(x = dataset, y = proportion, fill = dataset)) +
#         geom_boxplot() +
#         ggpubr::stat_compare_means(
#           comparisons = list(c("Spatial", "scRNAseq"), c("Spatial", "snRNAseq")),
#           label = "p.signif"
#         ) +
#         scale_fill_brewer(palette = "Set3", name = "Dataset") +
#         scale_y_log10(breaks=c(0.001,.01,.1,1)) +        
#         facet_wrap(~ Type, scale = "fixed") +
#         theme_minimal(base_size = 18) +
#         labs(title = "Immune Cell Proportions by Dataset", x = "Dataset", y = "Proportion")+
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.text = element_text(size = 18),       # facet labels
#     axis.title = element_text(size = 18),       # axis labels
#     legend.position = "none"
#   )) 
# dev.off()

ampCounts <- rbind(immuneProps %>% select(tissue, broad, Type, count) %>% mutate(dataset = "Spatial"), 
  allAmp <- ampProps %>% mutate(tissue = sample) %>%
      ungroup() %>% 
      select(tissue, broad, Type, count, dataset)) %>% filter(broad == "immune")

ampCounts$dataset <- factor(ampCounts$dataset, levels = c("Spatial", "scRNAseq", "snRNAseq"))

# pdf("medianCells.pdf", width = 8, height = 6)
#     print(ggplot(ampCounts, aes(x = dataset, y = count, fill = dataset)) +
#         geom_boxplot() +
#         geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +
#         # ggpubr::stat_compare_means(
#         #   comparisons = list(c("Spatial", "scRNAseq"), c("Spatial", "snRNAseq")),
#         #   label = "p.signif"
#         # ) +
#         scale_fill_brewer(palette = "Set3", name = "Dataset") +
#         scale_y_log10() +        
#         facet_wrap(~ Type, scale = "fixed") +
#         theme_minimal(base_size = 18) +
#         labs(title = "Immune Cell Counts per Biopsy", x = "Dataset", y = "Cell Counts")+
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.text = element_text(size = 18),       # facet labels
#     axis.title = element_text(size = 18),       # axis labels
#     legend.position = "none"
#   )) 
# dev.off()


# Same plot, but stratified by cell type
immuneProps <- cellProportions %>%
  group_by(tissue, cellType, Type) %>%
  summarise(count = n()) %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count))

broad_levels <- amp %>%
  distinct(cellType) %>%
  pull(cellType)

ampProps <- amp %>%
  group_by(sample, cellType, dataset, Type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample, dataset, Type) %>%
  complete(
    cellType = broad_levels,
    fill = list(count = 0)
  ) %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count))

allDanaher <- immuneProps %>% select(tissue, cellType, Type, proportion) %>% mutate(dataset = "Spatial")
allAmp <- ampProps %>% mutate(tissue = sample) %>%
    ungroup() %>% 
    select(tissue, cellType, Type, proportion, dataset)

allSamples <- rbind(allDanaher, allAmp)

allSamplesImmune <- allSamples %>% filter(cellType != "nonImmune") %>%
mutate(proportion = if_else(proportion == 0,
                            1e-04,
                            proportion))

allSamplesImmune$dataset <- factor(allSamplesImmune$dataset, levels = c("Spatial", "scRNAseq", "snRNAseq")) 

pdf("immuneVsNonImmune_faceted.pdf", width = 12, height = 4)

print(
  ggplot(allSamplesImmune, aes(x = dataset, y = proportion, fill = dataset)) +
    geom_boxplot() +
    ggpubr::stat_compare_means(
      comparisons = list(c("Spatial", "snRNAseq"), c("Spatial", "scRNAseq")),
      label = "p.signif"
    ) +
    scale_fill_brewer(palette = "Set3", name = "Dataset") +
    scale_y_log10() +
    facet_wrap(~ cellType + Type, nrow = 1, scales = "free") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Immune Cell Proportions by Dataset",
      x = "Dataset",
      y = "Proportion"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.position = "none"
    )
)
dev.off()

# Check glomerular sparsity
cellProportions2 <- data.frame(cells = names(clust), cluster = clust, tissue = annot$tissuename)

cellProportions2$broad <- case_when(
  cellProportions2$cluster %in% c("Podocyte", 
    "Glomerular.endothelium",
    "Mesangial.cell",
    "Parietal.epithelium") ~ "glomerular",
    TRUE ~ "nonGlomerular"
)

cellProportions2$Type <- case_when(
  str_detect(cellProportions2$tissue, "SLE") ~ "LN",
  TRUE ~ "Control"
)

glomProps <- cellProportions2 %>%
  group_by(tissue, broad, Type) %>%
  summarise(count = n()) %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count))

print(glomProps)

# Combine AMP dataset
myeloid <- readRDS("../../finalObjects_v3/myeloid_tissue_meta.Rds")  %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonGlomerular")
b <- readRDS("../../finalObjects_v3/b_tissue_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonGlomerular")
tnk <- readRDS("../../finalObjects_v3/t_nk_tissue_meta.Rds")  %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonGlomerular")
glom <- readRDS("../../finalObjects_v3/glom_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "glomerular")
loh <- readRDS("../../finalObjects_v3/loh_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonGlomerular")
pt <- readRDS("../../finalObjects_v3/pt_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonGlomerular")
intl <- readRDS("../../finalObjects_v3/intl_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonGlomerular")
dn <- readRDS("../../finalObjects_v3/dn_meta.Rds") %>% select(cell, sample, dataset, Type) %>% mutate(broad = "nonGlomerular")

amp <- rbind(myeloid, b, tnk, glom, loh, pt, intl, dn)

ampProps <- amp %>%
  group_by(sample, broad, dataset, Type) %>%
  summarise(count = n()) %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count))

allDanaher <- glomProps %>% select(tissue, broad, Type, proportion) %>% mutate(dataset = "Spatial")
allAmp <- ampProps %>% mutate(tissue = sample) %>%
    ungroup() %>% 
    select(tissue, broad, Type, proportion, dataset)

allSamples <- rbind(allDanaher, allAmp)

allSamplesImmune <- allSamples %>% filter(broad == "glomerular")
allSamplesImmune$dataset <- factor(allSamplesImmune$dataset, levels = c("Spatial", "scRNAseq", "snRNAseq"))

pdf("glomVsNonGlom.pdf", width = 8, height = 6)
    print(ggplot(allSamplesImmune, aes(x = dataset, y = proportion, fill = dataset)) +
        geom_boxplot() +
        ggpubr::stat_compare_means(
          comparisons = list(c("Danaher", "scRNAseq"), c("Danaher", "snRNAseq")),
          label = "p.signif"
        ) +
        scale_fill_brewer(palette = "Set3", name = "Dataset") +
        scale_y_log10(breaks=c(.01,.1,1)) +        
        facet_wrap(~ Type, scale = "fixed") +
        theme_minimal(base_size = 18) +
        labs(title = "Glomerular Cell Proportions by Dataset (All Cells)", x = "Dataset", y = "Proportion")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 18),       # facet labels
    axis.title = element_text(size = 18),       # axis labels
    legend.position = "none"
  )) 
dev.off()

# # As a proportion of non-immune cells
# amp <- rbind(glom, loh, pt, intl, dn)

# ampProps <- amp %>%
#   group_by(sample, broad, dataset, Type) %>%
#   summarise(count = n()) %>%
#   group_by(sample) %>%
#   mutate(proportion = count / sum(count))

# cellProportions3 <- data.frame(cells = names(clust), cluster = clust, tissue = annot$tissuename)

# cellProportions3 <- cellProportions3[!(cellProportions3$cluster %in% 
#   c("monocyte", "mast", "macrophage", "neutrophil", "mDC", "pDC",
#   "B-cell", "plasmablast",
#   "T CD4 naive", "T CD8 memory", "T CD8 naive", "T CD4 memory", "NK", "Treg")), ]

# cellProportions3$Type <- case_when(
#   str_detect(cellProportions3$tissue, "SLE") ~ "LN",
#   TRUE ~ "Control"
# )

# cellProportions3$broad <- case_when(
#   cellProportions3$cluster %in% c("Podocyte", 
#     "Glomerular.endothelium",
#     "Mesangial.cell",
#     "Parietal.epithelium") ~ "glomerular",
#     TRUE ~ "nonGlomerular"
# )

# glomProps <- cellProportions3 %>%
#   group_by(tissue, broad, Type) %>%
#   summarise(count = n()) %>%
#   group_by(tissue) %>%
#   mutate(proportion = count / sum(count))

# print(glomProps)

# allDanaher <- glomProps %>% select(tissue, broad, Type, proportion) %>% mutate(dataset = "Danaher")

# allSamples <- rbind(allDanaher, allAmp)

# allAmp <- ampProps %>% mutate(tissue = sample) %>%
#     ungroup() %>% 
#     select(tissue, broad, Type, proportion, dataset)

# allSamples <- rbind(allDanaher, allAmp)

# allSamplesImmune <- allSamples %>% filter(broad == "glomerular")

# pdf("glomVsNonGlom_nonImmune.pdf", width = 8, height = 6)
#     print(ggplot(allSamplesImmune, aes(x = dataset, y = proportion, fill = dataset)) +
#         geom_boxplot() +
#         ggpubr::stat_compare_means(
#           comparisons = list(c("Danaher", "scRNAseq"), c("Danaher", "snRNAseq")),
#           label = "p.signif"
#         ) +
#         scale_fill_brewer(palette = "Set3", name = "Dataset") +
#         scale_y_log10(breaks=c(.01,.1,1)) +        
#         facet_wrap(~ Type, scale = "fixed") +
#         theme_minimal(base_size = 18) +
#         labs(title = "Glomerular Cell Proportions by Dataset (Nonimmune)", x = "Dataset", y = "Proportion")+
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.text = element_text(size = 18),       # facet labels
#     axis.title = element_text(size = 18),       # axis labels
#     legend.position = "none"
#   )) 
# dev.off()
