suppressPackageStartupMessages({
    library(dplyr)
    library(scales)
    library(pals)
    library(ggplot2)
    library(cowplot)
    library(tidyverse)
    library(extrafont)
    library(ggrastr)
})

args = commandArgs(trailingOnly=TRUE)

name <- args[1]
smallName <- args[2]
meta <- read.csv(
    paste0("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/", name, "_combinedMeta.csv"),
    header = TRUE
)

umap <- readRDS(paste0("../../harmonyIntegration_finalFinal/", name, "/harmony_0/", name, "_umap.Rds"))

counters <- read.csv(paste0("../../harmonyIntegration_finalFinal/", name, "/harmony_0/pbmcVsTissue.csv"))

umap_table <- data.frame(barcode = row.names(data.frame(umap$embedding)), 
    tissueScore = counters$tissue / (counters$tissue + counters$PBMC * (length(meta$origin[meta$origin == "Tissue"]) / 
        length(meta$origin[meta$origin == "PBMC"]))),
    annotation = meta$annotation,
    origin = meta$origin
)

smallMeta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', smallName, '/activity/sc_meta.csv'))
ncorr <- read.csv("../cna_results/myeloid/Final_Activity_SiteFBChron_ncorr.csv", header = FALSE)
fdrs <- read.csv("../cna_results/myeloid/Final_Activity_SiteFBChron_fdrs.csv", header = FALSE)

colnames(fdrs) <- c('threshold', 'fdr', 'ncells')
smallMeta$ncorr <- ncorr$V1

fdr <- fdrs %>% filter(fdr < 0.1) %>% 
            mutate(fdr = round(fdr, 4)) %>% 
            filter(fdr == max(fdr)) %>% 
            pull(threshold)

fullMeta <- merge(smallMeta, umap_table, by.x = "cell", by.y = "barcode")

fullMeta <- fullMeta %>% mutate(status = case_when(
    ncorr > fdr ~ "Expanded",
    ncorr < -1 * fdr ~ "Depleted",
    TRUE ~ "Unchanged"
  ))

fullMeta$status <- factor(fullMeta$status, levels = c("Depleted", "Unchanged", "Expanded"))

p<-ggplot(fullMeta, aes(x=status, y=tissueScore, fill=status)) +
  geom_violin(scale = "width") +
  scale_fill_manual(values = c("steelblue1", "lightgrey", "tomato")) + 
  geom_boxplot(width=0.1) + theme_classic() + 
  theme(
    legend.text = element_text(size = 14), 
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 14)
  )

pdf(paste0(smallName, "_vlns.pdf"), width = 6, height = 5)
    print(p)
dev.off()