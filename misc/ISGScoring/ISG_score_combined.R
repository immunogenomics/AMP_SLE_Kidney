suppressPackageStartupMessages({
    source("/data/srlab/anathan/scripts/scseq_utils.R")
    library(Matrix)
    library(MASS)
    library(tidyverse)
    library(ggplot2)
})

clinData <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects/clinicalData.rds")

# Immune Populations
bmeta <- readRDS("../../reorderedClusters/b_tissue_meta_reOrdered.Rds")
b_isg <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/misc/bp_tissue_ISGscoreseurat.rds")
bmeta$ISG <- b_isg$ISG
bmeta <- bmeta %>% select(sample, ISG)

tmeta <- readRDS("../../reorderedClusters/t_nk_tissue_meta_reOrdered.Rds")
t_isg <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/misc/tnk_tissue_ISGscoreseurat.rds")
tmeta$ISG <- t_isg$ISG
tmeta <- tmeta %>% select(sample, ISG)

mmeta <- readRDS("../../reorderedClusters/myeloid_tissue_meta_reOrdered.Rds")
m_isg <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/misc/myeloid_tissue_ISGscoreseurat.rds")
mmeta$ISG <- m_isg$ISG
mmeta <- mmeta %>% select(sample, ISG)

# Merge all meta data
allMeta <- rbind(bmeta, tmeta, mmeta)

pseudo <- allMeta %>%
  group_by(sample) %>%
  summarize(mean_isg = mean(ISG, na.rm = TRUE))

cases_immune <- pseudo[pseudo$sample %in% clinData$sample[clinData$Type == "LN"], ]
controls_immune <- pseudo[pseudo$sample %in% clinData$sample[clinData$Type == "Control"], ]

cases_immune <- merge(cases_immune, clinData[, c("sample", "Final_Chronicity", "Final_Activity")], by = "sample") %>%
  filter(!is.na(Final_Chronicity))

controls_immune$Final_Chronicity <- -1
controls_immune$Final_Activity <- -1

# Tissue
ptmeta <- readRDS("../../finalObjects/pt_meta_v2.Rds")
pt_isg <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/misc/pt_tissue_ISGscoreseurat.rds")
ptmeta$ISG <- ptmeta$ISG
ptmeta <- bmeta %>% select(sample, ISG)

glommeta <- readRDS("../../finalObjects/glom_meta.Rds")
glom_isg <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/misc/glom_tissue_ISGscoreseurat.rds")
glommeta$ISG <- glommeta$ISG
glommeta <- tmeta %>% select(sample, ISG)

dnmeta <- readRDS("../../finalObjects/dn_meta.Rds")
dn_isg <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/misc/dn_tissue_ISGscoreseurat.rds")
dnmeta$ISG <- dnmeta$ISG
dnmeta <- tmeta %>% select(sample, ISG)

intlmeta <- readRDS("../../finalObjects/intl_meta.Rds")
intl_isg <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/misc/intl_tissue_ISGscoreseurat.rds")
intlmeta$ISG <- intlmeta$ISG
intlmeta <- tmeta %>% select(sample, ISG)

lohmeta <- readRDS("../../finalObjects/loh_meta.Rds")
loh_isg <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/misc/loh_tissue_ISGscoreseurat.rds")
lohmeta$ISG <- lohmeta$ISG
lohmeta <- tmeta %>% select(sample, ISG)

# Merge all meta data
allMeta <- rbind(ptmeta, glommeta, dnmeta, intlmeta, lohmeta)

pseudo <- allMeta %>%
  group_by(sample) %>%
  summarize(mean_isg = mean(ISG, na.rm = TRUE))

cases_tissue <- pseudo[pseudo$sample %in% clinData$sample[clinData$Type == "LN"], ]
controls_tissue <- pseudo[pseudo$sample %in% clinData$sample[clinData$Type == "Control"], ]

cases_tissue <- merge(cases_tissue, clinData[, c("sample", "Final_Chronicity", "Final_Activity")], by = "sample") %>%
  filter(!is.na(Final_Chronicity))

controls_tissue$Final_Chronicity <- -1
controls_tissue$Final_Activity <- -1

# Combine immune and tissue data
cases_immune$Type <- "Immune"
controls_immune$Type <- "Immune"
cases_tissue$Type <- "Tissue"
controls_tissue$Type <- "Tissue"

cases <- rbind(cases_immune, cases_tissue)
controls <- rbind(controls_immune, controls_tissue)

p <- ggplot() + 
    geom_hline(yintercept = 0, color = "darkgrey", size = 0.5) +
    geom_vline(xintercept = 0, color = "darkgrey", size = 0.5) +
    # Points
    geom_point(data = cases, aes(x = Final_Chronicity, y = mean_isg, color = Type), shape = 21) +
    # Immune line
    geom_smooth(data = cases[cases$Type == "Immune", ], aes(x = Final_Chronicity, y = mean_isg), method = "lm", color = "#309343", se = FALSE) +
    # Tissue line
    geom_smooth(data = cases[cases$Type == "Tissue", ], aes(x = Final_Chronicity, y = mean_isg), method = "lm", color = "#ffc156", se = FALSE) +
    geom_boxplot(data = controls, aes(x = Final_Chronicity, y = mean_isg, fill = Type), alpha = 0.5) +
    xlab("Chronicity") + ylab("ISG Score") +
    scale_color_manual(values = c("#309343","#ffc156")) +
    scale_fill_manual(values = c("#309343","#ffc156")) +
    scale_x_continuous(breaks = seq(-2, 11, by = 1)) + 
    theme_minimal()

pdf("ISG_score_chronicity_combined.pdf", width = 7, height = 5)
    print(p)
dev.off()

p <- ggplot() + 
    geom_hline(yintercept = 0, color = "darkgrey", size = 0.5) +
    geom_vline(xintercept = 0, color = "darkgrey", size = 0.5) +
    # Points
    geom_point(data = cases, aes(x = Final_Activity, y = mean_isg, color = Type), shape = 21) +
    # Immune line
    geom_smooth(data = cases[cases$Type == "Immune", ], aes(x = Final_Activity, y = mean_isg), method = "lm", color = "#309343", se = FALSE) +
    # Tissue line
    geom_smooth(data = cases[cases$Type == "Tissue", ], aes(x = Final_Activity, y = mean_isg), method = "lm", color = "#ffc156", se = FALSE) +
    geom_boxplot(data = controls, aes(x = Final_Activity, y = mean_isg, fill = Type), alpha = 0.5) +
    xlab("Activity") + ylab("ISG Score") +
    scale_color_manual(values = c("#309343","#ffc156")) +
    scale_fill_manual(values = c("#309343","#ffc156")) +
    scale_x_continuous(breaks = seq(-2, 20, by = 1)) + 
    theme_minimal()

pdf("ISG_score_activity_combined.pdf", width = 7, height = 5)
    print(p)
dev.off()