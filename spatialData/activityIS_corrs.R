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
library(wesanderson)

infilScores <- readRDS("MacrophagesOnly/macrophagesOnly_InfilScore.Rds")
cnaCorrs <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/CNA_activityRerun/cna_results/myeloid/Final_Activity_SiteFBChron_ncorr.csv", header = FALSE)
fdr <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/CNA_activityRerun/cna_results/myeloid/Final_Activity_SiteFBChron_fdrs.csv", header = FALSE)
cnaMeta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/myeloid/activity/sc_meta.csv")

colnames(fdr) <- c('threshold', 'fdr', 'ncells')
threshold <- fdr %>% filter(fdr < 0.1) %>% 
          mutate(fdr = round(fdr, 4)) %>% 
          filter(threshold == min(threshold)) %>% 
          pull(threshold)

cnaMeta$ncorr <- cnaCorrs$V1

macsOnlyCNA <- cnaMeta %>%
  filter(cell %in% infilScores$cell) %>% select(cell, Final_Activity, ncorr)

plotDf <- merge(infilScores %>% select(cell, final_annotation, cell_group, SAMorNot, InfilScore), macsOnlyCNA, by = "cell")

plotDf <- plotDf %>% mutate(status = case_when(
    ncorr > threshold ~ "Expanded",
    ncorr < -1 * threshold ~ "Depleted",
    TRUE ~ "Unchanged"
))

# Calculate correlation
cor_test <- cor.test(plotDf$InfilScore, plotDf$ncorr)
cor_val <- signif(cor_test$estimate, 3)
p_val <- signif(cor_test$p.value, 3)

# Prepare subtitle
subtitle_text <- paste0("Pearson r = ", cor_val, ", p = ", p_val)

pdf("MacrophagesOnly/activityIS_corrs.pdf", width = 8, height = 6)
    ggplot(plotDf, aes(x = ncorr, y = InfilScore, color = status)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", color = "black", se = FALSE) +
        labs(title = "Correlation between CNA ncorr and InfilScore in Macrophages",
        x = "CNA ncorr",
        y = "InfilScore",
        subtitle = subtitle_text) +
        geom_vline(xintercept = threshold, linetype="dotted", color = "grey") + 
        geom_vline(xintercept = -1 * threshold, linetype="dotted", color = "grey") + 
        scale_color_manual(values = c("Expanded" = "tomato", "Depleted" = "steelblue1", "Unchanged" = "lightgrey")) +
        theme_minimal()
dev.off()

# See if the same holds true with chronicity
cnaCorrs <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/CNA_chronicityCorrections_siteModified/cna_results/myeloid/Final_Chronicity_SiteFirstBiopsyResponse_ncorr.csv", header = FALSE)
fdr <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/CNA_chronicityCorrections_siteModified/cna_results/myeloid/Final_Chronicity_SiteFirstBiopsyResponse_fdrs.csv", header = FALSE)
cnaMeta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/myeloid/chronicity/sc_meta.csv")

colnames(fdr) <- c('threshold', 'fdr', 'ncells')
threshold <- fdr %>% filter(fdr < 0.1) %>% 
          mutate(fdr = round(fdr, 4)) %>% 
          filter(threshold == min(threshold)) %>% 
          pull(threshold)

cnaMeta$c_corr <- cnaCorrs$V1
macsOnlyCNA <- cnaMeta %>%
  filter(cell %in% infilScores$cell) %>% select(cell, Final_Activity, c_corr)

plotDf2 <- merge(plotDf, macsOnlyCNA, by = "cell")
plotDf2 <- plotDf2 %>% mutate(status = case_when(
    c_corr > threshold ~ "Expanded",
    c_corr < -1 * threshold ~ "Depleted",
    TRUE ~ "Unchanged"
))

# Calculate correlation
cor_test <- cor.test(plotDf2$InfilScore, plotDf2$c_corr)
cor_val <- signif(cor_test$estimate, 3)
p_val <- signif(cor_test$p.value, 3)

# Prepare subtitle
subtitle_text <- paste0("Pearson r = ", cor_val, ", p = ", p_val)

pdf("MacrophagesOnly/chronicityIS_corrs.pdf", width = 8, height = 6)
    ggplot(plotDf2, aes(x = c_corr, y = InfilScore, color = status)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", color = "black", se = FALSE) +
        labs(title = "Correlation between Chronicity CNA and InfilScore in Macrophages",
        x = "Chronicity CNA",
        y = "InfilScore",
        subtitle = subtitle_text) +
        geom_vline(xintercept = threshold, linetype="dotted", color = "grey") + 
        geom_vline(xintercept = -1 * threshold, linetype="dotted", color = "grey") + 
        scale_color_manual(values = c("Expanded" = "tomato", "Depleted" = "steelblue1", "Unchanged" = "lightgrey")) +
        theme_minimal()
dev.off()
