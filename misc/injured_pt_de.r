library(Seurat)
library(dplyr) 
library(tidyr)
library(viridis)
library(stringr)
library(pheatmap)
library(presto)
library(pals)
library(harmony)
library(cowplot)
library(singlecellmethods)  
library(ggplot2)
library(lisi)
source("/data/srlab/ik936/Foxxy/utils/utils.R")
source("/data/srlab/anathan/scripts/scseq_utils.R")
library(parallel)
library('lme4', lib.loc = "/PHShome/ssg34/.conda/envs/jupy/lib/R/library")
library(glmmTMB, lib.loc ="/PHShome/ssg34/.conda/envs/jupy/lib/R/library")
.libPaths("/PHShome/ssg34/.conda/envs/plswork/lib/R/library")
set.seed(0)

norm <- readRDS('/data/srlab2/qxiao/AMP-SLE/sc_nuc_data/normalized_data_ScNuc_500nGene_1000nUMI_3pctnontargetMTwdoubletandsampleFINAL-8-10-22.rds')
all_meta <- readRDS('/data/srlab2/qxiao/AMP-SLE/sc_nuc_data/20230215_Tissue/meta_Tissue_500nGene_1000nUMI_3pctnontargetMTwdoubletandsampleFINAL-02-22-23.rds')
pt_meta <- readRDS('/data/srlab2/qxiao/AMP-SLE/sc_nuc_data/20230215_Tissue/PT/2023-06-30_meta_500nGene_1000nUMI_3.rds')
pt_meta <- pt_meta %>% left_join(all_meta %>% select(cell, cluster_names))
injured_meta <- pt_meta %>% filter(cluster_names == 'injuredPT/DTL', dataset == 'scRNAseq')

pseudobulk <- function(ind, meta, norm) {
    cells <- meta %>% filter(sample == ind) %>% pull(cell)
    pb <- rowMeans(norm[, cells])
    pb <- c(pb, ind)
    names(pb) <- c(rownames(norm), 'sample')
    return(pb)
}



out <- lapply(unique(injured_meta$sample), pseudobulk, injured_meta, norm)
out <- bind_rows(out)
out <- left_join(out, injured_meta %>% select(sample, Type) %>% unique())

sample_stats <- injured_meta %>% 
                        select(sample) %>% 
                        table() %>% data.frame() %>% 
                            rename('sample' = '.',ncells = Freq) %>% 
                            left_join(injured_meta %>% group_by(sample) %>% 
                                        summarize(mean_ncount = mean(nCount_RNA)))

out <- left_join(out, sample_stats)


de <- function(feature, df) {
    model_df <- df %>% select(feature, sample, Type) %>% rename(Exp = feature)
    if (sum(model_df$Exp > 0) > 0.05 * nrow(df)) { 
        m_0 <- lm(Exp ~ ncells + mean_ncount, data = model_df)
        m_1 <- lm(Exp ~ ncells + mean_ncount + Type, data = model_df)
        ANNO <- anova(m_0, m_1)
        LRP <- ANNO[2,6]
        F <- ANNO[2,5]
        Beta <- summary(m_1)$coefficients['TypeLN', 'Estimate']
        SE <- summary(m_1)$coefficients['TypeLN', 'Std. Error']
        res <- c(gene = feature, LRP = LRP, F = F, Beta = Beta, SE = SE)
    }
    }
    else {
        res <- c(gene = gene, LRP = NA, F = NA, Beta = NA, SE = NA)
    }
    return(res)
}

out <- lapply(rownames(norm), de, out)
saveRDS(out, '/data/srlab/ssg34/SLE_kidney_v2/data/tissue/injured_pt_de.rds')


