#@Script: T/NK Case Control MASC PBMC
#@Author: Sid Gurajala
#@Date:04/19/2024

set.seed(0)
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
library(ggrastr, lib.loc = "/PHShome/ssg34/.conda/envs/plswork/lib/R/library")
#library(Cairo)
source("/data/srlab/anathan/scripts/scseq_utils.R")
library(parallel)
library('lme4', lib.loc = "/PHShome/ssg34/.conda/envs/jupy/lib/R/library")
library(glmmTMB, lib.loc ="/PHShome/ssg34/.conda/envs/jupy/lib/R/library")
.libPaths("/PHShome/ssg34/.conda/envs/plswork/lib/R/library")

meta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/data/pbmc/t_nk/case_control/meta.csv")
model_df <- meta %>% select(Unified_Visit, Type, nCount_RNA, percent.mt, RNA_snn_res.1)
model_df <- model_df %>% 
                mutate(
                       nCount_RNA = log(nCount_RNA))

res <- MASC.me(model_df, factor(model_df$RNA_snn_res.1),
                contrast = "Type",
                random_effects = c("Unified_Visit"),
                fixed_effects = c("nCount_RNA", "percent.mt"),
                verbose = TRUE,
                save_models = F) %>%
            dplyr::mutate(bonferroni = p.adjust(model.pvalue, method = "bonferroni")) %>%
            dplyr::arrange(model.pvalue)

saveRDS(res, '/data/srlab/ssg34/SLE_kidney_v2/data/masc/tnk_pbmc/case_control.rds')