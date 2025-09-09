#@Script: SC T/NK Case ncorr Gene Correlation Permutation
#@Date: 05/10/2023
#@Author: Sid Gurajala 

library(dplyr) 
library(tidyr)
library(Matrix)
library(parallel)


qcd_norm <- readRDS('/data/srlab/ssg34/SLE_kidney_v2/data/qcd/T_NK_clusterQCd_norm_09092022.rds')
keep_genes <- readRDS('/data/srlab/ssg34/SLE_kidney_v2/data/DE/keep_genes_TNK_5.rds')

meta <- read.csv('/data/srlab/ssg34/SLE_kidney_v2/data/cna/t_nk/sc_t_nk_casecontrol_meta.csv')
ncorr <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/data/cna/t_nk/sc_t_nk_casecontrol.ncorrs.csv", header = FALSE)
ncorr <- ncorr[, 1]

meta_tmp_1 <- meta
meta_tmp_1$ncorr <- as.numeric(ncorr)

gene_ncorr_perm <- function(gene) {
    df <- data.frame(cbind(ncorr = meta_tmp_1$ncorr, exp = qcd_norm[gene, meta_tmp_1$cell]))
    cors <- vector()
    i <- 1
    set.seed(0)
    while (i < 10001) {
        perm_expr <- df$exp[sample(nrow(df))]
        cors <- c(cors, cor.test(perm_expr, df$ncorr)$estimate)
        i <- i + 1
    }
    statistic <- cor.test(df$exp, df$ncorr)$estimate
    if(statistic > min(cors) & statistic < max(cors)) {
        cors_meanshift <- cors - mean(cors)
        statistic_meanshift <- cor.test(df$exp, df$ncorr)$estimate - mean(cors)
        ppval <- sum(abs(cors_meanshift) >= abs(statistic_meanshift)) / length(cors_meanshift)
    }
    else {
        ppval <- 1/10001
    }
    return(ppval)
}


res <- mclapply(keep_genes, FUN = gene_ncorr_perm, mc.cores = 20)

out <- data.frame(gene = keep_genes, ppval = do.call(rbind, res))

saveRDS(out, "/data/srlab/ssg34/SLE_kidney_v2/data/DE/t_nk_DE_expr_ncorr_perm.rds")
