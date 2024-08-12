#@Script: Myeloid Chronicity DE Poisson Single-Cell Model
#@Author: Sid Gurajala
#@Date:05/26/2023


#Load in Libraries
library(dplyr)
library(parallel)
library('lme4', lib.loc = "/PHShome/ssg34/.conda/envs/jupy/lib/R/library")
library('glmmTMB', lib.loc ="/PHShome/ssg34/.conda/envs/jupy/lib/R/library")


#Assign Start and End Index 
args <- commandArgs(trailingOnly=TRUE)

start <- as.numeric(args[1])
end <- as.numeric(args[2])



myeloid_rawcounts <- readRDS('/data/srlab/ssg34/SLE_kidney_v2/data/qcd/Myeloid_clusterQCd_rawcounts_03102023.rds')
meta <- read.csv('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/myeloid/chronicity/sc_meta.csv')

if(end > nrow(myeloid_rawcounts)) {
    end <- nrow(myeloid_rawcounts)
}

mito_genes_subset <- c("MT-ND5", "MT-ND6")
percent_mito <- (colSums(myeloid_rawcounts[mito_genes_subset, ]) / colSums(myeloid_rawcounts)) * 100
nCount_RNA <- colSums(myeloid_rawcounts)
cellstats <- data.frame(cell = names(percent_mito), percent.mt = percent_mito, nCount_RNA = nCount_RNA)
meta <- meta %>% left_join(cellstats)

poisson_model_stats <- function(gene, meta, rawcounts) {
df <- meta %>% select(nCount_RNA, percent.mt, sample, First_biop, Final_Chronicity, Responder_Status) %>% 
                mutate(Exp = rawcounts[gene, meta$cell])
if(sum(df$Exp > 0) > 0.05 * nrow(df)) { 
H0 <- glmmTMB(Exp ~ log(nCount_RNA) + percent.mt + Responder_Status + First_biop + (1|sample), 
                     data = df, family=poisson(link='log'), 
                     control = glmmTMBControl(parallel = 10))
H1 <- glmmTMB(Exp ~ log(nCount_RNA) + percent.mt + Responder_Status + First_biop + (1|sample) + Final_Chronicity, 
                     data = df, family=poisson(link='log'), 
                     control = glmmTMBControl(parallel = 10))
ANNO <- anova(H0, H1)
LRP <- ANNO[["Pr(>Chisq)"]][2]
LRchisq<- ANNO[["Chisq"]][2]
Beta <- summary(H1)$coefficients$cond["Final_Chronicity", "Estimate"]
SE <- summary(H1)$coefficients$cond["Final_Chronicity", "Std. Error"]

stats <- c(gene = gene, LRP = LRP, LRChisq = LRchisq, Beta = Beta, SE = SE)
}
else {
stats <- c(gene = gene, LRP = NA, LRChisq = NA, Beta = NA, SE = NA)
}
return(stats) }

res <- lapply(rownames(myeloid_rawcounts)[start:end], FUN=poisson_model_stats, meta, myeloid_rawcounts)

#Save Results 
saveRDS(res, paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/myeloid/chronicity/Chronicity_DE_', start, '_', end, '.rds'))