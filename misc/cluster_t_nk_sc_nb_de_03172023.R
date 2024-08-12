#@Script:T_NK DE Poisson Single-Cell Model
#@Author: Sid Gurajala
#@Date:03/06/2023


#Load in Libraries
library(dplyr)
library(parallel)
library('lme4', lib.loc = "/PHShome/ssg34/.conda/envs/jupy/lib/R/library")
library('glmmTMB', lib.loc ="/PHShome/ssg34/.conda/envs/jupy/lib/R/library")

#Assign Start and End Index 
args <- commandArgs(trailingOnly=TRUE)

start <- as.numeric(args[1])
end <- as.numeric(args[2])
cluster <- as.numeric(args[3])

#load in files
qcd_meta <- readRDS('/data/srlab/ssg34/SLE_kidney_v2/data/qcd/T_NK_clusterQCd_meta_harmonizedPCUMAPclusters_annotations09212022.rds')
tnk_rawcounts <- readRDS('/data/srlab/ssg34/SLE_kidney_v2/data/qcd/T_NK_clusterQCd_rawcounts_03062023.rds')

#compute mito percentage
mito_genes_subset <- c("MT-ND5", "MT-ND6")
percent_mito <- colSums(tnk_rawcounts[mito_genes_subset, ]) / colSums(tnk_rawcounts) * 100
percent_mito <- data.frame(cell = names(percent_mito), percent.mt = percent_mito)


qcd_meta <- left_join(qcd_meta, percent_mito)


#subset to singlecell
sc_meta <- qcd_meta %>% filter(dataset == "scRNAseq") %>% filter(new_cluster_number == cluster)

#For last gene subset
if(end > nrow(tnk_rawcounts)) {
    end <- nrow(tnk_rawcounts)
}

#Poisson Single Cell Model Function
poisson_model_stats <- function(gene, meta, rawcounts) {
df <- meta %>% select(nCount_RNA, percent.mt, sample, Type) %>% 
                mutate(numeric_type = ifelse(Type == "Control", 0, 1)) %>% 
                mutate(Exp = rawcounts[gene, meta$cell])
if(sum(df$Exp > 0) > 0.025 * nrow(df)) { 
H0 <- glmmTMB(Exp ~ log(nCount_RNA) + percent.mt + (1|sample), 
                     data = df, family=poisson(link='log'), 
                     control = glmmTMBControl(parallel = 10))
H1 <- glmmTMB(Exp ~ log(nCount_RNA) + percent.mt + (1|sample) + numeric_type, 
                     data = df, family=poisson(link='log'), 
                      control = glmmTMBControl(parallel = 10))
ANNO <- anova(H0, H1)
LRP <- ANNO[["Pr(>Chisq)"]][2]
LRchisq<- ANNO[["Chisq"]][2]
Beta <- summary(H1)$coefficients$cond["numeric_type", "Estimate"]
SE <- summary(H1)$coefficients$cond["numeric_type", "Estimate"]

stats <- c(gene = gene, LRP = LRP, LRChisq = LRchisq, Beta = Beta, SE = SE)
}
else {
stats <- c(gene = gene, LRP = NA, LRChisq = NA, Beta = NA, SE = NA)
}
return(stats) }

#Run gene indices
res <- lapply(rownames(tnk_rawcounts)[start:end], FUN=poisson_model_stats, sc_meta, tnk_rawcounts)

#Save Results 
saveRDS(res, paste0('/data/srlab/ssg34/SLE_kidney_v2/data/DE/T_NK_DE_cluster', cluster, '_', start, '_', end, '.rds'))
