suppressPackageStartupMessages({
    library(Matrix)
    library(ggplot2)
    library(dplyr)
    library(tidyverse)
    library(wesanderson)
    library(pals)
    library(presto)
    library(fedmatch)
    library(plotly)
    library(gridExtra)
})

meta <- readRDS('/data/srlab/ssg34/SLE_kidney_v2/data/qcd/Myeloid_clusterQCd_meta_harmonizedPCUMAPCellStateClusters_10042022.rds')
norm <- readRDS('/data/srlab/ssg34/SLE_kidney_v2/data/qcd/Myeloid_clusterQCd_norm_10042022.rds')

de <- wilcoxauc(norm, meta$final_annotation)
df <- de %>% filter(abs(logFC) >= 0.5 & padj < 0.01)

cna_meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/myeloid/activity/sc_meta.csv'), header = TRUE)
cna_meta$corrs <- read.csv("cna_results/myeloid_tissue/Final_Activity_ChronicityBiopAndResponse_ncorr.csv", header = FALSE)$V1

fdrs <- read.csv("cna_results/myeloid_tissue/Final_Activity_ChronicityBiopAndResponse_fdrs.csv", header = FALSE)
activityThreshold <- min(fdrs[fdrs$V2 < 0.1, ]$V1)
cna_meta <- cna_meta[abs(cna_meta$corrs) > activityThreshold, ]

testCorr <- function() {
    gex <- norm[, colnames(norm) %in% cna_meta$cell]
    calculate_correlations <- function(gene_expression, activity) {
        cor_result <- cor.test(gene_expression, activity)
        c(correlation = cor_result$estimate, p_value = cor_result$p.value)
    }
    cellTypeCorr <- cna_meta$corrs[cna_meta$cell %in% colnames(gex)]
    correlations <- apply(gex, 1, calculate_correlations, 
                      activity = cellTypeCorr)
    rownames(correlations) <- c("corr", "p")
    correlations <- as.data.frame(t(correlations))
    correlations <- correlations[!(is.na(correlations$corr)) & abs(correlations$corr) > 0.3, ] 
    correlations_sorted <- correlations[order(-correlations$corr), ]
    write.csv(correlations_sorted, "activity_globalGeneCorrelations.csv")
}

testCorr()

# Now do the same for chronicity

cna_meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/myeloid/activity/sc_meta.csv'), header = TRUE)
cna_meta$chronicity <- read.csv("cna_results/myeloid_tissue/Final_Chronicity_justChecking_ncorr.csv", header = FALSE)$V1
chronicity_fdrs <- read.csv("cna_results/myeloid_tissue/Final_Chronicity_justChecking_fdrs.csv", header = FALSE)

chronThreshold <- min(chronicity_fdrs[chronicity_fdrs$V2 < 0.1, ]$V1)

cna_meta <- cna_meta[abs(cna_meta$chronicity) > chronThreshold, ]

testChronCorr <- function() {
    gex <- norm[, colnames(norm) %in% cna_meta$cell]
    calculate_correlations <- function(gene_expression, chronicity) {
        cor_result <- cor.test(gene_expression, chronicity)
        c(correlation = cor_result$estimate, p_value = cor_result$p.value)
    }
    cellTypeCorr <- cna_meta$chronicity[cna_meta$cell %in% colnames(gex)]
    correlations <- apply(gex, 1, calculate_correlations, 
                      chronicity = cellTypeCorr)
    rownames(correlations) <- c("corr", "p")
    correlations <- as.data.frame(t(correlations))
    correlations <- correlations[!(is.na(correlations$corr)) & abs(correlations$corr) > 0.3, ] 
    correlations_sorted <- correlations[order(-correlations$corr), ]
    write.csv(correlations_sorted, "chronicity_globalGeneCorrelations.csv")
}

testChronCorr()