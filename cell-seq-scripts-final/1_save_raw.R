#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")

#-------------------------------------------------------

filepath1 <- "/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/"
filepath2 <- "/outs/filtered_feature_bc_matrix/"
sample <- commandArgs(trailingOnly=TRUE) 

#-------------------------------------------------------

## Read in raw data and rename columns to include sample name
    raw_data <- Read10X(data.dir = paste(filepath1, sample, filepath2, sep = ""))
    colnames(raw_data) <- paste(colnames(raw_data), rep(sample, length(colnames(raw_data))), sep = "-")

## Make Raw Meta
    meta_raw <- data.frame(cells = colnames(raw_data),
                       nGene = colSums(raw_data > 0),
                       nUMI = colSums(raw_data),
                       sample = sample
                       )
## Identify mitochondrial genes
    mito_genes <- rownames(raw_data)[grep("^Mt-", rownames(raw_data),
                                                        value = FALSE, ignore.case = TRUE)]
## Calculate % Mitochondrial UMI
    meta_raw$percent.mito <- Matrix::colSums(raw_data[mito_genes, ])/Matrix::colSums(raw_data)

## Calculate $ Mitochondrial excluding DASHed MT-genes
    mito_genes_include <- c("MT-ND5", "MT-ND6")
    mito_genes_exclude <- c('MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP6','MT-ATP8','MT-CO3','MT-ND3','MT-ND4','MT-ND4L','MT-CYB')
    tmp <- raw_data[!(rownames(raw_data) %in% mito_genes_exclude),]
    meta_raw$percent.mito.sub <- Matrix::colSums(tmp[mito_genes_include, ])/Matrix::colSums(tmp)
    tmp <- raw_data[!(rownames(raw_data) %in% mito_genes),]
    meta_raw$nUMInonMT <- colSums(tmp)

## Save metadata
    saveRDS(meta_raw, paste("../cell-seq-output/raw_metadata/", sample, "_raw_meta.rds", sep = ""))

## Save raw data
    saveRDS(raw_data, paste("../cell-seq-output/raw_data/", sample, "_raw_data.rds", sep = ""))
