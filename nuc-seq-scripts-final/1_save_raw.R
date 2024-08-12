#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")

#-------------------------------------------------------

filepath1 <- "/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/"
filepath2 <- "/outs/raw_feature_bc_matrix/"
samples <- read.delim("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/list_of_nuc_samples.txt",
                            sep = "\n", header = F)
samples <- as.character(samples[,1])
path <- "/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis"

#-------------------------------------------------------

for (i in 1:length(samples)) {

## Read in raw data and rename columns to include sample name
    raw_data <- Read10X(data.dir = paste(filepath1, samples[i], filepath2, sep = ""))
    colnames(raw_data) <- paste(colnames(raw_data), rep(samples[i], length(colnames(raw_data))), sep = "-")

## Make Raw Meta
    meta_raw <- data.frame(cells = colnames(raw_data),
                       nGene = colSums(raw_data > 0),
                       nUMI = colSums(raw_data),
                       gene.to.mol.ratio = colSums(raw_data > 0)/colSums(raw_data),
		       sample = samples[i]
                       )
## Identify mitochondrial genes
    mito_genes <- rownames(raw_data)[grep("^Mt-", rownames(raw_data),
                                                        value = FALSE, ignore.case = TRUE)]
## Calculate % Mitochondrial UMI
    meta_raw$percent.mito <- Matrix::colSums(raw_data[mito_genes, ])/Matrix::colSums(raw_data)

## Save metadata
    saveRDS(meta_raw, paste(path, "/nuc-seq-output/raw_metadata/", samples[i], "_raw_meta.rds", sep = ""))
}
