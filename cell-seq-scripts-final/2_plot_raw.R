#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")

#-------------------------------------------------------

filepath1 <- "/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/"
filepath2 <- "/outs/raw_feature_bc_matrix/"
samples <- read.delim("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/list_of_cell_samples.txt",
                            sep = "\n", header = F)
samples <- as.character(samples[,1])

#-------------------------------------------------------

for (i in 1:length(samples)) {
   
## Read in metadata
    meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/raw_metadata/",
                          samples[i], "_raw_meta.rds", sep = ""))

## Plot raw data and save   
    plot_raw(meta, samples[i])

}
