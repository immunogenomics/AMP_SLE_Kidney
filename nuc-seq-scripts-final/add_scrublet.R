#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")
library(RcppCNPy)
#-------------------------------------------------------

filepath1 <- "/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/"
filepath2 <- "/outs/raw_feature_bc_matrix/"
samples <- read.delim("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/list_of_nuc_samples.txt",
                            sep = "\n", header = F)
samples <- as.character(samples[,1])

#-------------------------------------------------------

## FILTERING PARAMETERS
### NGENE : 500
### MT    : NONE

nGenepar <- 500
pctmtpar <- 0.01
nGeneupperpar <- 7500
nUMIpar <- 1000
nUMIupperpar <- 40000

#-------------------------------------------------------


for (i in 1:length(samples)) {

## Read in metadata
    meta <- readRDS(paste("../nuc-seq-output/full_metadata/",
                            samples[i],
                            "_meta_filtered_",
                            nGenepar,
                            "nGene_",
                            pctmtpar,
                            "pctmtUMI",
			    nGeneupperpar,
			    "uppernGene_",
			    nUMIpar,
			    "nUMI_",
			    nUMIupperpar,
			    "uppernUMI.rds",
                            sep = ""))

    doublet_score <- npyLoad(paste("../nuc-seq-output/python_data/doublet_score_", samples[i], ".npy", sep = ""))

    meta$doublet_score_scrublet <- doublet_score


    map <- readRDS("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/scrublet_output/scrublet_map.rds")

    thresh <- map[map$samples == samples[i],]$thresholds
    meta$isscrubletdoublet <- "no"
    meta[meta$doublet_score_scrublet > thresh,]$isscrubletdoublet <- "yes"

## Save everything
saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_metadata/",
                            samples[i],
                            "_meta_filtered_",
                            nGenepar,
                                 "nGene_",
                                 pctmtpar,
                                 "pctmtUMI",
                                 nGeneupperpar,
                                 "uppernGene_",
                                 nUMIpar,
                                 "nUMI_",
                                 nUMIupperpar,
                                 "uppernUMI.rds",
                            sep = "")
       )

rm(meta)
rm(doublet_score)
gc()

}


