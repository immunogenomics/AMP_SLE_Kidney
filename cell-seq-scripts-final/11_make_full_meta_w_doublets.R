#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")

#-------------------------------------------------------

samples <- read.delim("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/list_of_cell_samples.txt",
                            sep = "\n", header = F)
samples <- as.character(samples[,1])

#-------------------------------------------------------

## FILTERING PARAMETERS
### NGENE : 500
### NUMI : 1000
### MT : 3% non-targeted

nGenepar <- 500
nUMIpar <- 1000
pctmt <- 3

#-------------------------------------------------------


meta_full <- c()


for (i in 1:length(samples)) {

## Read in metadata
    meta <- readRDS(paste("../cell-seq-output/filtered_metadata/",
                                 samples[i],
                                 "_meta_filtered_",
                                 nGenepar,
                                 "nGene_",
                                 nUMIpar,
                                 "nUMI_",
                                 pctmt,
                                 "pctnontargetMTwdoublets.rds",
                                 sep = "")
                   )

## Add to full meta and data
    
    meta_full <- rbind(meta_full, meta)
    
}

saveRDS(meta_full, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                                 "meta_filtered_all",
                                 nGenepar,
                                 "nGene_",
                                 nUMIpar,
                                 "nUMI_",
                                 pctmt,
                                 "pctnontargetMTwdoublets.rds",
                                 sep = "")
       )


rm(meta_full)
rm(meta)
rm(samples)
rm(nGenepar)

gc()


