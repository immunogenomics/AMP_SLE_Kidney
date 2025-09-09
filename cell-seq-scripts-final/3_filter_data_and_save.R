#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")

#-------------------------------------------------------

sample <- commandArgs(trailingOnly=TRUE) 

#-------------------------------------------------------

## FILTERING PARAMETERS
### NGENE : 500
### NUMI : 1000
### MT    :3% non-targeted

nGenepar <- 500
nUMIpar <- 1000
pctmt <- 3

#-------------------------------------------------------



## Read in raw data
    raw_data <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/raw_data/",
                              sample, "_raw_data.rds", sep = ""))

## Read in metadata
    meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/raw_metadata/",
                          sample, "_raw_meta.rds", sep = ""))

## Filter meta and data
    filtered <- filter(meta, nGene > nGenepar, nUMI > nUMIpar, percent.mito.sub < 0.03)
    filtered_data <- raw_data[,colnames(raw_data) %in% filtered$cell]
    

## Save meta and data

    saveRDS(filtered, paste("../cell-seq-output/filtered_metadata/",
                            sample,
                            "_meta_filtered_",
                            nGenepar,
                            "nGene_",
			    nUMIpar,
			    "nUMI_",
                            pctmt,
                            "pctnontargetMT.rds",
                            sep = "")
           )
    saveRDS(filtered_data, paste("../cell-seq-output/filtered_data/",
                                 sample,
                                 "_meta_filtered_",
                                 nGenepar,
                                 "nGene_",
                                 nUMIpar,
                                 "nUMI_",
                                 pctmt,
                                 "pctnontargetMT.rds",
                                 sep = "")           
           )


rm(raw_data)
rm(meta)
rm(filtered)
rm(filtered_data)
gc()

