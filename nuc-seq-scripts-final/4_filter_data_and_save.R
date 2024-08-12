#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")

#-------------------------------------------------------

filepath1 <- "/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/"
filepath2 <- "/outs/filtered_feature_bc_matrix/"
samples <- read.delim("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/list_of_nuc_samples.txt",
                            sep = "\n", header = F)
samples <- as.character(samples[,1])

path <- "/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis"

#-------------------------------------------------------

## FILTERING PARAMETERS
### NGENE : 500
### MT    : 1%

nGenepar <- 500
pctmtpar <- 0.01
nGeneupperpar <- 7500
nUMIpar <- 1000
nUMIupperpar <- 40000

#-------------------------------------------------------


for (i in 1:length(samples)) {

## Read in raw data
    raw_data <- Read10X(data.dir = paste(filepath1, samples[i], filepath2, sep = ""))
    colnames(raw_data) <- paste(colnames(raw_data), rep(samples[i], length(colnames(raw_data))), sep = "-")

## Read in metadata
    meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/raw_metadata/",
                          samples[i], "_filtered_meta.rds", sep = ""))

## Filter meta and data
    filtered <- filter(meta, nGene > nGenepar, percent.mito < pctmtpar, nGene < nGeneupperpar, nUMI > nUMIpar, nUMI < nUMIupperpar)
    filtered_data <- raw_data[,colnames(raw_data) %in% filtered$cell]
    
## Save meta and data

    saveRDS(filtered, paste("../nuc-seq-output/filtered_metadata/",
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
    saveRDS(filtered_data, paste("../nuc-seq-output/filtered_data/",
                                 samples[i],
                                 "_data_filtered_",
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


}

rm(raw_data)
rm(meta)
rm(filtered)
rm(filtered_data)
gc()

