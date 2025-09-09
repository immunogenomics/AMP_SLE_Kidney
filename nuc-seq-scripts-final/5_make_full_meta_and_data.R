#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")

#-------------------------------------------------------

filepath1 <- "/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/"
filepath2 <- "/outs/raw_feature_bc_matrix/"
samples <- read.delim("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/list_of_nuc_samples.txt",
                            sep = "\n", header = F)
samples <- as.character(samples[,1])
path <- "/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis"

#-------------------------------------------------------

## FILTERING PARAMETERS
### NGENE : 500
### MT    : 5%

nGenepar <- 500
pctmtpar <- 0.01
nGeneupperpar <- 7500
nUMIpar <- 1000
nUMIupperpar <- 40000

#-------------------------------------------------------


meta_full <- c()
data_full <- c()


for (i in 1:length(samples)) {

## Read in filtered data
    data <- readRDS(paste("../nuc-seq-output/filtered_data/",
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
                                 sep = ""))


## Read in metadata
    meta <- readRDS(paste("../nuc-seq-output/filtered_metadata/",
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

## Add to full meta and data
    
    data_full <- cbind(data_full, data)
    meta_full <- rbind(meta_full, meta)
    
}

print("Dimensions of full meta and data...")
print(dim(meta_full))
print(dim(data_full))

## Merge meta with map
cell_map <- readRDS("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/outside_meta/nuc_iso_map.rds")
meta_full <- merge(meta_full, cell_map, by = "sample", all.x = T)
print("Have we added columns?")
print(dim(meta_full))


saveRDS(meta_full, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_metadata/",
                            "meta_all_filtered_",
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
saveRDS(data_full, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_data/",
                           "data_all_filtered_",
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

rm(meta_full)
rm(data_full)
rm(meta)
rm(data)
rm(samples)
rm(nGenepar)

gc()


