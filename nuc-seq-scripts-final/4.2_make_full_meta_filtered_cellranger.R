#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")

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


meta_full <- c()
data_full <- c()


for (i in 1:length(samples)) {


## Read in metadata
    meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/raw_metadata/",
                            samples[i],
                            "_filtered_meta.rds",
                            sep = "")
                  )
## Add to full meta and data
    
    meta_full <- rbind(meta_full, meta)
    
}

print("Dimensions of full meta and data...")
print(dim(meta_full))

cell_map <- readRDS("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/outside_meta/nuc_iso_map.rds")
meta_full <- merge(meta_full, cell_map, by = "sample", all.x = T)


saveRDS(meta_full, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_metadata/",
                            "meta_all_filtered.rds",
                            sep = ""))

rm(meta_full)
rm(meta)

gc()


