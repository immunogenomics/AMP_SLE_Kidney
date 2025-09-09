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

#-------------------------------------------------------


for (i in 1:length(samples)) {

## Read in metadata
    meta <- readRDS(paste("../nuc-seq-output/full_metadata/",
                            samples[i],
                            "_meta_filtered_",
                            nGenepar,
                            "nGene.rds",
                            sep = ""))

# make umap mapping
	umap_res <- data.frame(X1 = meta$UMAP1,
                               X2 = meta$UMAP2)

# make variable for would pass 1000 gene threshold
	meta$pass <- "no"
	meta[meta$nGene > 1000,]$pass <- "yes"

# plot and save
	p1 <- plot_clusters_nuc(meta$res.0.6, meta$sample[1], umap_use = umap_res)
	p2 <- plot_clusters_nuc(meta$pass, meta$sample[1], umap_use = umap_res)

	ggsave(
      		paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/individual_plots/",
            		samples[i], "_res.0.6clusters.png", sep = ""),
      		plot = p1,
      		scale = 1,
      		width = 4,
      		height = 4,
      		dpi = 300,
    	)
	ggsave(
     	 	paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/individual_plots/",
            		samples[i], "_1000nGene.png", sep = ""),
      		plot = p1,
      		scale = 1,
      		width = 4,
      		height = 4,
      		dpi = 300,
    	)

rm(meta)
rm(umap)
gc()

}


