#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")

#-------------------------------------------------------

## FILTERING PARAMETERS
### NGENE : 500
### MT    : NONE

nGenepar <- 500

#-------------------------------------------------------

meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_metadata/",
                            "meta_all_filtered_",
                            nGenepar,
                            "nGene.rds",
                            sep = "")
               )

harmony <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/harmony_output/",
                            "harmony_filtered_",
                            nGenepar,
                            "nGene.rds",
                            sep = "")
                   )

umap_res <- umap$UMAP(n_neighbors = 30L, metric = "euclidean", min_dist = .3)$fit_transform(harmony)
meta$hUMAP1 <- umap_res[, 1]
meta$hUMAP2 <- umap_res[, 2]

saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_metadata/",
                            "meta_all_filtered_",
                            nGenepar,
                            "nGene.rds",
                            sep = "")
       )

rm(harmony)
rm(meta)

gc()
