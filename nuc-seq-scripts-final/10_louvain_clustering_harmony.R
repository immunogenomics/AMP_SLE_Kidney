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

## Build SNNGraph and calculate louvain clustering at various resolutions
SNNGraph <- BuildSNNSeurat(harmony, k.param = 30, nn.eps = 0)

resolution_list <- c(0.1, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0)
ids_ref_cca <- Reduce(cbind, mclapply(resolution_list, function(res_use) {
    Seurat:::RunModularityClustering(SNN = SNNGraph, modularity = 1, 
        resolution = res_use, algorithm = 1, n.start = 20, 
        n.iter = 20, random.seed = 100, print.output = FALSE, 
        temp.file.location = NULL, edge.file.name = NULL)
}, mc.preschedule = FALSE, mc.cores = min(20, length(resolution_list))))
ids_ref_cca %<>% data.frame()
print("Removing SNNGraph")
rm(SNNGraph)
gc()
colnames(ids_ref_cca) <- sprintf("res_%.2f", resolution_list)

meta$hres.0.1 <- ids_ref_cca$res_0.10
meta$hres.0.1 <- factor(meta$hres.0.1, levels = sort(unique(meta$hres.0.1)))
meta$hres.0.4 <- ids_ref_cca$res_0.40
meta$hres.0.4 <- factor(meta$hres.0.4, levels = sort(unique(meta$hres.0.4)))
meta$hres.0.6 <- ids_ref_cca$res_0.60
meta$hres.0.6 <- factor(meta$hres.0.6, levels = sort(unique(meta$hres.0.6)))
meta$hres.0.8 <- ids_ref_cca$res_0.80
meta$hres.0.8 <- factor(meta$hres.0.8, levels = sort(unique(meta$hres.0.8)))
meta$hres.1.2 <- ids_ref_cca$res_1.20
meta$hres.1.2 <- factor(meta$hres.1.2, levels = sort(unique(meta$hres.1.2)))
meta$hres.1.6 <- ids_ref_cca$res_1.60
meta$hres.1.6 <- factor(meta$hres.1.6, levels = sort(unique(meta$hres.1.6)))
meta$hres.2.0 <- ids_ref_cca$res_2.00
meta$hres.2.0 <- factor(meta$hres.2.0, levels = sort(unique(meta$hres.2.0)))

## Save everything
saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_metadata/",
                            "meta_all_filtered_",
                            nGenepar,
                            "nGene.rds",
                            sep = "")
       )

rm(harmony)
rm(meta)
rm(ids_ref_cca)
gc()
