#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")

#-------------------------------------------------------

## FILTERING PARAMETERS
### NGENE : 500
### NUMI : 1000
### MT : 3% non-targeted

nGenepar <- 500
nUMIpar <- 1000
pctmt <- 3

#-------------------------------------------------------

pca_res <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "pca_cell_embeddings_filtered_",
                            nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMT.rds",
                            sep = "")
)

meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_filtered_all_",
                            nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMT.rds",
                            sep = "")
)

## harmony

    harmony <- HarmonyMatrix(pca_res[,1:20], meta, vars_use = c("biopsy.ID", "processing.date", "processing.batch", "Site"), theta = c(0,0,0,0),
                         # lambda = 1, 
#                          tau = 0, 
#                          epsilon.cluster = -Inf,
#                          epsilon.harmony = -Inf,
                         max.iter.cluster = 30,
                         max.iter.harmony = 20,
                         plot_convergence = T,
                             npcs = 20,
                            do_pca = F)

colnames(harmony) <- c("hPC-1", "hPC-2", "hPC-3", "hPC-4",
                      "hPC-5", "hPC-6", "hPC-7", "hPC-8",
                      "hPC-9", "hPC-10", "hPC-11", "hPC-12",
                      "hPC-13", "hPC-14", "hPC-15", "hPC-16",
                      "hPC-17", "hPC-18", "hPC-19", "hPC-20")

print("Checking to make sure harmony PCs and metadata are in the same order...")
all(rownames(harmony) == meta$cells)

## Add harmony PCs to meta
meta <- cbind(meta, harmony)


saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_filtered_all_",
                            nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMT.rds",
                            sep = "")
)

saveRDS(harmony, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/harmony_output/",
                            "harmony_filtered_",
                            nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMT.rds",
			    sep = "")
       )

rm(pca_res)
rm(meta)
rm(harmony)

gc()
