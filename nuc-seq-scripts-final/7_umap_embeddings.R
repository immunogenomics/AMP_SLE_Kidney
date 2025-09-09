#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")
path <- "/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis"

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

pca_res <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            "pca_cell_embeddings_filtered_",
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

meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_metadata/",
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

umap_res <- umap$UMAP(n_neighbors = 30L, metric = "euclidean", min_dist = .3)$fit_transform(pca_res[,1:20])
# umap_res <- umap$UMAP(n_neighbors = 30L, metric = "correlation", min_dist = .3)$fit_transform(harmony)
meta$UMAP1 <- umap_res[, 1]
meta$UMAP2 <- umap_res[, 2]

saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_metadata/",
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

rm(pca_res)
rm(meta)

gc()
