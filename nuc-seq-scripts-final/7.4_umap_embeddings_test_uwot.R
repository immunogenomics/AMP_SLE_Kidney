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
                            "uppernUMI_doubletremoval.rds",
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
                            "uppernUMI_doubletremoval.rds",
                            sep = ""))

umap_object <- uwot::umap(
    X = pca_res[,1:20],
    ret_extra = c('nn', 'fgraph', 'model')
) 
meta$final_UMAP1 <- umap_object$embedding[, 1]
meta$final_UMAP2 <- umap_object$embedding[, 2]

saveRDS(umap_object, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            "umap_all_filtered_",
                            nGenepar,
                            "nGene_",
                            pctmtpar,
                            "pctmtUMI",
                            nGeneupperpar,
                            "uppernGene_",
                            nUMIpar,
                            "nUMI_",
                            nUMIupperpar,
                            "uppernUMI_doubletremoval.rds",
                            sep = ""))

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
                            "uppernUMI_doubletremoval.rds",
                            sep = ""))

rm(pca_res)
rm(meta)
rm(umap_object)
gc()
