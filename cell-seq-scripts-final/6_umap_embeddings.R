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

umap_object <- uwot::umap(
    X = pca_res[,1:20],
    ret_extra = c('nn', 'fgraph', 'model')
) 
meta$uwotUMAP1 <- umap_object$embedding[, 1]
meta$uwotUMAP2 <- umap_object$embedding[, 2]

saveRDS(umap_object, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "umap_all_filtered_",
                            nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMT.rds",
                            sep = "")
)

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

rm(pca_res)
rm(meta)

gc()
