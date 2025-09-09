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

meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_all_filtered_",
			    nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMTwdoubletandsampleFINAL.rds",
                            sep = "")
)

harmony <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/harmony_output/",
                            "harmony_filtered_",
                            nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMTwdoubletandsampleFINAL.rds",
			    sep = "")
)

umap_object <- uwot::umap(
    X = harmony,
    ret_extra = c('nn', 'fgraph', 'model')
) 
meta$huwotUMAP1 <- umap_object$embedding[, 1]
meta$huwotUMAP2 <- umap_object$embedding[, 2]

saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_all_filtered_",
			    nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMTwdoubletandsampleFINAL.rds",
                            sep = "")
)

saveRDS(umap_object, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "umap_harmony_all_filtered_",
                            nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMTwdoubletandsampleFINAL.rds",
                            sep = "")
)

rm(harmony)
rm(meta)

gc()
