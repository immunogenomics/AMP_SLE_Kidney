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

## Read in data

data <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_data/",
                           "data_filtered_all_",
                            nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMTwdoubletandsampleFINAL.rds",
                            sep = "")
		)


sc_data <- CreateSeuratObject(data)
rm(data);gc()
sc_data <- NormalizeData(sc_data)
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 5000)
sc_data <- ScaleData(sc_data)
sc_data <- RunPCA(sc_data, features = VariableFeatures(object = sc_data))

saveRDS(sc_data@assays$RNA@data, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "normalized_data_all_filtered_",
                            nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMTwdoubletandsampleFINAL.rds",
                            sep = ""),
        version = 2
       )

saveRDS(sc_data@reductions$pca@cell.embeddings, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "pca_cell_embeddings_filtered_",
                            nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMTwdoubletandsampleFINAL.rds",
                            sep = ""),
        version = 2
       )

saveRDS(sc_data@assays$RNA@var.features, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "variable_genes_filtered_",
                            nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMTwdoubletandsampleFINAL.rds",
                            sep = ""),
        version = 2
       )

saveRDS(sc_data@reductions$pca@feature.loadings, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "feature_loadings_filtered_",
                            nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMTwdoubletandsampleFINAL.rds",
                            sep = ""),
        version = 2
       )

saveRDS(sc_data, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "seurat_object_filtered_",
                            nGenepar,
                            "nGene_",
                            nUMIpar,
                            "nUMI_",
                            pctmt,
                            "pctnontargetMTwdoubletandsampleFINAL.rds",
                            sep = ""), 
        version = 2
       )

rm(sc_data);gc()

