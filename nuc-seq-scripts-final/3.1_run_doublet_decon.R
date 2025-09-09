#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")
library(DoubletDecon)


#-------------------------------------------------------

filepath1 <- "/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/"
filepath2 <- "/outs/filtered_feature_bc_matrix/"
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

## Read in data
    data <- readRDS(paste("../nuc-seq-output/filtered_data/",
                                 samples[i],
                                 "_data_filtered_",
                                 nGenepar,
                                 "nGene.rds",
                                 sep = "")
		)

## Read in metadata
    meta <- readRDS(paste("../nuc-seq-output/filtered_metadata/",
                            samples[i],
                            "_meta_filtered_",
                            nGenepar,
                            "nGene.rds",
                            sep = "")
		)


## make seurat object and preprocess
    sc_data <- CreateSeuratObject(data, project = "test", meta.data = meta)
    sc_data <- NormalizeData(sc_data)
    sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(sc_data)
    sc_data <- ScaleData(sc_data, features = all.genes)
    sc_data <- RunPCA(sc_data, features = VariableFeatures(object = sc_data), npcs = 50)

    sc_data <- FindNeighbors(sc_data, dims = 1:20)
    sc_data <- FindClusters(sc_data, resolution = 0.6)
    sc_data <- RunUMAP(sc_data, dims = 1:20)


## run doubletdecon
    deconfile <- Improved_Seurat_Pre_Process(sc_data, num_genes=50,
                                         write_files=F)

    res <- Main_Doublet_Decon(rawDataFile = deconfile$newExpressionFile,
                   groupsFile = deconfile$newGroupsFile,
                   filename = samples[i],
                   location = "/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/doubletdecon_output/",
                   rhop = 1)

## add doublet true/false to meta
    meta$barcode <- as.character(lapply(strsplit(meta$cells, "-"), "[[", 1))
    res$DRS_doublet_table$barcode <- as.character(lapply(strsplit(rownames(res$DRS_doublet_table), "[.]"), "[[", 1))
    all(meta$barcode == res$DRS_doublet_table$barcode)
    ordered_res <- res$DRS_doublet_table[order(match(res$DRS_doublet_table$barcode, meta$barcode)),]
    all(meta$barcode == ordered_res$barcode)
    meta$doublet <- ordered_res$isADoublet
    sc_data@meta.data$doublet <- ordered_res$isADoublet

## save everything
    saveRDS(meta, paste("../nuc-seq-output/filtered_metadata/",
                            samples[i],
                            "_meta_filtered_",
                            nGenepar,
                            "nGene_doubletdecon.rds",
                            sep = "")
                )

    saveRDS(sc_data, paste("../nuc-seq-output/seurat_output/",
                            samples[i],
                            "_seurat_obj_",
                            nGenepar,
                            "nGene_doubletdecon.rds",
                            sep = "")
		)

    saveRDS(res, paste("../nuc-seq-output/doubletdecon_output/",
                            samples[i],
                            "_doubletdecon_out_",
                            nGenepar,
                            "nGene_doubletdecon.rds",
                            sep = "")
		)

rm(data)
rm(meta)
rm(res)
rm(deconfile)
rm(sc_data)
gc()


}



