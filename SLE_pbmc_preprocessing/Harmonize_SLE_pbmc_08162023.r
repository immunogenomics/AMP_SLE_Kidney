#@SCRIPT: Harmonizing of SLE pbmc Data
#@DATE: 08/15/2023
#@AUTHOR: SID GURAJALA

library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(stringr)
library(pheatmap)
library(ggrepel)
library(presto)
library(pals)
library(Seurat)
library(harmony)
library(singlecellmethods)
library(lisi)
source("/data/srlab/ik936/Foxxy/utils/utils.R")
source("/data/srlab/anathan/scripts/scseq_utils.R")
library(parallel)
set.seed(0)


all_obj <- readRDS('/data/srlab/ssg34/SLE_pbmc_analysis/data/gene_expression/demux_cell_qcd_gex_seurat_08152023.rds')

all_obj <- all_obj %>% 
                    RunHarmony(c("Batch", "Unified_Visit"), plot_convergence = TRUE)


all_obj <- all_obj %>% 
                RunUMAP(reduction = "harmony", dims = 1:20) %>% 
                FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
                FindClusters(all_obj, resolution = c(0, 0.25, 0.5, 0.75, 1.0))


all_obj <- DietSeurat(all_obj, counts = TRUE, data = TRUE, 
                        scale.data = FALSE, dimreducs = "pca")


downsample_cells <- Cells(all_obj)[sample(length(Cells(all_obj)), 50000)]

downsample_obj <- subset(all_obj, cells = downsample_cells)



saveRDS(downsample_obj, '/data/srlab/ssg34/SLE_pbmc_analysis/data/gene_expression/demux_cell_qcd_gex_seurat_sc_analysis_downsampled_08182023.rds')

saveRDS(all_obj@meta.data, '/data/srlab/ssg34/SLE_pbmc_analysis/data/metadata/demux_cell_qcd_gex_seurat_sc_analysis_metadata_08182023.rds')

saveRDS(Embeddings(all_obj, reduction = "harmony"), '/data/srlab/ssg34/SLE_pbmc_analysis/data/gene_expression/demux_cell_qcd_gex_seurat_sc_analysis_harmony_08182023.rds')

saveRDS(all_obj, '/data/srlab/ssg34/SLE_pbmc_analysis/data/gene_expression/demux_cell_qcd_gex_seurat_sc_analysis_08182023.rds')