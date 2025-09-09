d#@SCRIPT: SEURAT PREPROCESSING OF SLE
#@DATE: 08/15/2023
#@AUTHOR: SID GURAJALA

library(symphony)
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


working_dir = "/data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/"

flow_cell_list = c('220910_SL-NVA_0921_BHWM3JDSX3',
                  '220910_SL-NVA_0920_AHWMHWDSX3', 
                  '220910_SL-NVU_0582_AHWMKWDSX3',
                  '220910_SL-NVR_0633_AHWNYKDSX3',
                  '220910_SL-NVR_0632_BHWFFHDSX3',
                  '220910_SL-NVU_0581_BHWMNGDSX3')


format_seurat_channel <- function(channel, flow_dir) { 
    channel_10X_dir <- paste0(flow_dir, '/cellranger-6.1.1/GRCh38/', channel, '/outs/filtered_feature_bc_matrix')
    channel_10X_obj <- Read10X(channel_10X_dir)
    channel_seurat_obj <- CreateSeuratObject(counts = channel_10X_obj$`Gene Expression`)
    channel_seurat_obj@meta.data$BARCODE <- rownames(channel_seurat_obj@meta.data)
    channel_seurat_obj <- RenameCells(channel_seurat_obj, 
                                      new.names = paste0(channel, '_', Cells(channel_seurat_obj)))
    channel_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(channel_seurat_obj, pattern = "^MT-")
    channel_seurat_obj <- subset(channel_seurat_obj, subset = nFeature_RNA > 500  & percent.mt < 20)
    channel_demux <- read.table(paste0(flow_dir, '/genotyping/', channel, ".best"),
                                 header = TRUE) %>% 
                                filter(str_detect(BEST, "SNG"))
    cells.keep <- rownames(channel_seurat_obj@meta.data)[channel_seurat_obj@meta.data$BARCODE %in% channel_demux$BARCODE]
    channel_seurat_obj <- subset(channel_seurat_obj, cells = cells.keep)
    row.names <- rownames(channel_seurat_obj@meta.data)
    channel_seurat_obj@meta.data <- left_join(channel_seurat_obj@meta.data, channel_demux %>% select(BARCODE, SNG.1ST))
    rownames(channel_seurat_obj@meta.data) <- row.names
    channel_seurat_obj@meta.data$Channel <- channel
return(channel_seurat_obj)
}

format_seurat_flow <- function(flow_cell) {
    flow_dir <- paste0(working_dir, flow_cell) 
    channels <- read.table(paste0(flow_dir, '/lsf_params_count_feature'),  header = F)$V2
    flow_seurat_obj <- lapply(channels, format_seurat_channel, flow_dir)
    combined_seurat_obj <- flow_seurat_obj[[1]]
    for (i in 2:length(flow_seurat_obj)) {
        combined_seurat_obj <- merge(x = combined_seurat_obj, 
                                     y = flow_seurat_obj[[i]])

}
    rm(flow_seurat_obj)
    return(combined_seurat_obj)
}


all_flow_cells_obj <- lapply(flow_cell_list, format_seurat_flow)
saveRDS(all_flow_cells_obj, '/data/srlab/ssg34/SLE_pbmc_analysis/data/gene_expression/demux_cell_qcd_gex_seurat_list_ 08152023.rds')


all_obj <- all_flow_cells_obj[[1]]
for (i in 2:length(all_flow_cells_obj)) {
        all_obj <- merge(x = all_obj, 
                         y = all_flow_cells_obj[[i]])
}

saveRDS(all_obj, '/data/srlab/ssg34/SLE_pbmc_analysis/data/gene_expression/demux_cell_qcd_gex_seurat_sc_analysis_08152023.rds')