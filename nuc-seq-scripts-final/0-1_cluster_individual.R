#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")

#-------------------------------------------------------

filepath1 <- "/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/"
filepath2 <- "/outs/raw_feature_bc_matrix/"
samples <- read.delim("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/list_of_nuc_samples.txt",
                            sep = "\n", header = F)
samples <- as.character(samples[,1])

samples <- c("AMPSLEkid_nuc_0137")
#-------------------------------------------------------

## FILTERING PARAMETERS
### NGENE : 500
### MT    : NONE

nGenepar <- 500

#-------------------------------------------------------


for (i in 1:length(samples)) {

## Read in filtered_data
    data <- readRDS(paste("../nuc-seq-output/filtered_data/",
                                 samples[i],
                                 "_data_filtered_",
                                 nGenepar,
                                 "nGene_clusterremoval.rds",
                                 sep = ""))
## Read in metadata
    meta <- readRDS(paste("../nuc-seq-output/full_metadata/",
                            samples[i],
                            "_meta_filtered_",
                            nGenepar,
                            "nGene_clusterremoval.rds",
                            sep = ""))


sc_data <- CreateSeuratObject(data)
rm(data);gc()
sc_data <- NormalizeData(sc_data)
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 2000)
sc_data <- ScaleData(sc_data)
sc_data <- RunPCA(sc_data, features = VariableFeatures(object = sc_data))

saveRDS(sc_data@assays$RNA@data, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            samples[i],
                            "_normalized_data_filtered_",
                            nGenepar,
                            "nGene_clusterremoval.rds",
                            sep = ""),
        version = 2
       )

saveRDS(sc_data@reductions$pca@cell.embeddings, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            samples[i],
                            "_pca_cell_embeddings_filtered_",
                            nGenepar,
                            "nGene_clusterremoval.rds",
                            sep = ""),
        version = 2
       )

saveRDS(sc_data@assays$RNA@var.features, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            samples[i],
                            "_variable_genes_filtered_",
                            nGenepar,
                            "nGene_clusterremoval.rds",
                            sep = ""),
        version = 2
       )

saveRDS(sc_data@reductions$pca@feature.loadings, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            samples[i],
                            "_feature_loadings_filtered_",
                            nGenepar,
                            "nGene_clusterremoval.rds",
                            sep = ""),
        version = 2
       )

saveRDS(sc_data, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            samples[i],
                            "_seurat_object_filtered",
                            nGenepar,
                            "nGene_clusterremoval.rds",
                            sep = ""), 
        version = 2
       )

rm(sc_data);gc()


pca_res <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            samples[i],
                            "_pca_cell_embeddings_filtered_",
                            nGenepar,
                            "nGene_clusterremoval.rds",
                            sep = "")
                  )


umap_res <- umap$UMAP(n_neighbors = 30L, metric = "euclidean", min_dist = .3)$fit_transform(pca_res[,1:20])
# umap_res <- umap$UMAP(n_neighbors = 30L, metric = "correlation", min_dist = .3)$fit_transform(harmony)
meta$pUMAP1 <- umap_res[, 1]
meta$pUMAP2 <- umap_res[, 2]

saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_metadata/",
                            samples[i],
                            "_meta_filtered_",
                            nGenepar,
                            "nGene_clusterremoval.rds",
                            sep = "")
       )


## Build SNNGraph and calculate louvain clustering at various resolutions
SNNGraph <- BuildSNNSeurat(pca_res[,1:20], k.param = 30, nn.eps = 0)

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

meta$pres.0.1 <- ids_ref_cca$res_0.10
meta$pres.0.1 <- factor(meta$pres.0.1, levels = sort(unique(meta$pres.0.1)))
meta$pres.0.4 <- ids_ref_cca$res_0.40
meta$pres.0.4 <- factor(meta$pres.0.4, levels = sort(unique(meta$pres.0.4)))
meta$pres.0.6 <- ids_ref_cca$res_0.60
meta$pres.0.6 <- factor(meta$pres.0.6, levels = sort(unique(meta$pres.0.6)))
meta$pres.0.8 <- ids_ref_cca$res_0.80
meta$pres.0.8 <- factor(meta$pres.0.8, levels = sort(unique(meta$pres.0.8)))
meta$pres.1.2 <- ids_ref_cca$res_1.20
meta$pres.1.2 <- factor(meta$pres.1.2, levels = sort(unique(meta$pres.1.2)))
meta$pres.1.6 <- ids_ref_cca$res_1.60
meta$pres.1.6 <- factor(meta$pres.1.6, levels = sort(unique(meta$pres.1.6)))
meta$pres.2.0 <- ids_ref_cca$res_2.00
meta$pres.2.0 <- factor(meta$pres.2.0, levels = sort(unique(meta$pres.2.0)))

## Save everything
saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_metadata/",
                            samples[i],
                            "_meta_filtered_",
                            nGenepar,
                            "nGene_clusterremoval.rds",
                            sep = "")
       )

rm(pca_res)
rm(meta)
rm(ids_ref_cca)
rm(umap_res)
gc()

}


