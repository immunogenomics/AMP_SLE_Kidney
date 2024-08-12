#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")

data <- readRDS("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_data/data_all_corrected_overlap.rds",
		)


sc_data <- CreateSeuratObject(data)
rm(data);gc()
sc_data <- NormalizeData(sc_data)
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 5000)
sc_data <- ScaleData(sc_data)
sc_data <- RunPCA(sc_data, features = VariableFeatures(object = sc_data))

saveRDS(sc_data@assays$RNA@data, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "normalized_data_all_corrected_overlap.rds",
                            sep = ""),
        version = 2
       )

saveRDS(sc_data@reductions$pca@cell.embeddings, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "pca_cell_embeddings_filtered_all_corrected_overlap.rds",
                            sep = ""),
        version = 2
       )

saveRDS(sc_data@assays$RNA@var.features, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "variable_genes_filtered_all_corrected_overlap.rds",
                            sep = ""),
        version = 2
       )

saveRDS(sc_data@reductions$pca@feature.loadings, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "feature_loadings_filtered_all_corrected_overlap.rds",
                            sep = ""),
        version = 2
       )

saveRDS(sc_data, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "seurat_object_filtered_all_corrected_overlap.rds",
                            sep = ""), 
        version = 2
       )

rm(sc_data);gc()

pca_res <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "pca_cell_embeddings_filtered_all_corrected_overlap.rds",
                            sep = "")
)

meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_all_corrected_overlap.rds",
                            sep = "")
)

umap_object <- uwot::umap(
    X = pca_res[,1:20],
    ret_extra = c('nn', 'fgraph', 'model')
) 
meta$uwotUMAP1 <- umap_object$embedding[, 1]
meta$uwotUMAP2 <- umap_object$embedding[, 2]

saveRDS(umap_object, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "umap_all_filtered_all_corrected_overlap.rds",
                            sep = "")
)

saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_all_corrected_overlap.rds",
                            sep = "")
)

rm(pca_res)
rm(meta)

gc()

umap_object <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "umap_all_filtered_all_corrected_overlap.rds",
                            sep = "")
)

meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_all_corrected_overlap.rds",
                            sep = "")
)

## Build SNNGraph and calculate louvain clustering at various resolutions

resolution_list <- c(0.1, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0)
ids_ref_cca <- Reduce(cbind, mclapply(resolution_list, function(res_use) {
    Seurat:::RunModularityClustering(SNN = umap_object$fgraph, modularity = 1, 
        resolution = res_use, algorithm = 1, n.start = 20, 
        n.iter = 20, random.seed = 100, print.output = FALSE, 
        temp.file.location = NULL, edge.file.name = NULL)
}, mc.preschedule = FALSE, mc.cores = min(20, length(resolution_list))))
ids_ref_cca %<>% data.frame()
print("Removing SNNGraph")
rm(SNNGraph)
gc()
colnames(ids_ref_cca) <- sprintf("res_%.2f", resolution_list)

meta$res.0.1 <- ids_ref_cca$res_0.10
meta$res.0.1 <- factor(meta$res.0.1, levels = sort(unique(meta$res.0.1)))
meta$res.0.4 <- ids_ref_cca$res_0.40
meta$res.0.4 <- factor(meta$res.0.4, levels = sort(unique(meta$res.0.4)))
meta$res.0.6 <- ids_ref_cca$res_0.60
meta$res.0.6 <- factor(meta$res.0.6, levels = sort(unique(meta$res.0.6)))
meta$res.0.8 <- ids_ref_cca$res_0.80
meta$res.0.8 <- factor(meta$res.0.8, levels = sort(unique(meta$res.0.8)))
meta$res.1.2 <- ids_ref_cca$res_1.20
meta$res.1.2 <- factor(meta$res.1.2, levels = sort(unique(meta$res.1.2)))
meta$res.1.6 <- ids_ref_cca$res_1.60
meta$res.1.6 <- factor(meta$res.1.6, levels = sort(unique(meta$res.1.6)))
meta$res.2.0 <- ids_ref_cca$res_2.00
meta$res.2.0 <- factor(meta$res.2.0, levels = sort(unique(meta$res.2.0)))

## Save everything
saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_all_corrected_overlap.rds",
                            sep = "")
)

rm(pca_res)
rm(meta)
rm(ids_ref_cca)
gc()


pca_res <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "pca_cell_embeddings_filtered_all_corrected_overlap.rds",
                            sep = "")
)

meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_all_corrected_overlap.rds",
                            sep = "")
)

## harmony

    harmony <- HarmonyMatrix(pca_res[,1:20], meta, vars_use = c("Sample","processing.batch", "Site"), theta = c(1,1,1),
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
# meta <- cbind(meta, harmony)


saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_all_corrected_overlap.rds",
                            sep = "")
)

saveRDS(harmony, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/harmony_output/",
                            "harmony_filtered_all_corrected_overlap.rds",
			    sep = "")
       )

rm(pca_res)
rm(meta)
rm(harmony)

gc()

meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_all_corrected_overlap.rds",
                            sep = "")
)

harmony <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/harmony_output/",
                            "harmony_filtered_all_corrected_overlap.rds",
			    sep = "")
)

umap_object <- uwot::umap(
    X = harmony,
    ret_extra = c('nn', 'fgraph', 'model')
) 
meta$huwotUMAP1 <- umap_object$embedding[, 1]
meta$huwotUMAP2 <- umap_object$embedding[, 2]

saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_all_corrected_overlap.rds",
                            sep = "")
)

saveRDS(umap_object, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "umap_harmony_all_filtered_all_corrected_overlap.rds",
                            sep = "")
)

rm(harmony)
rm(meta)

gc()

meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_all_corrected_overlap.rds",
                            sep = "")
)

umap_object <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "umap_harmony_all_filtered_all_corrected_overlap.rds",
                            sep = "")
)

## Build SNNGraph and calculate louvain clustering at various resolutions

resolution_list <- c(0.1, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0)
ids_ref_cca <- Reduce(cbind, mclapply(resolution_list, function(res_use) {
    Seurat:::RunModularityClustering(SNN = umap_object$fgraph, modularity = 1, 
        resolution = res_use, algorithm = 1, n.start = 20, 
        n.iter = 20, random.seed = 100, print.output = FALSE, 
        temp.file.location = NULL, edge.file.name = NULL)
}, mc.preschedule = FALSE, mc.cores = min(20, length(resolution_list))))
ids_ref_cca %<>% data.frame()
gc()
colnames(ids_ref_cca) <- sprintf("res_%.2f", resolution_list)

meta$hres.0.1 <- ids_ref_cca$res_0.10
meta$hres.0.1 <- factor(meta$hres.0.1, levels = sort(unique(meta$hres.0.1)))
meta$hres.0.4 <- ids_ref_cca$res_0.40
meta$hres.0.4 <- factor(meta$hres.0.4, levels = sort(unique(meta$hres.0.4)))
meta$hres.0.6 <- ids_ref_cca$res_0.60
meta$hres.0.6 <- factor(meta$hres.0.6, levels = sort(unique(meta$hres.0.6)))
meta$hres.0.8 <- ids_ref_cca$res_0.80
meta$hres.0.8 <- factor(meta$hres.0.8, levels = sort(unique(meta$hres.0.8)))
meta$hres.1.2 <- ids_ref_cca$res_1.20
meta$hres.1.2 <- factor(meta$hres.1.2, levels = sort(unique(meta$hres.1.2)))
meta$hres.1.6 <- ids_ref_cca$res_1.60
meta$hres.1.6 <- factor(meta$hres.1.6, levels = sort(unique(meta$hres.1.6)))
meta$hres.2.0 <- ids_ref_cca$res_2.00
meta$hres.2.0 <- factor(meta$hres.2.0, levels = sort(unique(meta$hres.2.0)))

## Save everything
saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_all_corrected_overlap.rds",
                            sep = "")
)

rm(umap_object)
rm(meta)
rm(ids_ref_cca)
gc()


