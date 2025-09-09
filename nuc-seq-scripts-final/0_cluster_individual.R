#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")

#-------------------------------------------------------

filepath1 <- "/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/"
filepath2 <- "/outs/raw_feature_bc_matrix/"
samples <- read.delim("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/list_of_nuc_samples.txt",
                            sep = "\n", header = F)
samples <- as.character(samples[,1])

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


for (i in 1:length(samples)) {

## Read in filtered_data
    data <- readRDS(paste("../nuc-seq-output/filtered_data/",
                                 samples[i],
                                 "_data_filtered_",
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
## Read in metadata
    meta <- readRDS(paste("../nuc-seq-output/filtered_metadata/",
                            samples[i],
                            "_meta_filtered_",
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


sc_data <- CreateSeuratObject(data)
rm(data);gc()
sc_data <- NormalizeData(sc_data)
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 2000)
sc_data <- ScaleData(sc_data)
sc_data <- RunPCA(sc_data, features = VariableFeatures(object = sc_data), npcs = 40)

saveRDS(sc_data@assays$RNA@data, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            samples[i],
                            "_normalized_data_filtered_",
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
                            sep = ""),
        version = 2
       )

saveRDS(sc_data@reductions$pca@cell.embeddings, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            samples[i],
                            "_pca_cell_embeddings_filtered_",
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
                            sep = ""),
        version = 2
       )

saveRDS(sc_data@assays$RNA@var.features, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            samples[i],
                            "_variable_genes_filtered_",
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
                            sep = ""),
        version = 2
       )

saveRDS(sc_data@reductions$pca@feature.loadings, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            samples[i],
                            "_feature_loadings_filtered_",
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
                            sep = ""),
        version = 2
       )

saveRDS(sc_data, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            samples[i],
                            "_seurat_object_filtered",
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
                            sep = ""), 
        version = 2
       )

rm(sc_data);gc()


pca_res <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            samples[i],
                            "_pca_cell_embeddings_filtered_",
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
                            sep = "")
                  )



umap_object <- uwot::umap(
    X = pca_res[,1:20],
    ret_extra = c('nn', 'fgraph', 'model')
) 
meta$uwotUMAP1 <- umap_object$embedding[, 1]
meta$uwotUMAP2 <- umap_object$embedding[, 2]

saveRDS(umap_object, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
                            samples[i],
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
                            "uppernUMI.rds",
                            sep = ""))

saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_metadata/",
                            samples[i],
                            "_meta_filtered_",
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
                            sep = "")
       )



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
saveRDS(meta, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_metadata/",
                            samples[i],
                            "_meta_filtered_",
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
                            sep = "")
       )

rm(pca_res)
rm(meta)
rm(ids_ref_cca)
rm(umap_res)
gc()

}


