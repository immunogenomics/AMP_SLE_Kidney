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

umap_object <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/",
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
                            "uppernUMI.rds",
                            sep = ""))


resolution_list <- c(0.1, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0)
ids_ref_cca <- Reduce(cbind, mclapply(resolution_list, function(res_use) {
    Seurat:::RunModularityClustering(SNN = umap_object$fgraph, modularity = 1, 
        resolution = res_use, algorithm = 1, n.start = 20, 
        n.iter = 20, random.seed = 100, print.output = FALSE, 
        temp.file.location = NULL, edge.file.name = NULL)
}, mc.preschedule = FALSE, mc.cores = min(20, length(resolution_list))))
ids_ref_cca %<>% data.frame()
colnames(ids_ref_cca) <- sprintf("res_%.2f", resolution_list)

meta$test_res.0.1 <- ids_ref_cca$res_0.10
meta$test_res.0.1 <- factor(meta$test_res.0.1, levels = sort(unique(meta$test_res.0.1)))
meta$test_res.0.4 <- ids_ref_cca$res_0.40
meta$test_res.0.4 <- factor(meta$test_res.0.4, levels = sort(unique(meta$test_res.0.4)))
meta$test_res.0.6 <- ids_ref_cca$res_0.60
meta$test_res.0.6 <- factor(meta$test_res.0.6, levels = sort(unique(meta$test_res.0.6)))
meta$test_res.0.8 <- ids_ref_cca$res_0.80
meta$test_res.0.8 <- factor(meta$test_res.0.8, levels = sort(unique(meta$test_res.0.8)))
meta$test_res.1.2 <- ids_ref_cca$res_1.20
meta$test_res.1.2 <- factor(meta$test_res.1.2, levels = sort(unique(meta$test_res.1.2)))
meta$test_res.1.6 <- ids_ref_cca$res_1.60
meta$test_res.1.6 <- factor(meta$test_res.1.6, levels = sort(unique(meta$test_res.1.6)))
meta$test_res.2.0 <- ids_ref_cca$res_2.00
meta$test_res.2.0 <- factor(meta$test_res.2.0, levels = sort(unique(meta$test_res.2.0)))

## Save everything
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
                            "uppernUMI.rds",
                            sep = ""))

rm(meta)
rm(umap_object)
rm(ids_ref_cca)
gc()


