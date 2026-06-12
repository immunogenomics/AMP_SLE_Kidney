# To be used with the "monocole3" environment
suppressPackageStartupMessages({
    library(slingshot)
    library(RColorBrewer)
    library(ggplot2)
    library(dplyr)
    library(igraph)
    library(ggraph)
    library(wesanderson)
})

dPCs <- readRDS("B_cells_diffusionMap_eigenvectors.Rds")
meta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/B_combinedMeta.csv")
meta <- meta[startsWith(meta$annotation, "B") | startsWith(meta$annotation, "bl-B"), ]

print("Slingin'")
sds <- slingshot(dPCs, clusterLabels = meta$annotation, start.clus = "bl−B0. CXCR5high Naive", approx_points = 150)
print("Finished slangin")
saveRDS(sds, "B_destiny_slingshot.Rds")