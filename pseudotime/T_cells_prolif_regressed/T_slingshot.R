# To be used with the "monocole3" environment
library(slingshot)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)
library(wesanderson)

dPCs <- readRDS("T_cells_diffusionMap_eigenvectors.Rds")
meta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/TNK_combinedMeta.csv")

tnk_meta <- meta[startsWith(meta$annotation, "T") | startsWith(meta$annotation, "bl-T"), ]

sds <- slingshot(dPCs, clusterLabels = tnk_meta$annotation, start.clus = "bl-T12. TRGC1+ Gamma/Delta", approx_points = 150)
saveRDS(sds, "T_destiny_slingshot.Rds")