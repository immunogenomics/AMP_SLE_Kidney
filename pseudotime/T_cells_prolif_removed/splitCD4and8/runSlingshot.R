# To be used with the "monocole3" environment
library(slingshot)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)
library(wesanderson)

# CD4 Cells
dPCs <- readRDS("pt_outputs/CD4_diffusionMap_eigenvectors.Rds")
meta <- read.csv("CD4_inputs/CD4_T_cells_meta.csv")

sds <- slingshot(dPCs, clusterLabels = meta$annotation, start.clus = "bl-T10. CD4+ Central Memory/Naive", approx_points = 150)
saveRDS(sds, "pt_outputs/CD4_destiny_slingshot.Rds")

# CD8 Cells
dPCs <- readRDS("pt_outputs/CD8_diffusionMap_eigenvectors.Rds")
meta <- read.csv("CD8_inputs/CD8_T_cells_meta.csv")

sds <- slingshot(dPCs, clusterLabels = meta$annotation, start.clus = "bl-T8. CD8+ Central Memory/Naive", approx_points = 150)
saveRDS(sds, "pt_outputs/CD8_destiny_slingshot.Rds")