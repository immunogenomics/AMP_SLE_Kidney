# Use with the SLE environment
library(harmony)
library(dplyr)

set.seed(0)

args <- commandArgs(trailingOnly=TRUE)
name <- as.character(args[1])
theta <- as.numeric(args[2])

# Read in data
combined_pcs <- readRDS(paste0("PCs/", name, "_PCs.Rds"))
combined_meta <- read.csv(paste0("PCs/", name, "_combinedMeta.csv"), header = TRUE)

harmony <- RunHarmony(combined_pcs$embeddings, combined_meta, 
    c("sample", "batch", "origin"), theta = c(0, 0, theta),  
    max_iter=50, plot_convergence = TRUE)

saveRDS(harmony, paste0(name, "/harmony_", theta, "/", name, "_combined_hPCs.Rds"))