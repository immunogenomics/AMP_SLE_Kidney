library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggbeeswarm)

args = commandArgs(trailingOnly=TRUE)

name <- "myeloid"
sampleType <- "tissue"
cellTypes <- '/data/srlab/ssg34/SLE_kidney_v2/data/qcd/Myeloid_clusterQCd_meta_harmonizedPCUMAPCellStateClusters_10042022.rds'

ncorrs <- read.csv("cna_results/myeloid_tissue/Final_Activity_ChronicityBiopAndResponse_ncorr.csv", header = FALSE)
fdrs <- read.csv("cna_results/myeloid_tissue/Final_Activity_ChronicityBiopAndResponse_fdrs.csv", header = FALSE)

new_ncorrs <- read.csv("cna_results/myeloid_tissue/Final_Activity_withISN_V_ncorr.csv", header = FALSE)
new_fdrs <- read.csv("cna_results/myeloid_tissue/Final_Activity_withISN_V_fdrs.csv", header = FALSE)

threshold <- min(new_fdrs[new_fdrs$V2 < 0.1, ]$V1)
oldThreshold <- min(fdrs[fdrs$V2 < 0.1, ]$V1)

folderName = "cna_new"
umap <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/', folderName, '/', name, '/activity/sc_umap.csv'), header = TRUE)
colnames(umap) <- c("hUMAP1", "hUMAP2")
meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/', folderName, '/', name, '/activity/sc_meta.csv'), header = TRUE)
umap$oldCorr <- ncorrs$V1
umap$newCorr <- new_ncorrs$V1
meta <- cbind(meta, umap)

plotTable <- meta[, c("cell", "final_annotation", "hUMAP1", "hUMAP2", "oldCorr", "newCorr")]

colnames(plotTable) <- c("barcode", "cellType", "hUMAP1", "hUMAP2", "oldCorr", "newCorr")

title <- "beforeAndAfter"

myPalette <- colorRampPalette(c(RColorBrewer::brewer.pal(11, "RdBu")[11:7], "#DCDCDC", RColorBrewer::brewer.pal(11, "RdBu")[5:1]))

p <- ggplot(plotTable, 
       aes(x = reorder(cellType, oldCorr), y = newCorr)) +
    ggbeeswarm::geom_quasirandom(aes(color = newCorr), width = 0.5, size = 0.75) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "darkgrey") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
    geom_hline(yintercept = -1 * threshold, linetype = "dashed", color = "darkgrey") +
    scale_colour_gradientn(colours = myPalette(20)) +
   labs(title = title, x= "", y = "Correlation") +
    theme_classic() +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(color="black", size=22, face = "bold"),
        axis.text.x = element_text(color = "black", size = 15, angle = 90, hjust=1, face = "bold"),
        axis.title = element_text(size=20, face = "bold"))

pdf(paste0("cna_plots/", name, "_", sampleType, "/", title, "_cellType.pdf"), width = 8, height = 8)
    print(p)
dev.off()

oldThreshold <- min(fdrs[fdrs$V2 < 0.1, ]$V1)

p <- ggplot(plotTable, 
       aes(x = oldCorr, y = newCorr)) +
    geom_point(aes(color = newCorr)) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "darkgrey") +
    geom_hline(yintercept = -1 * threshold, linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = oldThreshold, linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = -1 * oldThreshold, linetype = "dashed", color = "darkgrey") +
    scale_colour_gradientn(colours = myPalette(20)) +
   labs(title = title, x= "noCorrection", y = "withCorrection") +
    theme_classic() +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(color="black", size=22, face = "bold"))

pdf(paste0("cna_plots/", name, "_", sampleType, "/", title, "_scatter.pdf"), width = 8, height = 8)
    print(p)
dev.off()