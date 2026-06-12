library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggbeeswarm)

args = commandArgs(trailingOnly=TRUE)

name <- as.character(args[1])
sampleType <- as.character(args[2])
cellTypes <- readRDS(args[4])

ncorrs <- read.csv(as.character(args[3]), header = FALSE)
fdrs <- read.csv(sub("_ncorr.csv", "_fdrs.csv", args[3]), header = FALSE)

threshold <- min(fdrs[fdrs$V2 < 0.1, ]$V1)

if (sampleType == "tissue"){
    folderName = "cna_new"
    umap <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/', folderName, '/', name, '/activity/sc_umap.csv'), header = TRUE)
    colnames(umap) <- c("hUMAP1", "hUMAP2")
    meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/', folderName, '/', name, '/activity/sc_meta.csv'), header = TRUE)
    umap$corr <- ncorrs$V1
    meta <- cbind(meta, umap)
    plotTable <- meta[, c("cell", "final_annotation", "hUMAP1", "hUMAP2", "corr")]
} else {
    folderName = "pbmc"
    umap <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/', folderName, '/', name, '/activity/umap.csv'), header = TRUE)
    meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/', folderName, '/', name, '/activity/meta.csv'), header = TRUE)
    umap$corr <- ncorrs$V1
    meta <- cbind(meta, umap)
    meta <- merge(meta, cellTypes[cellTypes$Cell %in% meta$Cell,c("Cell", "annotation")], by = "Cell")
    plotTable <- meta[, c("Cell", "annotation", "UMAP_1", "UMAP_2", "corr")]
}

colnames(plotTable) <- c("barcode", "cellType", "hUMAP1", "hUMAP2", "corr")

title <- tools::file_path_sans_ext(basename(as.character(args[3])))

tempTable <- plotTable

tempTable$corr[abs(tempTable$corr) < threshold] <- NA
p <- ggplot(
      data = tempTable %>% arrange(desc(is.na(corr)), abs(corr)), 
      aes(x = hUMAP1, y = hUMAP2)) + 
    geom_point(mapping = aes(color = corr), alpha = 5, size = 1) + 
      scale_color_gradient2(low = 'blue', mid = 'lightgray', high = 'red', na.value = "lightgray") + 
      labs(x="UMAP1", y="UMAP2") +
      theme_bw(base_size = 15) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(color="black", hjust = 0.5, face = "bold", size=20)
    ) + ggtitle(title, subtitle = paste0("corr > ", threshold, "; FDR < 0.10"))

pdf(paste0("cna_plots/", name, "_", sampleType, "/", title, ".pdf"), width = 8, height = 8)
    print(p)
dev.off()

myPalette <- colorRampPalette(c(RColorBrewer::brewer.pal(11, "RdBu")[11:7], "#DCDCDC", RColorBrewer::brewer.pal(11, "RdBu")[5:1]))

p <- ggplot(plotTable, 
       aes(x = reorder(cellType, corr), y = corr)) +
    ggbeeswarm::geom_quasirandom(aes(color = corr), width = 0.5, size = 0.75) +
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