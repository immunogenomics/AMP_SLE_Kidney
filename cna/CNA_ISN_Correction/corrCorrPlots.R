# Correlation Plots of Correlation
library(ggplot2)
library(dplyr)
library(pals)
library(ggthemes)
library(ggrepel)

args = commandArgs(trailingOnly=TRUE)

name <- as.character(args[1])
sampleType <- as.character(args[2])
cellTypes <- readRDS(args[4])

activity <- read.csv(as.character(args[3]), header = FALSE)
activity_fdrs <- read.csv(sub("_ncorr.csv", "_fdrs.csv", args[3]), header = FALSE)

chronicity <- read.csv(paste0("cna_results/", name, "_", sampleType, "/Final_Chronicity_justChecking_ncorr.csv"), header = FALSE)
chronicity_fdrs <- read.csv(paste0("cna_results/", name, "_", sampleType, "/Final_Chronicity_justChecking_fdrs.csv"), header = FALSE)

activityThreshold <- min(activity_fdrs[activity_fdrs$V2 < 0.1, ]$V1)
chronicityThreshold <- min(chronicity_fdrs[chronicity_fdrs$V2 < 0.1, ]$V1)

title <- tools::file_path_sans_ext(basename(as.character(args[3])))

if ( sampleType == "pbmc") {
    folderName = "pbmc"
    colorScale <- as.vector(ggthemes_data$tableau$'color-palettes'$regular$`Classic 20`[1])$value
    meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/', folderName, '/', name, '/activity/meta.csv'), header = TRUE)
    meta$chronicityCorrs <- chronicity$V1
    meta$activityCorrs <- activity$V1
    meta <- merge(meta, cellTypes[cellTypes$Cell %in% meta$Cell, c("Cell", "annotation")], by = "Cell")
    plotTable <- data.frame(chronicityCorrs = meta$chronicityCorrs, 
        activityCorrs = meta$activityCorrs, 
        cellTypes = meta$annotation)

} else if ( sampleType == "tissue" ) {
    folderName = "cna_new"
    colorScale <- as.vector(rev(polychrome(26)))
    meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/', folderName, '/', name, '/activity/sc_meta.csv'), header = TRUE)
    plotTable <- data.frame(chronicityCorrs = chronicity$V1, 
        activityCorrs = activity$V1, 
        cellTypes = meta$final_annotation)
}

p <- ggplot(plotTable, aes(x = chronicityCorrs, y = activityCorrs, color = cellTypes, fill = cellTypes)) + 
    geom_hline(yintercept = 0, color = "black", size = 1) +
    geom_vline(xintercept = 0, color = "black", size = 1) +
    ggrastr::rasterise(geom_point(size = 0.01)) + 
    theme_classic() +
    geom_vline(xintercept = chronicityThreshold, linetype = "dashed", color = "black", size = 0.8) +
    geom_vline(xintercept = -1 * chronicityThreshold, linetype = "dashed", color = "black", size = 0.8) +
    geom_hline(yintercept = activityThreshold, linetype = "dashed", color = "black", size = 0.8) +
    geom_hline(yintercept = -1 * activityThreshold, linetype = "dashed", color = "black", size = 0.8) +
    scale_x_continuous(limits = range(plotTable$chronicityCorrs)) + 
    scale_y_continuous(limits = range(plotTable$activityCorrs)) + 
    scale_colour_manual(values = colorScale) + 
    labs(x = "Chronicity Correlation", y = "Activity Correlation") +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 15))) +
    ggtitle(title)

legend <- cowplot::get_legend(p)
p <- p + theme(legend.position = "none")
legend <- cowplot::plot_grid(legend)

cowplot::save_plot(paste0("cna_plots/", name, "_", sampleType, "/correlations/correlation_", title, "_legend.pdf"),
       legend,
       base_height = 7,
       base_width = 14)

dir.create(paste0("cna_plots/", name, "_", sampleType, "/correlations"))
pdf(paste0("cna_plots/", name, "_", sampleType, "/correlations/correlation_", title, ".pdf"), width = 10, height = 10)
    print(p)
dev.off()

centroids <- plotTable%>%
  group_by(cellTypes) %>%
  summarise(
    centroid_x = mean(chronicityCorrs),
    centroid_y = mean(activityCorrs),
    .groups = 'drop'
  ) %>%
  mutate(cellTypes = sub("\\. .*", "", cellTypes))

p <- ggplot(plotTable, aes(x = chronicityCorrs, y = activityCorrs, color = cellTypes, fill = cellTypes)) + 
    geom_hline(yintercept = 0, color = "black", size = 1) +
    geom_vline(xintercept = 0, color = "black", size = 1)  + 
    theme_classic() +
    geom_vline(xintercept = chronicityThreshold, linetype = "dashed", color = "black") +
    geom_vline(xintercept = -1 * chronicityThreshold, linetype = "dashed", color = "black") +
    geom_hline(yintercept = activityThreshold, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -1 * activityThreshold, linetype = "dashed", color = "black") +
    scale_x_continuous(limits = range(plotTable$chronicityCorrs)) + 
    scale_y_continuous(limits = range(plotTable$activityCorrs)) + 
    scale_colour_manual(values = colorScale) + 
    scale_fill_manual(values = colorScale) + 
    stat_ellipse(geom = "polygon", level = 0.68, alpha = 0.25) + 
    labs(x = "Chronicity Correlation", y = "Activity Correlation") +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 15))) +
    ggtitle(title) + theme(legend.position = "none") + 
    geom_text(data = centroids, aes(x = centroid_x, y = centroid_y, label = cellTypes), color = "black", inherit.aes = FALSE)

pdf(paste0("cna_plots/", name, "_", sampleType, "/correlations/correlation_", title, "_centroids.pdf"), width = 10, height = 10)
    print(p)
dev.off()

dir.create(paste0("cna_plots/", name, "_", sampleType, "/correlations/", sub("_ncorr", "", title)))

onePlot <- function(i) {
    cellType <- unique(plotTable$cellTypes)[order(unique(plotTable$cellTypes))][i]
    print(cellType)
    p <- ggplot(plotTable[plotTable$cellTypes == cellType, ], aes(x = chronicityCorrs, y = activityCorrs)) + 
        geom_hline(yintercept = 0, color = "black", size = 1) +
        geom_vline(xintercept = 0, color = "black", size = 1) +
        stat_ellipse(geom = "polygon", level = 0.68, alpha = 0.25, color = colorScale[i], fill = colorScale[i]) + 
        geom_point(size = 0.01, color = colorScale[i], alpha = 1) + 
        theme_classic() +
        geom_vline(xintercept = chronicityThreshold, linetype = "dashed", color = "black", size = 0.8) +
        geom_vline(xintercept = -1 * chronicityThreshold, linetype = "dashed", color = "black", size = 0.8) +
        geom_hline(yintercept = activityThreshold, linetype = "dashed", color = "black", size = 0.8) +
        geom_hline(yintercept = -1 * activityThreshold, linetype = "dashed", color = "black", size = 0.8) + 
        labs(x = "Chronicity Correlation", y = "Activity Correlation") +
        scale_x_continuous(limits = range(plotTable$chronicityCorrs)) + 
        scale_y_continuous(limits = range(plotTable$activityCorrs)) + 
        guides(color = guide_legend(ncol = 2, override.aes = list(size = 15))) +
        ggtitle(title, subtitle = cellType)
    pdf(paste0("cna_plots/", name, "_", sampleType, "/correlations/", sub("_ncorr", "", title), "/", i, ".pdf"), width = 10, height = 10)
        print(p)
    dev.off()
}
lapply(c(1:length(unique(plotTable$cellTypes))), onePlot)