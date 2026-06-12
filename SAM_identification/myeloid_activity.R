suppressPackageStartupMessages({
    library(dplyr)
    library(scales)
    library(pals)
    library(ggplot2)
    library(cowplot)
    library(tidyverse)
    library(extrafont)
    library(ggrastr)
})

labelfontsize = 20
tickfontsize = 16

meta <- read.csv(
    "/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/Myeloid_combinedMeta.csv",
    header = TRUE
)

umap <- readRDS("../harmonyIntegration_finalFinal/Myeloid/harmony_0/Myeloid_umap.Rds")

counters <- read.csv(paste0("../harmonyIntegration_finalFinal/Myeloid/harmony_0/pbmcVsTissue.csv"))

umap_table <- data.frame(barcode = row.names(data.frame(umap$embedding)), 
    UMAP_1 = data.frame(umap$embedding[, 1]), 
    UMAP_2 = data.frame(umap$embedding[, 2]),
    tissueScore = counters$tissue / (counters$tissue + counters$PBMC * (length(meta$origin[meta$origin == "Tissue"]) / 
        length(meta$origin[meta$origin == "PBMC"]))),
    annotation = meta$annotation,
    origin = meta$origin
)

smallMeta <- read.csv('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/myeloid/chronicity/sc_meta.csv')
ncorr <- read.csv(paste0("../CNA_activityRerun/cna_results/myeloid/Final_Activity_SiteFBChron_ncorr.csv"), header = FALSE)
fdrs <- read.csv(paste0("../CNA_activityRerun/cna_results/myeloid/Final_Activity_SiteFBChron_fdrs.csv"), header = FALSE)

colnames(fdrs) <- c('threshold', 'fdr', 'ncells')
smallMeta$ncorr <- ncorr$V1

fullMeta <- merge(smallMeta, umap_table, by.x = "cell", by.y = "barcode")

fdr <- fdrs %>% filter(fdr < 0.1) %>% 
            mutate(fdr = round(fdr, 4)) %>% 
            filter(fdr == max(fdr)) %>% 
            pull(threshold)

pos_fdr_thresh <- fdr
neg_fdr_thresh <- -1 * fdr

vmax = quantile(abs(fullMeta$ncorr), .95)[[1]]
if (vmax < pos_fdr_thresh) {
  vmax <- pos_fdr_thresh + 0.03
}
vmin = -vmax
custom_palette <- c("blue", "lightgray", "lightgray","lightgray",  "red") 
break_points <- c(vmin, neg_fdr_thresh, 0, pos_fdr_thresh, vmax)  

breaks <- seq(round(-max(abs(fullMeta$ncorr)), 2), round(max(abs(fullMeta$ncorr)), 2), round(max(abs(fullMeta$ncorr)), 2))

sig_cmap <- scale_color_gradientn(
  colours = custom_palette,
  values = rescale(break_points, to = c(0, 1)),  # Normalize thresholds
  limits = c(vmin, vmax),  # Ensure full scale
  na.value = "pink",
  guide = guide_colorbar(direction = "horizontal"),  # Horizontal legend
    breaks = c(-floor(abs(vmin)*100)/100, floor(vmax*100)/100),
oob = scales::squish
)

plotMeta <- fullMeta %>% select(cell, ncorr, annotation)
plotMeta <- plotMeta %>% 
  mutate(type = case_when(
    grepl("^M7|^M9|^M11|^M5", annotation) & grepl("Macrophage", annotation) ~ "SAM",
    grepl("Macrophage", annotation) ~ "non-SAM",
    TRUE ~ "Monocyte"
  ))

vlnPlot <- ggplot(plotMeta, aes(x = type, y = ncorr)) +
                    ggrastr::rasterise(ggbeeswarm::geom_quasirandom(aes(color = ncorr), width = 0.25, size = 0.5), dpi = 400) +
                    geom_boxplot(width = 0.25, alpha = 0) +
                    geom_hline(yintercept = neg_fdr_thresh, linetype = "dashed", color = "darkgrey") +
                    geom_hline(yintercept = pos_fdr_thresh, linetype = "dashed", color = "darkgrey") +
                    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
                    sig_cmap + 
                    scale_y_continuous(breaks = breaks) +
                    labs( y= "Activity CNA Correlation", x = "", title = "") +
                    theme_classic(base_size = tickfontsize) +
                    theme(
                        text=element_text(family="Arial Narrow"),
                        legend.position = "none",
                        panel.grid = element_blank(),
                        axis.text.x = element_text(color = "black", size = labelfontsize, hjust=0.5),
                        axis.text.y = element_text(color = "black", size = labelfontsize, hjust=1),
                        axis.title = element_text(size=labelfontsize, hjust = 0.5),
                        axis.line.x.bottom = element_line(color = 'black'),
                        axis.line.y.left   = element_line(color = 'black'),
                        plot.margin = margin(10, 40, 10, 10)
                    )

pdf("SAMs_Activity_vlns.pdf", width = 5, height = 5)
    print(vlnPlot)
dev.off()