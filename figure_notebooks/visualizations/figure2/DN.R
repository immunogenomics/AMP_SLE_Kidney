library(singlecellmethods)
library(Seurat)
library(parallel)
library(ggplot2)
library(dplyr)
library(stringr)

meta <- readRDS("../finalObjects/dn_meta.Rds")

# Taken from some pals color palette that I forgot the name of...
palette <- c("#782AB6", "#C075A6", "#F7E1A0", "#FC1CBF", "#FA0087", "#1CBE4F", "#FBE426",
    "#B10DA1", "#85660D", "#1C8356", "#C4451C", "#325A9B", "#F8A19F", "#AA0DFE",
    "#DEA0FD", "#2ED9FF", "#90AD1C", "#1CFFCE", "#B00068", "#FEAF16", "#3283FE",
    "#16FF32")

meta$new_cluster_number <- substr(meta$final_annotation, 3, 3)

cluster_center <- meta %>%
                  group_by(new_cluster_number, final_annotation) %>%
                  summarise_at(vars(huwotUMAP1, huwotUMAP2), funs(median(., na.rm=TRUE)))

cluster_center <- as.data.frame(cluster_center)
cluster_center$final_annotation <- as.factor(cluster_center$final_annotation)

# New clustering
p <- ggplot(meta, aes(x = huwotUMAP1, y = huwotUMAP2, color = final_annotation)) +
        ggrastr::rasterise(geom_point(size = 0.05), dpi = 400) +
        ggrepel::geom_label_repel(
          data = cluster_center,
          aes(x = huwotUMAP1, y = huwotUMAP2, 
              label = new_cluster_number, color = as.factor(final_annotation)),
          size = 7,  fontface = "bold",
          box.padding = unit(0.5, "lines"),
          point.padding = unit(0.01, "lines"),
          show.legend = FALSE) +
        scale_color_manual(values = palette) +
        theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), plot.title = element_text(size = 15)) +
        theme_bw() +
        theme(
        legend.position = "right",
              axis.title = element_text(hjust = 0.75, 
                                        size = 10, 
                                        face = "bold"), 
              axis.text = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              legend.text = element_text(size = 20),
              legend.title = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) + 
        labs(title = "", x = 'UMAP1', y = 'UMAP2') +
        guides(color = guide_legend(ncol = 1, override.aes = list(size = 15)))

legend <- cowplot::get_legend(p)

p <- p + theme(legend.position = "none")

cowplot::save_plot("DN_Tissue_legend.pdf",
       legend,
       base_height = 10,
       base_width = 10)

pdf("DN_Tissue.pdf", width = 5, height = 5)
    print(p)
dev.off()