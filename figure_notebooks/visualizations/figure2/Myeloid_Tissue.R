library(singlecellmethods)
library(Seurat)
library(parallel)
library(ggplot2)
library(dplyr)
library(stringr)

meta <- readRDS("../reorderedClusters/myeloid_tissue_meta_reOrdered.Rds")

# Taken from some pals color palette that I forgot the name of...
palette <- c("#782AB6", "#C075A6", "#F7E1A0", "#FC1CBF", "#FA0087", "#1CBE4F", "#FBE426",
    "#B10DA1", "#85660D", "#1C8356", "#C4451C", "#325A9B", "#F8A19F", "#AA0DFE",
    "#DEA0FD", "#2ED9FF", "#90AD1C", "#1CFFCE", "#B00068", "#FEAF16", "#3283FE",
    "#16FF32")

meta <- meta %>%
  mutate(new_cluster_number = str_extract(final_annotation, "^[^.]+"),  # Extract everything before the period
         new_cluster_number = ifelse(str_starts(new_cluster_number, "DC"), substr(new_cluster_number, 3, nchar(new_cluster_number)), substr(new_cluster_number, 2, nchar(new_cluster_number))))

cluster_center <- meta %>%
                  group_by(new_cluster_number, final_annotation) %>%
                  summarise_at(vars(hUMAP1, hUMAP2), funs(median(., na.rm=TRUE)))

cluster_center <- as.data.frame(cluster_center)
cluster_center$final_annotation <- as.factor(cluster_center$final_annotation)

levels(cluster_center$final_annotation) <- c(
      'M0. CD16+ CXC3CR1+ Monocyte', 
      'M1. CD14+ CD16+ CCL2+ CX3CR1+ Monocyte',
      'M2. CD14+ CCR2+ Monocyte',
      'M3. CCL2+ CCL3+ Monocyte',
      'M4. TPSB2+ MAST cell',
      'M5. GPNMBhigh NUPR1high Macrophage',
      'M6. SELENOPinter ISGhigh Macrophage',
      'M7. SPP1high FABP5high Macrophage',
      'M8. SPP1low FABP5high Macrophage',
      'M9. MERTKhigh FABP5high Macrophage',
      'M10. SELENOPinter LYVE1inter Resident Macrophage',
      'M11. GPMNBhigh NUPR1low Macrophage',
      'M12. SELENOPhigh LYVE1high Resident Macrophage',
      'M13. APOChigh C3high Macrophage',
      'M14. APOClow C3high Macrophage',
      'M15. CENPF+ MKI67+ Proliferating',
      'DC16. CCR7+ LAMP3+ DC2', 
      'DC17. CLEC10Alow cDC2',
      'DC18. CLEC10Ahigh cDC2',
      'DC19. cDC1',
      'DC20. pDC'
)

meta$final_annotation <- as.factor(meta$final_annotation)
levels(meta$final_annotation) <- levels(cluster_center$final_annotation)

# New clustering
p <- ggplot(meta, aes(x = hUMAP1, y = hUMAP2, color = factor(new_cluster_number, levels = as.character(seq(0, 21))))) +
        ggrastr::rasterise(geom_point(size = 0.05), dpi = 400) +
        ggrepel::geom_label_repel(
          data = cluster_center,
          aes(x = hUMAP1, y = hUMAP2, 
              label = new_cluster_number, color = factor(new_cluster_number, levels = as.character(seq(0, 21)))),
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

cowplot::save_plot("Myeloid_Tissue_legend.pdf",
       legend,
       base_height = 10,
       base_width = 10)

pdf("Myeloid_Tissue.pdf", width = 5, height = 5)
    print(p)
dev.off()