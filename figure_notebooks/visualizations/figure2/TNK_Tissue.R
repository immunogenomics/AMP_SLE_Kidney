library(singlecellmethods)
library(Seurat)
library(parallel)
library(ggplot2)
library(dplyr)
library(stringr)

meta <- readRDS("../reorderedClusters/t_nk_tissue_meta_reOrdered.Rds")
meta2 <- readRDS("../finalObjects/t_nk_tissue_meta.Rds")
# Taken from some pals color palette that I forgot the name of...
palette <- c("#782AB6", "#C075A6", "#F7E1A0", "#FC1CBF", "#FA0087", "#1CBE4F", "#FBE426",
    "#B10DA1", "#85660D", "#1C8356", "#C4451C", "#325A9B", "#F8A19F", "#AA0DFE",
    "#DEA0FD", "#2ED9FF", "#90AD1C", "#1CFFCE", "#B00068", "#FEAF16", "#3283FE",
    "#16FF32")

meta <- meta %>%
  mutate(new_cluster_number = str_extract(as.character(final_annotation), "^[^.]+"),  # Extract everything before the period
         new_cluster_number = ifelse(str_starts(new_cluster_number, "NK"), substr(new_cluster_number, 3, nchar(new_cluster_number)), substr(new_cluster_number, 2, nchar(new_cluster_number))))

cluster_center <- meta %>%
                  group_by(new_cluster_number, final_annotation) %>%
                  summarise_at(vars(hUMAP1, hUMAP2), funs(median(., na.rm=TRUE)))

cluster_center <- as.data.frame(cluster_center)
cluster_center$final_annotation <- as.factor(cluster_center$final_annotation)

levels(cluster_center$final_annotation) <- c(
  "NK0. CD56dim NK" ,
  "NK1. CD56bright NK" ,
  "T2. CD8+ GZMBhigh GZMHhigh CTL" ,
  "T3. CD8+ GZMBlow GZMHhigh CTL" ,
  "T4. CENPF+ MKI67+ Proliferating" ,
  "T5. GZMK+ CD8+ CCL5high" ,
  "T6. GZMK+ CD8+ CCL5low" ,
  "T7. GZMK+ CD8+ Effector Memory" ,
  "T8. GZMK+ CD8+ NEAT1+",
  "T9. GZMK+ CD8+ Resident Memory" ,
  "T10. GZMK+ CD8+ ITGAE+ ITGA1+",
  "T11. CD4+ Effector Memory" ,
  "T12. CD8+ GMZK+ CD69+",
  "T13. CD4+ JUNlow Resident Memory" ,
  "T14. CD4+ JUNhigh Resident Memory" ,
  "T15. CD4+ S1PR1+ Central memory/Naive" ,
  "T16. KLRB1+ KIT+ ILC" ,
  "T17. CD4+ RORC+ CCR6+ Th17" ,
  "T18. CD4+ Central Memory/Naive" ,
  "T19. CD4+ IL2RA++ FOXP3++ Treg" ,
  "T20. CD4+ FOXP3+ CXCR5+ Central Memory/Naive" ,
  "T21. CD4+ PDCD1+ CXCR5+ TFH/TPH"
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

cowplot::save_plot("TNK_Tissue_legend.pdf",
       legend,
       base_height = 10,
       base_width = 10)

pdf("TNK_Tissue.pdf", width = 5, height = 5)
    print(p)
dev.off()