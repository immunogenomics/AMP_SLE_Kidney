library(singlecellmethods)
library(Seurat)
library(parallel)
library(ggplot2)
library(dplyr)
library(stringr)

meta <- readRDS("../finalObjects/t_nk_pbmc_meta.Rds")

# Taken from some pals color palette that I forgot the name of...
palette <- c("#782AB6", "#C075A6", "#F7E1A0", "#FC1CBF", "#FA0087", "#1CBE4F", "#FBE426",
    "#B10DA1", "#85660D", "#1C8356", "#C4451C", "#325A9B", "#F8A19F", "#AA0DFE",
    "#DEA0FD", "#2ED9FF", "#90AD1C", "#1CFFCE", "#B00068", "#FEAF16", "#3283FE",
    "#16FF32")

meta <- meta %>%
  mutate(new_cluster_number = str_extract(annotation, "^[^.]+"),  # Extract everything before the period
         new_cluster_number = ifelse(str_starts(new_cluster_number, "NK"), substr(new_cluster_number, 3, nchar(new_cluster_number)), substr(new_cluster_number, 2, nchar(new_cluster_number))))

cluster_center <- meta %>%
                  group_by(new_cluster_number, annotation) %>%
                  summarise_at(vars(UMAP_1, UMAP_2), funs(median(., na.rm=TRUE)))

cluster_center <- as.data.frame(cluster_center)
cluster_center$annotation <- as.factor(cluster_center$annotation)

levels(cluster_center$annotation) <- c(
        "NK0. CD56dim CD16high FGFBP2high NK", 
        "T1. CD8+ GZMHhigh FGFBP2high",
        "T2. CD8+ GZMK+ CD74+ HLA-DR+",
        "T3. CD4+ IL7Rhigh VIMhigh",
        "T4. CD8+ MT-high",
        "NK5. CD56bright NEAT1high PRF1high NK",
        "NK6. CD56bright XCL2high IL2RBhigh NK",
        "T7. CD8+ GZMK+ TEMRA", 
        "T8. CD8+ Central Memory/Naive",
        "T9. CD4+ Effector Memory",
        "T10. CD4+ Central Memory/Naive", 
        "T11. TRDC+ Gamma/Delta",
        "T12. TRGC1+ Gamma/Delta", 
        "T13. CD4+ MAF+ IT2MA+ Effector Memory",
        "T14. CD8+ GZMK+ CD74high HLA-DRhigh",
        "T15. CD8+ GZMB+ DNMT1+ HELLS+ Proliferating",
        "T16. CD4+ T-reg",
        "T17. ISGhigh",
        "T18. CD8+ GZMB+ PCNAhigh Proliferating",
        "T19. CD8+ GZMB+ CENPFhigh Proliferating"
)

meta$annotation <- as.factor(meta$annotation)
levels(meta$annotation) <- c(
        "NK0. CD56dim CD16high FGFBP2high NK", 
        "T1. CD8+ GZMHhigh FGFBP2high",
        "T2. CD8+ GZMK+ CD74+ HLA-DR+",
        "T3. CD4+ IL7Rhigh VIMhigh",
        "T4. CD8+ MT-high",
        "NK5. CD56bright NEAT1high PRF1high NK",
        "NK6. CD56bright XCL2high IL2RBhigh NK",
        "T7. CD8+ GZMK+ TEMRA", 
        "T8. CD8+ Central Memory/Naive",
        "T9. CD4+ Effector Memory",
        "T10. CD4+ Central Memory/Naive", 
        "T11. TRDC+ Gamma/Delta",
        "T12. TRGC1+ Gamma/Delta", 
        "T13. CD4+ MAF+ IT2MA+ Effector Memory",
        "T14. CD8+ GZMK+ CD74high HLA-DRhigh",
        "T15. CD8+ GZMB+ DNMT1+ HELLS+ Proliferating",
        "T16. CD4+ T-reg",
        "T17. ISGhigh",
        "T18. CD8+ GZMB+ PCNAhigh Proliferating",
        "T19. CD8+ GZMB+ CENPFhigh Proliferating"
)

# New clustering
p <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = annotation)) +
        ggrastr::rasterise(geom_point(size = 0.05), dpi = 400) +
        ggrepel::geom_label_repel(
          data = cluster_center,
          aes(x = UMAP_1, y = UMAP_2, 
              label = new_cluster_number, color = as.factor(annotation)),
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