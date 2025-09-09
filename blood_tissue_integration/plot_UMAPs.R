# Use with the sc environment
suppressPackageStartupMessages({
    library(uwot)
    library(ggplot2)
    library(tidyverse)
    library(dplyr)
    library(ggrastr)
    library(pals)
    library(ggthemes)
    library(cowplot)
})

set.seed(0)

args <- commandArgs(trailingOnly=TRUE)
name <- as.character(args[1])
theta <- as.numeric(args[2])

# Read in data
harmony <- readRDS(paste0(name, "/harmony_", theta, "/", name, "_combined_hPCs.Rds"))
combined_meta <- read.csv(paste0("PCs/", name, "_combinedMeta.csv"), header = TRUE)

# umap_res <- umap(harmony, n_neighbors = 30L, metric = "euclidean", min_dist = .1, pca = 20, 
#   ret_model = TRUE, seed = 1)

# saveRDS(umap_res, paste0(name, "/harmony_", theta, "/", name, "_umap.Rds"))

umap_res <- readRDS(paste0(name, "/harmony_", theta, "/", name, "_umap.Rds"))

umap_table <- data.frame(barcode = row.names(data.frame(umap_res$embedding)), 
    UMAP_1 = data.frame(umap_res$embedding[, 1]), 
    UMAP_2 = data.frame(umap_res$embedding[, 2]),
    cellType = combined_meta$annotation,
    origin = combined_meta$origin
)
colnames(umap_table) <- c("barcode", "UMAP_1", "UMAP_2", "cellType", "origin")

if (name == "TNK") {
    umap_table$cellType <- factor(umap_table$cellType, 
        levels = c(
        # PBMC
        "bl-NK0. CD56dim CD16high FGFBP2high NK",
        "bl-NK1. CD56bright NEAT1high PRF1high NK",
        "bl-NK2. CD56bright XCL2high IL2RBhigh NK",
        "bl-T3. CD8+ GZMHhigh FGFBP2high",
        "bl-T4. CD8+ GZMK+ CD74+ HLA-DR+",
        "bl-T5. CD4+ IL7Rhigh VIMhigh",
        "bl-T6. CD8+ MT-high",
        "bl-T7. CD8+ GZMK+ TEMRA",
        "bl-T8. CD8+ Central Memory/Naive",
        "bl-T9. CD4+ Effector Memory",
        "bl-T10. CD4+ Central Memory/Naive",
        "bl-T11. TRDC+ Gamma/Delta",
        "bl-T12. TRGC1+ Gamma/Delta",
        "bl-T13. CD4+ MAF+ IT2MA+ Effector Memory",
        "bl-T14. CD8+ GZMK+ CD74high HLA-DRhigh",
        "bl-T15. CD8+ GZMB+ DNMT1+ HELLS+ Proliferating",
        "bl-T16. CD4+ T-reg",
        "bl-T17. ISGhigh",
        "bl-T18. CD8+ GZMB+ PCNAhigh Proliferating",
        "bl-T19. CD8+ GZMB+ CENPFhigh Proliferating",
        # Tissue
        "NK0. CD56dim NK",
        "NK1. CD56bright NK",
        "T2. CD8+ GZMB+ CTL",
        "T3. CD8+ GZMB+ SYNE2bright CTL",
        "T4. CENPF+ MKI67+ Proliferating",
        "T5. GZMK+ CD8+ NKG7high",
        "T6. GZMK+ CD8+ NKG7low",
        "T7. GZMK+ CD8+ Effector Memory",
        "T8. GZMK+ CD8+ NEAT1+",
        "T9. GZMK+ CD8+ Resident Memory",
        "T10. GZMK+ CD8+ ITGAE",
        "T11. CD4+ Effector Memory",
        "T12. CD8+ GMZK+ CD69+",
        "T13. CD4+ JUNlow Resident Memory",
        "T14. CD4+ JUNhigh Resident Memory",
        "T15. CD4+ S1PR1+ Central memory/Naive",
        "T16. KLRB1+ KIT+ ILC",
        "T17. CD4+ RORC+ CCR6+ Th17",
        "T18. CD4+ Central Memory/Naive",
        "T19. CD4+ IL2RA++ FOXP3++ Treg",
        "T20. T20. CD4+ FOXP3+ Central Memory/Naive",
        "T21. CD4+ PDCD1+ CXCR5+ TFH/TPH"
        ))
        tissueBreak <- 22
} else if (name == "B") {
    umap_table$cellType <- factor(umap_table$cellType, 
                                            levels = c(
                                            # PBMCs
                                            "bl-B0. CXCR5high Naive",
                                            "bl-B1. CRIP1+ Mature",
                                            "bl-B2. FCRL5+ ITGAX+ ABC",
                                            "bl-B3. CD79A+ VREPB3+ pre B Cell",
                                            "bl-B4. ISGhigh Naive",
                                            "bl-B5. MThigh",
                                            "bl-B6. JUN+ NFKB1+ Activated",
                                            "bl-B7. IGKC+ XBP1+ Plasma-like",
                                            "bl-B8. FCRL5+ IGHG1+ B Cell",
                                            # Tissue
                                            "B0. FOXO1+ BCL6+ GC",
                                            "B1. CD28+ IGM- activated B Cell",
                                            "B2. IGHD+ FCER2+ Naïve B Cell",
                                            "B3. BCL2+ CD27+ MCL1+ Unswitched Memory B Cell",
                                            "B4. BCL2+ BCL11A+ IGHE+ Memory B Cell",
                                            "B5. FCRL5+ ITGAX+ TBX21+ ABC-like", 
                                            "P6. CD38++ MKI67+ Plasmablast",
                                            "P7. IGG+ Plasma Cell",
                                            "P8. IGA+ Plasma Cell"))
} else if (name == "Myeloid") {
    umap_table$cellType <- factor(umap_table$cellType, 
                                            levels = c(
                                            # PBMCs
                                            "bl-M0. CD14+ CD16- S100Ahigh",
                                            "bl-M1. CD14+ CD16- CXCL8+",
                                            "bl-M2. CD14+ CD16- CCR2high",
                                            "bl-M3. CD16++ CD14dim CDKN1C+",
                                            "bl-M4. CD14+ CD16+ MHC2higher",
                                            "bl-M5. CD14+ CD16- LGALS2+",
                                            "bl-M6. CD14+ CD16- ISGhigh",
                                            "bl-M7. MThigh",
                                            "bl-M8. CD14+ CD16+ MHC2lower",
                                            "bl-P9. PPBP+ GP1BB+ Platelet",
                                            "bl-DC10. CLEC10A+ CD1C+ DC2",
                                            "bl-DC11. TCF4+ CLEC4C+ pDC",
                                            "bl-DC12. CLEC9A+ XCR1+ DC1",
                                            # Tissue
                                            "M0. CD16+ CXC3CR1+ Monocyte",
                                            "M1. CD14+ CD16+ CCL2+ CX3CR1+ Monocyte",
                                            "M2. CD14+ CCR2+ Monocyte",
                                            "M3. CCL2+ CCL3+ Monocyte",
                                            "M4. TPSB2+ MAST cell",
                                            "M5. GPNMBhigh NUPR1high Macrophage",
                                            "M6. SELENOPinter ISGhigh Macrophage",
                                            "M7. SPP1high FABP5high Macrophage",
                                            "M8. SPP1low FABP5high Macrophage",
                                            "M9. MERTKhigh FABP5high Macrophage",
                                            "M10. SELENOPinter LYVE1inter Resident Macrophage",
                                            "M11. GPMNBhigh NUPR1low Macrophage",
                                            "M12. SELENOPhigh LYVE1high Resident Macrophage",
                                            "M13. APOChigh C3high Macrophage",
                                            "M14. APOClow C3high Macrophage",
                                            "M15. CENPF+ MKI67+ Proliferating",
                                            "DC16. CCR7+ LAMP3+ DC2",
                                            "DC17. CLEC10Alow cDC2",
                                            "DC18. CLEC10Ahigh cDC2",
                                            "DC19. cDC1",
                                            "DC20. pDC"))
}

# Tissue Plot
p <- ggplot() + 
    ggrastr::rasterise(geom_point(data = umap_table[umap_table$origin == "Tissue", ], 
        mapping = aes(UMAP_1, UMAP_2, colour = cellType), size = 0.1), dpi=400) + 
    scale_colour_manual(values = as.vector(rev(polychrome(26)))) + 
    theme_classic() + 
    theme(
        legend.position = "right",
              axis.title = element_text(hjust = 0.75, 
                                        size = 20, 
                                        face = "bold"), 
              axis.text = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              legend.text = element_text(size = 20),
              legend.title = element_blank()
          
    ) +
    labs(x = "UMAP1", y = "UMAP2") +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 15)))

legend <- cowplot::get_legend(p)
p <- p + theme(legend.position = "none")
legend <- cowplot::plot_grid(legend)

pdf(paste0(name, "/harmony_", theta, "/tissuePlot.pdf"), width = 10, height = 10)
    print(p)
dev.off()

cowplot::save_plot(paste0(name, "/harmony_", theta, "/tissuePlot_legend.pdf"),
       legend,
       base_height = 7,
       base_width = 14)

# PBMC Plot
p <- ggplot() + 
    ggrastr::rasterise(geom_point(data = umap_table[umap_table$origin == "PBMC", ], 
        mapping = aes(UMAP_1, UMAP_2, colour = cellType), size = 0.1), dpi=400) + 
    ggthemes::scale_color_tableau(palette = "Classic 20") + 
    theme_classic() + 
    theme(
        legend.position = "right",
              axis.title = element_text(hjust = 0.75, 
                                        size = 20, 
                                        face = "bold"), 
              axis.text = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              legend.text = element_text(size = 20),
              legend.title = element_blank(),
              plot.title = element_text(size = 20, face = "bold")
          
    ) +
    labs(x = "UMAP1", y = "UMAP2") +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 15)))

legend <- cowplot::get_legend(p)
p <- p + theme(legend.position = "none")
legend <- cowplot::plot_grid(legend)

pdf(paste0(name, "/harmony_", theta, "/pbmcPlot.pdf"), width = 10, height = 10)
    print(p)
dev.off()

cowplot::save_plot(paste0(name, "/harmony_", theta, "/pbmcPlot_legend.pdf"),
       legend,
       base_height = 7,
       base_width = 14)

# Origin Plot
p <- ggplot() + 
    ggrastr::rasterise(geom_point(data = umap_table, 
        mapping = aes(UMAP_1, UMAP_2, colour = origin), size = 0.1), dpi=400) + 
    ggthemes::scale_color_tableau(palette = "Superfishel Stone") + 
    theme_classic() + 
    theme(
        legend.position = "right",
              axis.title = element_text(hjust = 0.75, 
                                        size = 20, 
                                        face = "bold"), 
              axis.text = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              legend.text = element_text(size = 20),
              legend.title = element_blank()
          
    ) +
    labs(x = "UMAP1", y = "UMAP2") +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 15)))

legend <- cowplot::get_legend(p)
p <- p + theme(legend.position = "none")
legend <- cowplot::plot_grid(legend)

pdf(paste0(name, "/harmony_", theta, "/originPlot.pdf"), width = 10, height = 10)
    print(p)
dev.off()

cowplot::save_plot(paste0(name, "/harmony_", theta, "/originPlot_legend.pdf"),
       legend,
       base_height = 7,
       base_width = 14)