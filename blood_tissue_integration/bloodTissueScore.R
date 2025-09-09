library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyverse)
library(wesanderson)
library(gridExtra)
library(pals)

args <- commandArgs(trailingOnly=TRUE)
name <- as.character(args[1])
theta <- as.numeric(args[2])

meta <- read.csv(paste0("PCs/", name, "_combinedMeta.csv"), header = TRUE)

umap <- readRDS(paste0(name, "/harmony_", theta, "/", name, "_umap.Rds"))
counters <- read.csv(paste0(name, "/harmony_", theta, "/pbmcVsTissue.csv"))

umap_table <- data.frame(barcode = row.names(data.frame(umap$embedding)), 
    UMAP_1 = data.frame(umap$embedding[, 1]), 
    UMAP_2 = data.frame(umap$embedding[, 2]),
    tissueScore = counters$tissue / (counters$tissue + counters$PBMC * (length(meta$origin[meta$origin == "Tissue"]) / 
        length(meta$origin[meta$origin == "PBMC"]))),
    annotation = meta$annotation,
    annotationNum = meta$annotationNumber,
    origin = meta$origin
)

colnames(umap_table) <- c("barcode", "UMAP_1", "UMAP_2", "tissueScore", "annotation", "annotationNum", "origin")

if (name == "TNK") {
    umap_table$annotation <- factor(umap_table$annotation, 
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
    umap_table$annotationNum <- factor(umap_table$annotationNum, 
        levels = c(paste0("bl-", c(paste0("NK", 0:2), paste0("T", 3:19))), paste0("NK", 0:1), paste0("T", 2:21))
    )
} else if (name == "B") {
    umap_table$annotation <- factor(umap_table$annotation, 
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
    umap_table$annotationNum <- factor(umap_table$annotationNum, 
        levels = c(paste0("bl-", c(paste0("B", 0:8))), paste0("B", 0:5), paste0("P", 6:8))
    )
} else if (name == "Myeloid") {
    umap_table$annotation <- factor(umap_table$annotation, 
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
    umap_table$annotationNum <- factor(umap_table$annotationNum, 
        levels = c(paste0("bl-", c(paste0("M", 0:8), "P9", paste0("DC", 10:12))), paste0("M", 0:15), paste0("DC", 16:20))
    )
}


# Combined
p <- ggplot() + 
    ggrastr::rasterise(geom_point(data = umap_table %>% arrange(tissueScore), mapping = aes(UMAP_1, UMAP_2, colour = tissueScore), 
        alpha = 0.75, size = 0.1), dpi = 400) + 
    theme_classic() +
    ggtitle("Tissue Score") +
    scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous"))

# pdf(paste0(name, "/harmony_", theta,"/tissueScore.pdf"), width = 7, height = 5)
#     print(p)
# dev.off()

# Split by origin
plotScores <- function(origin) {
    p <- ggplot() + 
        ggrastr::rasterise(geom_point(data = umap_table[umap_table$origin != origin, ], 
            mapping = aes(UMAP_1, UMAP_2), colour = "lightgrey", size = 0.1), dpi = 400) + 
        ggrastr::rasterise(geom_point(data = umap_table[umap_table$origin == origin, ] %>% arrange(tissueScore), 
            mapping = aes(UMAP_1, UMAP_2, colour = tissueScore), size = 0.1), dpi = 400) + 
        scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous")) +
        theme_classic() + 
        guides(colour = guide_legend(override.aes = list(size=3))) +
        ggtitle(paste("Tissue Score:", origin, "cells"))
    return(p)
}
# pdf(paste0(name, "/harmony_", theta, "/tissueScore_split.pdf"), width = 11, height = 5)
#     grid.arrange(plotScores("PBMC"), plotScores("Tissue"), nrow = 1)
# dev.off()

colorScale <- c(
    as.vector(ggthemes_data$tableau$'color-palettes'$regular$`Classic 20`[1])$value[1:length(unique(meta$annotation[meta$origin == "PBMC"]))],
    as.vector(rev(polychrome(26)))[1:length(unique(meta$annotation[meta$origin == "Tissue"]))]
)

plotBoxScores <- function() {
    medians <- umap_table %>% mutate(origin = factor(origin, c("Tissue", "PBMC"))) %>% group_by(origin, annotationNum) %>% 
        summarise(median = median(tissueScore)) %>% 
        arrange(origin, desc(median))
    write.csv(medians, paste0(name, "/harmony_", theta, "/", name, "_medians.csv"))
    p <- ggplot(data = umap_table, aes(x=annotationNum, y = tissueScore, fill = annotationNum)) + 
        geom_boxplot(notch = TRUE) + theme_classic() + 
        ggtitle(paste("Tissue Score:", name, "cells")) +
        geom_hline(yintercept=0.5, linetype="dashed", color = "black") + 
        scale_fill_manual(values = as.vector(colorScale)) + 
        theme(
            legend.title = element_text(hjust = -0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position="none") +
        scale_x_discrete(limits = medians$annotationNum)
    return(p)
}

pdf(paste0(name, "/harmony_", theta, "/boxScore_split.pdf"), width = 24, height = 5)
    print(plotBoxScores())
dev.off()

dir.create(paste0(name, "/harmony_", theta, "/tissueScore/"))

# Color each cluter by score
splitScorePlot <- function(i) {
    cellType <- unique(umap_table$annotation)[i]
    origin <- umap_table$origin[(umap_table$annotation == cellType)][1]
    print(cellType)
    p <- ggplot() + 
        ggrastr::rasterise(geom_point(data = umap_table[umap_table$annotation != cellType, ], 
            mapping = aes(UMAP_1, UMAP_2), colour = "lightgrey", size = 0.1), dpi = 400) + 
        ggrastr::rasterise(geom_point(data = umap_table[umap_table$annotation == cellType, ] %>% arrange(tissueScore), 
            mapping = aes(UMAP_1, UMAP_2, colour = tissueScore), size = 0.1), dpi = 400) + 
        scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous")) +
        theme_classic() + 
        theme(legend.position="none") +
        ggtitle(paste0(origin, ": ", cellType))
    pdf(paste0(name, "/harmony_", theta, "/tissueScore/", i, "_score.pdf"), width = 5, height = 5)
        print(p)
    dev.off()
}
lapply(c(1:length(unique(umap_table$annotation))), splitScorePlot)

percentages <- umap_table %>% group_by(annotation) %>% summarise(
    underFirst = sum(tissueScore < 0.25) /  n(), 
    underSecond = sum(tissueScore < 0.50) /  n(), 
    overSecond = sum(tissueScore > 0.50) /  n(), 
    overThird = sum(tissueScore > 0.75) /  n()
)

percentages <- as.data.frame(percentages)

bloodBars <- function(category) {
    plotPerc_pbmcs <- percentages[percentages$annotation %in% unique(umap_table$annotation[umap_table$origin == "PBMC"]), ]
    plotPerc_pbmcs <- plotPerc_pbmcs[rev(order(plotPerc_pbmcs[[category]])), ]

    pb <- ggplot(data = plotPerc_pbmcs, aes(x = annotation, y = .data[[category]])) + 
        geom_bar(stat = "identity") + theme_classic() + 
    ggtitle(paste("PBMC:", category, "Percentile")) +
    theme(
        legend.title = element_text(hjust = -0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.margin = unit(c(0.5,0.5,0.5,3), "cm")) +
    scale_x_discrete(limits = plotPerc_pbmcs$annotation)

    plotPerc_tissue <- percentages[percentages$annotation %in% unique(umap_table$annotation[umap_table$origin == "Tissue"]), ]
    plotPerc_tissue <- plotPerc_tissue[rev(order(plotPerc_tissue[[category]])), ]

    pt <- ggplot(data = plotPerc_tissue, aes(x = annotation, y = .data[[category]])) + 
        geom_bar(stat = "identity") + theme_classic() + 
    ggtitle(paste("Tissue:", category, "Percentile")) +
    theme(
        legend.title = element_text(hjust = -0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_x_discrete(limits = plotPerc_tissue$annotation)

    pdf(paste0(name, "/harmony_", theta, "/", category, "_percentile.pdf"), width = 11, height = 5)
        grid.arrange(pb, pt, nrow = 1)
    dev.off()
}

lapply(c("underFirst", "underSecond", "overSecond", "overThird"), bloodBars)