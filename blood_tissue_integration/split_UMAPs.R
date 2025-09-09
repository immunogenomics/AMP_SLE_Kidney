# Use with the sc environment
suppressPackageStartupMessages({
    library(uwot)
    library(ggplot2)
    library(tidyverse)
    library(dplyr)
    library(ggrastr)
    library(pals)
    library(ggthemes)
})

set.seed(0)

args <- commandArgs(trailingOnly=TRUE)
name <- as.character(args[1])
theta <- as.numeric(args[2])

# Read in data
combined_meta <- read.csv(paste0("PCs/", name, "_combinedMeta.csv"), header = TRUE)
umap_res <- readRDS(paste0(name, "/harmony_", theta, "/", name, "_umap.Rds"))

umap_table <- data.frame(barcode = row.names(data.frame(umap_res$embedding)), 
    UMAP_1 = data.frame(umap_res$embedding[, 1]), 
    UMAP_2 = data.frame(umap_res$embedding[, 2]),
    cellType = combined_meta$annotation,
    origin = combined_meta$origin
)
colnames(umap_table) <- c("barcode", "UMAP_1", "UMAP_2", "cellType", "origin")

if (name == "TNK") {
    umap_table$cellType[umap_table$cellType == "T1. CD8+ GZMBhigh GZMHhigh CTL"] <- "T1. CD8+ GZMB+ CTL"
    umap_table$cellType[umap_table$cellType == "T6. GZMK+ CD8+ CCL5low"] <- "T6. GZMK+ CD8+ NKG7low"
    umap_table$cellType[umap_table$cellType == "T5. GZMK+ CD8+ CCL5high"] <- "T5. GZMK+ CD8+ NKG7high"
    umap_table$cellType[umap_table$cellType == "T20. CD4+ FOXP3+ CXCR5+ Central Memory/Naive"] <- "T20. CD4+ FOXP3+ Central Memory/Naive"
    umap_table$cellType[umap_table$cellType == "T10. GZMK+ CD8+ ITGAE+ ITGA1+"] <- "T10. GZMK+ CD8+ ITGAE"
    umap_table$cellType[umap_table$cellType == "T2. CD8+ GZMBlow GZMHhigh CTL"] <- "T2. CD8+ GZMB+ SYNE2bright CTL"

    umap_table$cellType <- factor(umap_table$cellType, 
        levels = c(
        # Tissue
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
        "T19. CD8+ GZMB+ CENPFhigh Proliferating",
        # PBMCs
        "NK0. CD56dim NK",
        "NK3. CD56bright NK",
        "T1. CD8+ GZMB+ CTL",
        "T2. CD8+ GZMB+ SYNE2bright CTL",
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
        "T20. CD4+ FOXP3+ Central Memory/Naive",
        "T21. CD4+ PDCD1+ CXCR5+ TFH/TPH"
        ))
        tissueBreak <- 20
} else if (name == "B") {
    umap_table$cellType <- factor(umap_table$cellType, 
                                            levels = c(
                                            # PBMCs
                                            "B0. CXCR5high Naive",
                                            "B1. CRIP1+ Mature",
                                            "B2. FCRL5+ ITGAX+ ABC",
                                            "B3. CD79A+ VREPB3+ pre B Cell",
                                            "B4. ISGhigh Naive",
                                            "B5. MThigh",
                                            "B6. JUN+ NFKB1+ Activated",
                                            "B7. IGKC+ XBP1+ Plasma-like", 
                                            "B8. FCRL5+ IGHG1+ B Cell",
                                            # Tissue
                                            "B0. FOXO1+ BCL6+ GC",
                                            "B1. CD28+ IGM- activated B Cell",
                                            "B2. IGHD+ FCER2+ Naïve B Cell",
                                            "B3. BCL2+ CD27+ MCL1+ Unswitched Memory B Cell",
                                            "B4. BCL2+ BCL11A+ IGHE+ Memory B Cell",
                                            "B5. FCRL5+ ITGAX+ TBX21+ ABC-like", 
                                            "P6. CD38++ MKI67+ Plasmablast",
                                            "P7. IGG+ IGKC++ IGL- Plasma Cell",
                                            "P8. IGA+ IGL- Plasma Cell",
                                            "P9. IGG+ IGKC+ IGL- Plasma Cell",
                                            "P10. IGA+ IGL+ Plasma Cell",
                                            "P11. IGG+ IGKC- IGL+ Plasma Cell"))
    tissueBreak <- 9
} else if (name == "Myeloid") {
    umap_table$cellType <- factor(umap_table$cellType, 
                                            levels = c(
                                            # PBMCs
                                            "M0. CD14+ CD16- S100Ahigh", 
                                            "M1. CD14+ CD16- CXCL8+",
                                            "M2. CD14+ CD16- CCR2high",
                                            "M3. CD16++ CD14dim CDKN1C+",
                                            "M4. CD14+ CD16+ MHC2higher",
                                            "M5. CD14+ CD16- LGALS2+",
                                            "DC6. CLEC10A+ CD1C+ DC2",
                                            "M7. CD14+ CD16- ISGhigh", 
                                            "P8. PPBP+ GP1BB+ Platelet",
                                            "DC9. TCF4+ CLEC4C+ pDC",
                                            "M10. MThigh",
                                            "M11. CD14+ CD16+ MHC2lower",
                                            "DC12. CLEC9A+ XCR1+ DC1",
                                            # Tissue
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
                                            'DC13. CCR7+ LAMP3+ DC2', 
                                            'M14. APOChigh C3high Macrophage',
                                            'DC15. CLEC10Alow cDC2',
                                            'M16. APOClow C3high Macrophage',
                                            'DC17. CLEC10Ahigh cDC2',
                                            'M18. CENPF+ MKI67+ Proliferating',
                                            'DC19. cDC1',
                                            'DC20. pDC'))
    tissueBreak <- 13
}

# Tissue Plot

dir.create(paste0(name, "/harmony_", theta, "/tissue"))
dir.create(paste0(name, "/harmony_", theta, "/PBMC"))

pbmcPlot <- function(i) {
    cellType <- levels(umap_table$cellType)[i]
    print(cellType)
    p <- ggplot() + 
        ggrastr::rasterise(geom_point(data = umap_table[umap_table$cellType != cellType, ], 
            mapping = aes(UMAP_1, UMAP_2), colour = "lightgrey", size = 0.1), dpi = 400) + 
        ggrastr::rasterise(geom_point(data = umap_table[umap_table$cellType == cellType, ], 
            mapping = aes(UMAP_1, UMAP_2, colour = cellType), size = 0.1), dpi = 400) + 
        scale_colour_manual(values=as.vector(ggthemes_data$tableau$'color-palettes'$regular$`Classic 20`[1])$value[i]) +
        theme_classic() + 
        theme(legend.position="none") +
        ggtitle(cellType)
    pdf(paste0(name, "/harmony_", theta, "/PBMC/", i, ".pdf"), width = 5, height = 5)
        print(p)
    dev.off()
}

tissuePlot <- function(i) {
    cellType <- levels(umap_table$cellType)[i]
    print(cellType)
    p <- ggplot() + 
        ggrastr::rasterise(geom_point(data = umap_table[umap_table$cellType != cellType, ], 
            mapping = aes(UMAP_1, UMAP_2), colour = "lightgrey", size = 0.1), dpi = 400) + 
        ggrastr::rasterise(geom_point(data = umap_table[umap_table$cellType == cellType, ], 
            mapping = aes(UMAP_1, UMAP_2, colour = cellType), size = 0.1), dpi = 400) + 
        scale_colour_manual(values=as.vector(rev(polychrome(26)))[i - tissueBreak]) +
        theme_classic() + 
        theme(legend.position="none") +
        ggtitle(cellType)
    pdf(paste0(name, "/harmony_", theta, "/tissue/", i, ".pdf"), width = 5, height = 5)
        print(p)
    dev.off()
}
lapply(c(1:tissueBreak), pbmcPlot)


lapply(c((tissueBreak + 1):length(unique(umap_table$cellType))), tissuePlot)

if (name == "TNK") {
    pbmcPlot_manual <- function(i) {
        cellType <- levels(umap_table$cellType)[i]
        print(cellType)
        p <- ggplot() + 
            ggrastr::rasterise(geom_point(data = umap_table[umap_table$cellType != cellType, ], 
                mapping = aes(UMAP_1, UMAP_2), colour = "lightgrey", size = 0.1), dpi = 400) + 
            ggrastr::rasterise(geom_point(data = umap_table[umap_table$cellType == cellType, ], 
                mapping = aes(UMAP_1, UMAP_2, colour = cellType), size = 0.1), dpi = 400) + 
            scale_colour_manual(values="red") +
            theme_classic() + 
            theme(legend.position="none") +
            ggtitle(cellType)
        pdf(paste0(name, "/harmony_", theta, "/PBMC/", i, "_manual.pdf"), width = 5, height = 5)
            print(p)
        dev.off()
    }
    
    pbmcPlot_manual(16)
}