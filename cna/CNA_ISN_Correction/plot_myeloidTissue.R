suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyverse)
    library(wesanderson)
    library(pals)
    library(Matrix)
})

meta <- readRDS('/data/srlab/ssg34/SLE_kidney_v2/data/qcd/Myeloid_clusterQCd_meta_harmonizedPCUMAPCellStateClusters_10042022.rds')
norm <- readRDS('/data/srlab/ssg34/SLE_kidney_v2/data/qcd/Myeloid_clusterQCd_norm_10042022.rds')

dir.create("myeloid_features")

plotGene <- function(gene, folder) {
    dir.create(paste0("myeloid_features/", folder))
    plotTable <- data.frame(
        UMAP1 = meta$hUMAP1, 
        UMAP2 = meta$hUMAP2,
        geneExpr = norm[gene, ],
        cellType = meta$final_annotation
    )

    plotTable$cellType <- factor(plotTable$cellType, 
                                            levels = c(
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

    p <- ggplot(plotTable %>% arrange(geneExpr), aes(x = UMAP1, y = UMAP2, color = geneExpr)) + 
        ggrastr::rasterise(geom_point(size = 0.1), dpi = 400) +
        scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous")) + 
        theme_classic() + ggtitle(gene)
    
    pdf(paste0("myeloid_features/", folder, "/", gene, "_featurePlot.pdf"), height = 5, width = 6)
        print(p)
    dev.off()

    p <- ggplot(data = plotTable, aes(x=cellType, y = geneExpr, fill = cellType)) + 
        geom_boxplot(notch = TRUE) + theme_classic() + 
        ggtitle(gene) +
        scale_fill_manual(values = as.vector(rev(polychrome(26)))) + 
        theme(
            legend.title = element_text(hjust = -0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position="none")

    pdf(paste0("myeloid_features/", folder, "/", gene, "_boxPlots.pdf"))
        print(p)
    dev.off()
}

genes <- read.csv("activity_globalGeneCorrelations.csv", header = TRUE, row.names = 1)
lapply(rownames(genes), plotGene, folder = "activity")

genes <- read.csv("chronicity_globalGeneCorrelations.csv", header = TRUE, row.names = 1)
lapply(rownames(genes), plotGene, folder = "chronicity")