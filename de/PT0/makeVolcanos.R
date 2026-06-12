library(tidyverse)

tests <- readRDS("pairwise_geneExpresion_v2.rds")

completeVsNone_all <- tests$completeVsNone
completeVsPartial_all <- tests$completeVsPartial
partialVsNone_all <- tests$partialVsNone

completeVsNone_all$padj = p.adjust(completeVsNone_all$p, method = "BH")
completeVsPartial_all$padj = p.adjust(completeVsPartial_all$p, method = "BH")
partialVsNone_all$padj = p.adjust(partialVsNone_all$p, method = "BH")

# Volcano plots of genes
p <- ggplot(completeVsNone_all) +
    ggrastr::rasterise(geom_point(aes(x = beta, y = -log10(p)), color = 'grey'), dpi = 300) + 
    theme_classic(base_size = 14) +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 11, hjust=0.5),
        axis.text.y = element_text(color = "black", size = 11, hjust=1),
        axis.title = element_text(size=14, hjust = 0.5),
        plot.title = element_text(size=14, hjust = 0.5),
    ) +
    geom_hline(yintercept = -log10(0.05/nrow(completeVsNone_all)), color = 'gray55', linewidth = 1) +
    labs(x = 'Beta Complete', y = '-Log10(P-value)', title = "Complete vs Non Responders") +
    scale_color_identity()

ggsave(
  filename = "volcano_plot_completeVsNone_3.pdf",
  plot = p,
  width = 3,
  height = 3
)

p <- ggplot(completeVsPartial_all) +
    ggrastr::rasterise(geom_point(aes(x = beta, y = -log10(p)), color = 'grey'), dpi = 300) + 
    theme_classic(base_size = 14) +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 11, hjust=0.5),
        axis.text.y = element_text(color = "black", size = 11, hjust=1),
        axis.title = element_text(size=14, hjust = 0.5),
        plot.title = element_text(size=14, hjust = 0.5),
    ) +
    geom_hline(yintercept = -log10(0.05/nrow(completeVsPartial_all)), color = 'gray55', linewidth = 1) +
    labs(x = 'Beta Complete', y = '-Log10(P-value)', title = "Complete vs Partial Responders") +
    scale_color_identity()

ggsave(
  filename = "volcano_plot_completeVsPartial_3.pdf",
  plot = p,
  width = 3,
  height = 3
)

p <- ggplot(partialVsNone_all) +
    ggrastr::rasterise(geom_point(aes(x = beta, y = -log10(p)), color = 'grey'), dpi = 300) + 
    theme_classic(base_size = 14) +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 11, hjust=0.5),
        axis.text.y = element_text(color = "black", size = 11, hjust=1),
        axis.title = element_text(size=14, hjust = 0.5),
        plot.title = element_text(size=14, hjust = 0.5),
    ) +
    geom_hline(yintercept = -log10(0.05/nrow(partialVsNone_all)), color = 'gray55', linewidth = 1) +
    labs(x = 'Beta Partial', y = '-Log10(P-value)', title = "Partial vs Non Responders") +
    scale_color_identity()

ggsave(
  filename = "volcano_plot_partialVsNone_3.pdf",
  plot = p,
  width = 3,
  height = 3
)