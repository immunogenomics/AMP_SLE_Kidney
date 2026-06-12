# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(fgsea)
library(msigdbr)
library(ggrepel)

paths <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/data/tissue/myeloid_SAMvsNonSAM_gsea_20260216.rds")

# Filter significant pathways and restrict to HALLMARK + GOBP
sig <- paths %>%
  filter(
    padj < 0.05,
    NES > 0,
    str_detect(pathway, "^HALLMARK_|^GOBP_")
  ) %>%
  arrange(desc(NES))

unique(sig$pathway)
sig$leadingEdge[sig$pathway == "HALLMARK_CHOLESTEROL_HOMEOSTASIS"]
sig$leadingEdge[sig$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
sig$leadingEdge[sig$pathway == "HALLMARK_COAGULATION"]
sig$leadingEdge[sig$pathway == "HALLMARK_COMPLEMENT"]

pdf("samNonSAM_barplot.pdf", width = 10, height = 6)
    ggplot(sig, aes(x = reorder(pathway, NES), y = NES, fill = -log10(padj))) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_classic(base_size = 14) +
        labs(x = "Pathway", y = "Normalized Enrichment Score (NES)", title = "GSEA Results for SAM vs non-SAM") +
        scale_fill_viridis_c()
dev.off()

res = res_all_wide

min_p = res[res['p']!= 0, 'p'] %>% min # Replace zeros for visualization
res = res %>% mutate(p = ifelse(p==0, min_p, p))

write.csv(res, "myeloid_SAMvsOtherMyeloid_differential_expression_lm_20260216.csv", row.names = FALSE)

res <- read.csv("myeloid_SAMvsOtherMyeloid_differential_expression_lm_20260216.csv")

res <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/data/tissue/myeloid_SAMvsOtherMyeloid_differential_expression_lm_20260216.rds")

highlightedGenes <- unique(c(
  "SPP1", "FABP5", "TREM2", "CD9", "CD63",
  unlist(sig$leadingEdge[sig$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"])
))

res <- res %>%
  mutate(label_color = case_when(
    gene %in% highlightedGenes ~ "maroon",
    TRUE ~ "grey"
  ))

# Make maroon points plot last
res <- res %>%
  mutate(label_color = factor(label_color, levels = c("maroon", "grey")))

res <- as.data.frame(res)

print(highlightedGenes)

set.seed(1)
p <- ggplot() +
    ggrastr::rasterise(geom_point(data = res %>% filter(label_color == "grey"), aes(x = beta, y = -log10(p), color = label_color)), dpi = 300) + 
    ggrastr::rasterise(geom_point(data = res %>% filter(label_color == "maroon"), aes(x = beta, y = -log10(p), color = label_color)), dpi = 300) + 
    theme_classic(base_size = 14) +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 11, hjust=0.5),
        axis.text.y = element_text(color = "black", size = 11, hjust=1),
        axis.title = element_text(size=14, hjust = 0.5),
        plot.title = element_text(size=14, hjust = 0.5),
    ) +
    geom_hline(yintercept = -log10(0.05/nrow(res)), color = 'gray55', linewidth = 1) +
    labs(x = 'Beta SAM', y = '-Log10(P-value)', title = "SAMs vs nonSAMs") +
    geom_text_repel(data = res %>% filter(gene %in% highlightedGenes), 
                aes(x = beta, y = -log10(p), label = gene, color = label_color),
                size = 5, fontface = "bold",
                max.overlaps = 10, force = 4, min.segment.length = 0) +
    geom_text_repel(data = res %>% filter(gene %in% c("ITGB1", "TGFBI", "ITGAV")), 
                aes(x = beta, y = -log10(p), label = gene, color = label_color),
                size = 5, fontface = "bold",
                max.overlaps = 10, force = 4, min.segment.length = 0) +
    scale_color_identity()

ggsave(
  filename = "volcano_plot.pdf",
  plot = p,
  width = 5,
  height = 5
)

# For each gene, list which significant pathways it appears in as leading edge
gene_to_pathways <- tapply(
  rep(sig$pathway, lengths(sig$leadingEdge)),
  unlist(sig$leadingEdge),
  function(pw) paste(pw, collapse = "; ")
)

res$leadingEdgeIn <- gene_to_pathways[res$gene]
res$leadingEdgeIn[is.na(res$leadingEdgeIn)] <- ""

write.csv(res, "myeloid_SAMvsOtherMyeloid_DE_withPathways.csv", row.names = FALSE)
