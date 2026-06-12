# To be used with the "monocole3" environment
suppressPackageStartupMessages({
    library(slingshot)
    library(RColorBrewer)
    library(ggplot2)
    library(dplyr)
    library(igraph)
    library(ggraph)
    library(wesanderson)
    library(pals)
    library(ggnewscale)
    library(ggrepel)
})

dPCs <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/pseudotime_monoMacs/destiny/myeloid_DCs.Rds")
sds <- readRDS("destiny_slingshot_M0.Rds")
lineages <- slingLineages(sds)   # list of lineages

meta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/Myeloid_combinedMeta.csv", header = TRUE)

bloodNorm <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects_v3/myeloid_pbmc_norm.Rds")
tissueNorm <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects_v3/myeloid_tissue_norm.Rds")

allNorm <- cbind(tissueNorm, bloodNorm)

rm(bloodNorm)
rm(tissueNorm)
gc()

counters <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/Myeloid/harmony_0/pbmcVsTissue.csv")

tissueScore = counters$tissue / (counters$tissue + counters$PBMC * (length(meta$origin[meta$origin == "Tissue"]) / 
        length(meta$origin[meta$origin == "PBMC"])))

TSM <- data.frame(barcode = meta$barcode, TSM = tissueScore)

meta <- meta[meta$annotation %in% c(
  "M0. CD16+ CXC3CR1+ Monocyte",
  "M1. CD14+ CD16+ CCL2+ CX3CR1+ Monocyte",
  "M2. CD14+ CCR2+ Monocyte",
  "M3. CCL2+ CCL3+ Monocyte",
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
  "bl-M4. CD14+ CD16+ MHC2higher",
  "bl-M5. CD14+ CD16- LGALS2+",
  "bl-M3. CD16++ CD14dim CDKN1C+",
  "bl-M0. CD14+ CD16- S100Ahigh",
  "bl-M1. CD14+ CD16- CXCL8+",
  "bl-M2. CD14+ CD16- CCR2high",
  "bl-M8. CD14+ CD16+ MHC2lower",
  "bl-M6. CD14+ CD16- ISGhigh"
), ]

plot_df <- data.frame(DC1 = dPCs[,1], DC2 = dPCs[,2], TSM = TSM$TSM[TSM$barcode %in% rownames(dPCs)])

plot_df$final_annotation <- meta$annotationNumber

# dir.create("M0_PT_correlations/")
# pdf("M0_PT_correlations/testDestiny.pdf", width = 8, height = 6)
#     print(ggplot(plot_df, aes(DC1, DC2, colour = final_annotation)) +
#         ggrastr::rasterise(geom_point(size = 1), dpi = 300)) + 
#         theme_minimal(base_size = 18)
# dev.off()

# pdf("M0_PT_correlations/testDestiny_faceted.pdf", width = 12, height = 10)
#     print(ggplot(plot_df, aes(DC1, DC2, colour = final_annotation)) +
#         ggrastr::rasterise(geom_point(size = 1), dpi = 300) + 
#         theme_minimal(base_size = 18) + 
#         facet_wrap(~final_annotation)
#     )
# dev.off()


plot_df$type <- case_when(
    grepl("bl-", plot_df$final_annotation) ~ "Blood Monocyte",
    plot_df$final_annotation %in% c("M0", "M1", "M2", "M3") ~ "Tissue Monocyte",
    plot_df$final_annotation %in% c("M5", "M7", "M9", "M11") ~ "SAM",
    plot_df$final_annotation %in% c("M10", "M12") ~ "Resident Macrophage",
    TRUE ~ "Tissue Macrophage"
)

# pdf("M0_PT_correlations/testDestiny_broad.pdf", width = 8, height = 6)
#     print(ggplot(plot_df, aes(DC1, DC2, colour = type)) +
#         ggrastr::rasterise(geom_point()))
# dev.off()

# pdf("M0_PT_correlations/testDestiny_broad_faceted.pdf", width = 12, height = 8)
#   print(ggplot(plot_df, aes(DC1, DC2, colour = type)) +
#     ggrastr::rasterise(geom_point()) +
#     facet_wrap(~type))
# dev.off()

pdf("M0_PT_correlations/testDestiny_TSM.pdf", width = 6, height = 4.5)
    print(ggplot(plot_df, aes(DC1, DC2, colour = TSM)) +
        ggrastr::rasterise(geom_point(size = 2), dpi = 300) +
        theme_minimal(base_size = 14) + 
        scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous")))
dev.off()

# pdf("M0_PT_correlations/testDestiny_TSM_faceted.pdf", width = 8, height = 6)
#     print(ggplot(plot_df, aes(DC1, DC2, colour = TSM)) +
#         ggrastr::rasterise(geom_point(size = 1), dpi = 300) +
#         scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous")) +
#         facet_wrap(~type))
# dev.off()

dir.create("M0_PT_correlations/slingshot/")
mst <- slingMST(sds, as.df = TRUE)

type_cols <- brewer.pal(n = length(unique(plot_df$type)), name = "Set2")
names(type_cols) <- unique(plot_df$type)

for (i in seq_along(lineages)) {
    print(i)
    tree <- mst[mst$Lineage == i, c("DC1", "DC2", "Order", "Lineage", "Cluster")]
    tree$AnnotationNumber <- sub("\\..*$", "", tree$Cluster)

    tree$type <- case_when(
        grepl("bl-", tree$AnnotationNumber) ~ "Blood Monocyte",
        tree$AnnotationNumber %in% c("M0", "M1", "M2", "M3") ~ "Tissue Monocyte",
        tree$AnnotationNumber %in% c("M5", "M7", "M9", "M11") ~ "SAM",
        tree$AnnotationNumber %in% c("M10", "M12") ~ "Resident Macrophage",
        TRUE ~ "Tissue Macrophage"
    )

    dir.create(paste0("M0_PT_correlations/slingshot/lineage_", i))

    # Extract pseudotime
    pt <- sds@assays@data@listData$pseudotime[, i]

    # Make a data frame for ggplot
    df <- data.frame(
        DC1 = dPCs[,1],
        DC2 = dPCs[,2],
        Pseudotime = pt,
        Cluster = plot_df$final_annotation,
        TSM = TSM$TSM[TSM$barcode %in% rownames(dPCs)]
    ) %>% filter(!is.na(Pseudotime))

    stats <- cor.test(df$Pseudotime, df$TSM, method = "spearman")
    print(stats)
    
    print(cor.test(df$Pseudotime, df$TSM, method = "pearson"))
    
    # # Correlations between lineage and pseudotime
    # pdf(paste0("M0_PT_correlations/slingshot/lineage_", i, "/pseudotime_vs_TSM_lineage_", i, ".pdf"), width = 7, height = 7)
    #     print(ggplot(df, aes(x = Pseudotime, y = TSM)) +
    #         geom_point(alpha = 0.5) +
    #         geom_smooth(method = "lm", color = "blue") +
    #         theme_minimal(base_size = 18) +
    #         labs(title = paste("Lineage", i, "Pseudotime vs TSM")) +
    #         ggtitle(paste0("Lineage ", i, ";R = ", stats$estimate, "; p = ", round(stats$p.value, 2))) +
    #         xlab("Pseudotime") +
    #         ylab("TSM Score")) + 
    #         ylim(0, 1)
    # dev.off()

    # # Define a nice color palette
    # pal <- colorRampPalette(brewer.pal(11,'PRGn')[-6])(100)

    pdf(paste0("M0_PT_correlations/slingshot/lineage_", i, "/cellPlots", i, "_pseudotime.pdf"), width = 6, height = 5)
        print(
        ggplot() +
            # Pseudotime background (continuous color)
            ggrastr::rasterise(geom_point(
                data = df %>% filter(!is.na(Pseudotime)), aes(x = DC1, y = DC2, color = Pseudotime), size = 0.8, na.rm = TRUE), dpi = 300
            ) +
            scale_color_gradientn(colors = pal) +
            # Lineage path
            geom_path(
                data = tree %>% arrange(Order), aes(x = DC1, y = DC2, group = Lineage), size = 2, color = "black"
            ) +
            # 🔑 Reset color scale
            ggnewscale::new_scale_color() +
            # Tree points (discrete color)
            geom_point(
                data = tree, aes(x = DC1, y = DC2, color = type), size = 5, show.legend = TRUE
            ) +
            geom_label_repel(
                data = tree, aes(x = DC1, y = DC2, label = AnnotationNumber, color = type), fill = "white", label.size = 0.3,
                label.padding = unit(0.15, "lines"), size = 8, fontface = "bold", show.legend = FALSE
            ) + 
            scale_color_manual(values = type_cols, guide = "none") +
            theme_minimal(base_size = 14) +
            labs(title = paste("Lineage", i, "Pseudotime")) +
            theme()
        )
    dev.off()

    # # Look for genes that have the highest correlation with pseudotime
    # relevantCells <- allNorm[, colnames(allNorm) %in% rownames(df)[!is.na(pt)] ]
    # expressedGenes <- rowSums(relevantCells > 0) >= 0.1 * ncol(relevantCells)
    # geneCors <- apply(relevantCells[expressedGenes, ], 1, function(x) cor(x, df[rownames(df) %in% colnames(relevantCells), c("Pseudotime")], method = "pearson"))
    
    # pt_vals <- df[rownames(df) %in% colnames(relevantCells), c("Pseudotime")]

    # # 4. Compute p-values from t-statistic (vectorized)
    # n <- length(pt_vals)
    # t_stats <- geneCors * sqrt((n - 2) / (1 - geneCors^2))
    # pvals <- 2 * pt(-abs(t_stats), df = n - 2)

    # # 5. Adjust p-values for multiple testing (FDR)
    # fdr <- p.adjust(pvals, method = "BH")

    # # 6. Combine into a data.frame for easy ranking
    # geneStats <- data.frame(
    #     gene = names(geneCors),
    #     cor = geneCors,
    #     pval = pvals,
    #     FDR = fdr
    # )

    # # 7. Sort by correlation (decreasing)
    # geneStats <- geneStats[order(geneStats$cor, decreasing = TRUE), ]
    # write.csv(geneStats, paste0("M0_PT_correlations/slingshot/lineage_", i, "/correlatedGenes.csv"), row.names = FALSE)

    # # Plot top correlated genes
    # topPos <- head(geneStats$gene, 10)
    # topNeg <- tail(geneStats$gene, 10)
    # topGenes <- c(topPos, topNeg)

    # dir.create(paste0("M0_PT_correlations/slingshot/lineage_", i, "/topGenes/"))

    # for (gene in topGenes) {
    #     geneExp <- relevantCells[gene, ]
    #     geneDF <- df %>% filter(!is.na(Pseudotime)) %>% mutate(Expression = geneExp[rownames(df)[!is.na(pt)]])

    #     pdf(paste0("M0_PT_correlations/slingshot/lineage_", i, "/topGenes/", gene, "_expression.pdf"), width = 9, height = 7)
    #         print(ggplot() +
    #             # Gene expression background (continuous color)
    #             ggrastr::rasterise(geom_point(
    #                 data = geneDF %>% filter(!is.na(Expression)), aes(x = DC1, y = DC2, color = Expression), size = 0.8, na.rm = TRUE), dpi = 300
    #             ) +
    #             scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous")) + 
    #             # Lineage path
    #             geom_path(
    #                 data = tree %>% arrange(Order), aes(x = DC1, y = DC2, group = Lineage), size = 2, color = "black"
    #             ) +
    #             theme_minimal(base_size = 18) +
    #             labs(title = paste(gene, "Expression along Lineage", i, "(R=", paste0(round(geneStats$cor[geneStats$gene == gene], 3), ")")))
    #         )
    #     dev.off()
    # }
    # gc()
}