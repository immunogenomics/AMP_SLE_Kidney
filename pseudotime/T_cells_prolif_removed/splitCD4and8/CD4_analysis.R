# To be used with the "monocole3" environment
suppressPackageStartupMessages({
    library(slingshot)
    library(RColorBrewer)
    library(ggplot2)
    library(dplyr)
    library(igraph)
    library(ggraph)
    library(wesanderson)
    library(ggrepel)
})
dPCs <- readRDS("pt_outputs/CD4_diffusionMap_eigenvectors.Rds")
sds <- readRDS("pt_outputs/CD4_destiny_slingshot.Rds")

bloodNorm <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects_v3/t_nk_pbmc_norm.Rds")
tissueNorm <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects_v3/t_nk_tissue_norm.Rds")

allNorm <- cbind(tissueNorm, bloodNorm)

meta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/TNK_combinedMeta.csv")
meta_4 <- read.csv("CD4_inputs/CD4_T_cells_meta.csv", row.names = 1)

lineages <- slingLineages(sds)   # list of lineages
counters <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/TNK/harmony_0/pbmcVsTissue.csv")

tissueScore = counters$tissue / (counters$tissue + counters$PBMC * (length(meta$origin[meta$origin == "Tissue"]) / 
        length(meta$origin[meta$origin == "PBMC"])))

TSM <- data.frame(barcode = meta$barcode, TSM = tissueScore)

# Take CD8 cells only
plot_df <- data.frame(DC1 = dPCs[,1], DC2 = dPCs[,2], TSM = TSM$TSM[na.omit(match(meta_4$barcode, TSM$barcode))])

plot_df$final_annotation <- meta_4$annotationNumber

# Step 3: Reorder allNorm columns to match plot_df
allNorm <- allNorm[, meta_4$barcode]

dir.create("figures/")

# pdf("figures/CD4_destiny.pdf", width = 8, height = 6)
#     print(ggplot(plot_df, aes(DC1, DC2, colour = final_annotation)) +
#         ggrastr::rasterise(geom_point(size = 1), dpi = 300) + 
#         theme_minimal(base_size = 18)) 
# dev.off()

# pdf("figures/CD4_destiny_faceted.pdf", width = 12, height = 10)
#     print(ggplot(plot_df, aes(DC1, DC2, colour = final_annotation)) +
#         ggrastr::rasterise(geom_point(size = 1), dpi = 300) + 
#         theme_minimal(base_size = 18) + 
#         facet_wrap(~final_annotation)
#     )
# dev.off()

# pdf("figures/CD4_destiny_TSM.pdf", width = 8, height = 6)
#     print(ggplot(plot_df, aes(DC1, DC2, colour = TSM)) +
#         ggrastr::rasterise(geom_point(size = 2), dpi = 300) +
#         scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous")))
# dev.off()

mst <- slingMST(sds, as.df = TRUE)

dir.create("slingshot_CD4")
for (i in seq_along(lineages)) {
    tree <- mst[mst$Lineage == i, c("DC1", "DC2", "Order", "Lineage", "Cluster")]
    tree$AnnotationNumber <- sub("\\..*$", "", tree$Cluster)

    dir.create(paste0("slingshot_CD4/lineage_", i))

    # Extract pseudotime
    pt <- sds@assays@data@listData$pseudotime[, i]
    pt <- pt[meta_4$barcode]

    # Make a data frame for ggplot
    df <- data.frame(
        DC1 = plot_df$DC1,
        DC2 = plot_df$DC2,
        Pseudotime = pt,
        Cluster = plot_df$final_annotation,
        TSM = plot_df$TSM
    )

    # Define a nice color palette
    pal <- colorRampPalette(brewer.pal(11,'PRGn')[-6])(100)

    stats <- cor.test(df$Pseudotime, df$TSM, method = "spearman")
    # Correlations between lineage and pseudotime
    # Define a nice color palette
    pal <- colorRampPalette(brewer.pal(11,'PRGn')[-6])(100)

    pdf(paste0("slingshot_CD4/lineage_", i, "/cellPlots", i, "_pseudotime.pdf"), width = 6, height = 5)
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
            # Tree points (discrete color)
            geom_point(
                data = tree, aes(x = DC1, y = DC2), color = "black", size = 5, show.legend = TRUE
            ) +
            geom_label_repel(
                data = tree, aes(x = DC1, y = DC2, label = AnnotationNumber), color = "black", fill = "white", label.size = 0.3,
                label.padding = unit(0.15, "lines"), size = 5, fontface = "bold", show.legend = FALSE, max.overlaps = Inf, force = 4
            ) + 
            theme_minimal(base_size = 18) +
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
    # write.csv(geneStats, paste0("slingshot_CD4/lineage_", i, "/correlatedGenes.csv"), row.names = FALSE)

    # # Plot top correlated genes
    # topPos <- head(geneStats$gene, 10)
    # topNeg <- tail(geneStats$gene, 10)
    # topGenes <- c(topPos, topNeg)

    # dir.create(paste0("slingshot_CD4/lineage_", i, "/topGenes/"))

    # for (gene in topGenes) {
    #     geneExp <- relevantCells[gene, ]
    #     geneDF <- df %>% filter(!is.na(Pseudotime)) %>% mutate(Expression = geneExp[rownames(df)[!is.na(pt)]])

    #     pdf(paste0("slingshot_CD4/lineage_", i, "/topGenes/", gene, "_expression.pdf"), width = 9, height = 7)
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