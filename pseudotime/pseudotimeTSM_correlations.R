# To be used with the "monocole3" environment
library(slingshot)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)
library(wesanderson)
library(pals)

dPCs <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/pseudotime_monoMacs/destiny/myeloid_DCs.Rds")
meta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/Myeloid_combinedMeta.csv")

rightClusters <- c("M0. CD16+ CXC3CR1+ Monocyte",
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
  "bl-M0. CD14+ CD16- S100Ahigh",
  "bl-M1. CD14+ CD16- CXCL8+",
  "bl-M2. CD14+ CD16- CCR2high",
  "bl-M3. CD16++ CD14dim CDKN1C+",
  "bl-M4. CD14+ CD16+ MHC2higher",
  "bl-M5. CD14+ CD16- LGALS2+",
  "bl-M6. CD14+ CD16- ISGhigh",
  "bl-M8. CD14+ CD16+ MHC2lower")

meta <- meta[meta$annotation %in% rightClusters, ]

# sds <- slingshot(dPCs, clusterLabels = meta$annotation, start.clus = "bl-M3. CD16++ CD14dim CDKN1C+")
sds <- readRDS("destiny_slingshot.Rds")
lineages <- slingLineages(sds)   # list of lineages

meta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/Myeloid_combinedMeta.csv", header = TRUE)

meta <- meta[meta$annotation %in% rightClusters, ]

counters <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/Myeloid/harmony_0/pbmcVsTissue.csv")

tissueScore = counters$tissue / (counters$tissue + counters$PBMC * (length(meta$origin[meta$origin == "Tissue"]) / 
        length(meta$origin[meta$origin == "PBMC"])))

TSM <- data.frame(barcode = meta$barcode, TSM = tissueScore)

plot_df <- data.frame(DC1 = dPCs[,1], DC2 = dPCs[,2], TSM = TSM$TSM[TSM$barcode %in% rownames(dPCs)])

plot_df$final_annotation <- meta$annotationNumber

# pdf("testDestiny.pdf", width = 8, height = 6)
#     print(ggplot(plot_df, aes(DC1, DC2, colour = final_annotation)) +
#         ggrastr::rasterise(geom_point(size = 1), dpi = 300)) + 
#         theme_minimal(base_size = 18)
# dev.off()

# pdf("testDestiny_faceted.pdf", width = 12, height = 10)
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

# pdf("testDestiny_broad.pdf", width = 8, height = 6)
#     print(ggplot(plot_df, aes(DC1, DC2, colour = type)) +
#         ggrastr::rasterise(geom_point()))
# dev.off()

# pdf("testDestiny_broad_faceted.pdf", width = 12, height = 8)
#   print(ggplot(plot_df, aes(DC1, DC2, colour = type)) +
#     ggrastr::rasterise(geom_point()) +
#     facet_wrap(~type))
# dev.off()

pdf("testDestiny_TSM.pdf", width = 8, height = 6)
    print(ggplot(plot_df, aes(DC1, DC2, colour = TSM)) +
        ggrastr::rasterise(geom_point(size = 2), dpi = 300) +
        theme_minimal(base_size = 14) + 
        scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous")))
dev.off()

# pdf("testDestiny_TSM_faceted.pdf", width = 8, height = 6)
#     print(ggplot(plot_df, aes(DC1, DC2, colour = TSM)) +
#         ggrastr::rasterise(geom_point(size = 1), dpi = 300) +
#         scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous")) +
#         facet_wrap(~type))
# dev.off()

pal <- c("#782AB6", "#C075A6", "#F7E1A0", "#FC1CBF", "#1CBE4F", "#FBE426",
    "#B10DA1", "#85660D", "#1C8356", "#C4451C", "#325A9B", "#F8A19F", rev(polychrome(26))[c(1:7,9)])

cl <- max.col(sds@elementMetadata$clusterLabels)
cluster_names <- colnames(sds@elementMetadata$clusterLabels)

# Vector of cluster IDs in desired order
desired_ids <- match(rightClusters, cluster_names)  # indices in original matrix

# remap cl to this order
cl_ordered <- match(cl, desired_ids)

png("slingshot_plot.pdf", width = 7, height = 7)
    plot(
        sds@elementMetadata$reducedDim[,1:2],
            col = pal[cl_ordered],
            pch = 16,
            cex = 0.4,
            asp = 1
        )
    lines(SlingshotDataSet(sds), lwd = 2, type = 'lineages', col = 'black')[i]
dev.off()

for (i in seq_along(lineages)) {
    
    dir.create(paste0("slingshot/slingshot_UMAPs/lineage_", i))

    # Extract pseudotime
    pt <- sds@assays@data@listData$pseudotime[, i]

    # Make a data frame for ggplot
    df <- data.frame(
        DC1 = dPCs[,1],
        DC2 = dPCs[,2],
        Pseudotime = pt,
        Cluster = plot_df$final_annotation,
        Type = plot_df$type,
        TSM = plot_df$TSM
    )

    keep <- df %>%
        group_by(Cluster) %>%
        filter(mean(is.na(Pseudotime)) <= 0.8) %>%  # keep clusters with <= 80% NA
        ungroup()
    
    keep <- unique(keep$Cluster)

    # Define a nice color palette
    pal <- colorRampPalette(brewer.pal(11,'PRGn')[-6])(100)

    pdf(paste0("slingshot/slingshot_UMAPs/lineage_", i, "/cellPlots", i, "_pseudotime.pdf"), width = 7, height = 7)
        print(ggplot(df, aes(x = DC1, y = DC2, color = Pseudotime)) +
            ggrastr::rasterise(geom_point(data = df %>% filter(!is.na(Pseudotime)), size = 0.8, na.rm = TRUE), dpi = 300) +
                scale_color_gradientn(colors = pal) +
                theme_minimal() +
                labs(title = paste("Lineage", i, "Pseudotime")) + 
            theme(legend.position = "none"))

    dev.off()

    forBoxPlot <- df %>%
        filter(Cluster %in% keep & !is.na(Pseudotime)) 

    xOrder <- forBoxPlot %>%
        group_by(Cluster) %>%
        summarize(median_pt = median(Pseudotime, na.rm = TRUE)) %>%
        arrange(median_pt) %>%
        pull(Cluster)

    print(xOrder)

    forBoxPlot$Cluster <- factor(forBoxPlot$Cluster, levels = xOrder)

    # Create Label from Cluster, keeping the order
    forBoxPlot$Label <- factor(sub("\\..*", "", forBoxPlot$Cluster), 
                            levels = sub("\\..*", "", xOrder))

    pdf(paste0("slingshot/slingshot_UMAPs/lineage_", i, "/cellPlots", i, "_lineage.pdf"), width = 7, height = 7)
            print(ggplot(data = forBoxPlot, aes(x = Label, y = Pseudotime, fill = Type)) +
                geom_boxplot() +
                theme_minimal() +
                labs(title = paste("Lineage", i, "Pseudotime by Cluster")) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none"))
    dev.off()

    xOrder <- forBoxPlot %>%
        group_by(Type) %>%
        summarize(median_pt = median(Pseudotime, na.rm = TRUE)) %>%
        arrange(median_pt) %>%
        pull(Type)

    print(xOrder)

    forBoxPlot$Type <- factor(forBoxPlot$Type, levels = xOrder)

    pdf(paste0("slingshot/slingshot_UMAPs/lineage_", i, "/collapsedCellPlots", i, "_lineage.pdf"), width = 7, height = 7)
            print(ggplot(data = forBoxPlot, aes(x = Type, y = Pseudotime, fill = Type)) +
                geom_boxplot() +
                theme_minimal() +
                labs(title = paste("Lineage", i, "Pseudotime by Cluster")) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none"))
    dev.off()

    # Correlations with TSM
    pdf(paste0("slingshot/slingshot_UMAPs/lineage_", i, "/pseudotime_TSM_correlation_lineage_", i, ".pdf"), width = 7, height = 7)
        print(ggplot(data = df, aes(x = TSM, y = Pseudotime)) +
            ggrastr::rasterise(geom_point(size = 0.5), dpi = 300) +
            geom_smooth(method = "lm") +
            theme_minimal() +
            labs(title = paste("Lineage", i, "Pseudotime vs TSM")) +
            theme(legend.position = "none"))
    dev.off()
}