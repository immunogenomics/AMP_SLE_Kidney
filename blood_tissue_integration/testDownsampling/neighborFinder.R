library(RANN)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly=TRUE)
name <- as.character(args[1])

harmony <- readRDS(paste0("harmony_0/", name, "_combined_hPCs.Rds"))
meta <- read.csv(paste0("PCs/", name, "_combinedMeta.csv"))

if(name != "B") {
	lookup <- bind_rows(
		read.csv(paste0("../newClusterConversion/", name, "_PBMC.csv")), 
		read.csv(paste0("../newClusterConversion/", name, "_Tissue.csv")))
} else {
	lookup <- read.csv(paste0("../newClusterConversion/", name, "_PBMC.csv"))

	bMeta <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/reorderedClusters/b_tissue_meta_reOrdered.Rds")
	meta$annotation <- ifelse(meta$origin == "Tissue", 
		bMeta$final_annotation[match(meta$barcode, bMeta$cell)], 
		meta$annotation)
}

# neighbors <- nn2(harmony, treetype = "bd", k = 50)

# saveRDS(neighbors, paste0("harmony_0/", name, "_nn.Rds"))

neighbors <- readRDS(paste0("harmony_0/", name, "_nn.Rds"))

calledCellType <- apply(neighbors$nn.idx, 1:2, function(x) meta$annotation[x])
calledCellType <- cbind(meta$annotation, calledCellType)

tissueOfOrigin <- apply(neighbors$nn.idx, 1:2, function(x, row) meta$origin[x])
tissueOfOrigin <- cbind(meta$annotation, tissueOfOrigin)

counters <- apply(tissueOfOrigin, 1, function(x){
  type <- x[1]
  PBMC <- sum(x[2:ncol(tissueOfOrigin)] == "PBMC")
  Tissue <- (50 - PBMC)
  c(type, PBMC, Tissue)
})

counters <- as.data.frame(t(counters))
colnames(counters) <- c("annotation", "PBMC", "tissue")

counters$PBMC <- as.numeric(as.character(counters$PBMC))
counters$tissue <- as.numeric(as.character(counters$tissue))

counters$TSM <- counters$tissue / (counters$tissue + counters$PBMC * (length(meta$origin[meta$origin == "Tissue"]) / 
        length(meta$origin[meta$origin == "PBMC"])))

median_count <- counters %>%
  group_by(annotation) %>%
  summarise(median_TSM = median(TSM), stdev_TSM = sd(TSM))

median_count <- median_count %>%
  mutate(annotation = ifelse(annotation %in% lookup$oldAnnotation, lookup$newAnnotation[match(annotation, lookup$oldAnnotation)], annotation))

# Original dataset
meta2 <- read.csv(paste0("../harmonyIntegration_finalFinal/PCs/", name, "_combinedMeta.csv"), header = TRUE)
counters2 <- read.csv(paste0("../harmonyIntegration_finalFinal/", name, "/harmony_0/pbmcVsTissue.csv"))

umap_table <- data.frame(
    tissueScore = counters2$tissue / (counters2$tissue + counters2$PBMC * (length(meta2$origin[meta2$origin == "Tissue"]) / 
        length(meta2$origin[meta2$origin == "PBMC"]))),
    annotation = meta2$annotation,
    annotationNum = meta2$annotationNumber,
    origin = meta2$origin
)

median_count_original <- umap_table %>%
  group_by(annotation) %>%
  summarise(median_TSM_original = median(tissueScore), stdev_TSM_original = sd(tissueScore))

plot_data <- merge(median_count, median_count_original, by = "annotation")

plot_data$origin <- ifelse(
	startsWith(plot_data$annotation, "bl-"), "PBMC", "Tissue"
)

plot_data$clusterSize <- sapply(plot_data$annotation, function(x) sum(meta2$annotation == x))
plot_data$annotationNumber <- sub("\\..*", "", plot_data$annotation)

pdf(paste0(name, "_median_TSM_comparison.pdf"))
ggplot(plot_data, aes(x = median_TSM_original, y = median_TSM, color = origin, size = clusterSize)) +
  geom_point() +   
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
#   geom_label_repel(aes(label = annotationNumber), size = 5, max.overlaps = 100, min.segment.length = 0.0000001, force = 10) +
    scale_color_manual(values = c("PBMC" = "#295177", "Tissue" = "#aa5c56")) +
  labs(x = "Median TSM (Original)", y = "Median TSM (Downsampled)", title = paste0("Spearman R: ", round(cor(plot_data$median_TSM_original, plot_data$median_TSM, method = "spearman"), 2))) +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
dev.off()

rankPlot <- plot_data %>%
  arrange(median_TSM) %>%
  mutate(annotation = factor(annotation, levels = annotation))

library(dplyr)

rank_df <- plot_data %>%
  mutate(
    rank_TSM = rank(-median_TSM, ties.method = "average"),
    rank_TSM_original = rank(-median_TSM_original, ties.method = "average")
  ) %>%
  arrange(rank_TSM)

pdf(paste0(name, "_rank_comparison.pdf"), width = 5, height = 4)
ggplot(rank_df, aes(x = rank_TSM_original, y = rank_TSM, color = origin, size = clusterSize)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
#   geom_label_repel(aes(label = annotationNumber), size = 5, max.overlaps = 100, min.segment.length = 0.0000001, force = 10) +
  scale_color_manual(values = c("PBMC" = "#295177", "Tissue" = "#aa5c56")) +
  labs(x = "Rank of Median TSM	 (Original)", y = "Rank of Median TSM (Downsampled)", title = paste0("Spearman R: ", round(cor(rank_df$rank_TSM_original, rank_df$rank_TSM, method = "spearman"), 2))) +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = c(0, max(rank_df$rank_TSM_original)), ylim = c(0, max(rank_df$rank_TSM)))
dev.off()	