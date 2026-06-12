suppressPackageStartupMessages({
    library(pheatmap)
    library(RColorBrewer)
    library(ggplot2)
    library(dplyr)
    library(tidyverse)
    library(wesanderson)
    library(viridis)
})

meta <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects_v3/t_nk_pbmc_meta.Rds")

scores <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/TCAT/tcat_scores.csv")
usages <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/TCAT/tcat_usage.csv")

normalize_channel <- function(x) {
  x <- gsub("channel_", "channel", x)
  x <- gsub("_", "-", x)
  x <- gsub("channel([0-9]+)-", "channel\\1-", x)
  x
}

meta$normalized_channels <- normalize_channel(meta$Cell)
meta <- meta[meta$normalized_channels %in% scores$gene_symbol, ]

meta <- meta %>% mutate(annotation = paste0("bl-", annotation))

usages <- usages[match(meta$normalized_channels, usages$gene_symbol), ]
scores <- scores[match(meta$normalized_channels, scores$gene_symbol), ]

plot_df <- data.frame(
    hUMAP1=meta$UMAP_1,
    hUMAP2=meta$UMAP_2,
    CD4.Naive=usages$CD4.Naive,
    CD4.CM=usages$CD4.CM
)

plot_umap_continuous <- function(df, feature, title = NULL) {
    df <- df[order(df[[feature]]), ]

    p <- ggplot(df, aes(x = hUMAP1, y = hUMAP2, color = !!sym(feature))) +
        ggrastr::rasterise(geom_point(size = 0.5), dpi = 300) +
        scale_color_viridis(option = "viridis") +
        theme_minimal() +
        labs(color = feature, title = title) +
        theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    )
    return(p)
}

pdf("TCAT_CD4_Naive_UMAP.pdf", width = 6, height = 5)
    print(plot_umap_continuous(plot_df, "CD4.Naive", "TCAT CD4 Naive Usage"))
dev.off()

pdf("TCAT_CD4_CM_UMAP.pdf", width = 6, height = 5)
    print(plot_umap_continuous(plot_df, "CD4.CM", "TCAT CD4 Central Memory Usage"))
dev.off()

# Create heatmap of TCAT scores
## ---- NK remap (by original NK number) ----
labels <- unique(meta$annotation)
labels_new <- labels

nk_idx <- grep("^bl-NK", labels)
nk_old <- as.integer(sub("^bl-NK([0-9]+).*", "\\1", labels[nk_idx]))

# order by original NK number
nk_order <- order(nk_old)

# create mapping: old -> new
nk_map <- setNames(
  paste0("bl-NK", seq_along(nk_old) - 1),
  nk_old[nk_order]
)

# element-wise substitution (THIS is the key fix)
labels_new[nk_idx] <- mapply(
  function(lbl, old) sub("^bl-NK[0-9]+", nk_map[as.character(old)], lbl),
  labels[nk_idx],
  nk_old,
  USE.NAMES = FALSE
)

## ---- T remap (by original T number) ----
t_idx <- grep("^bl-T", labels_new)
t_old <- as.integer(sub("^bl-T([0-9]+).*", "\\1", labels_new[t_idx]))

# order by original T number
t_order <- order(t_old)

# map old T numbers -> new T numbers (start at T3)
t_map <- setNames(
  paste0("bl-T", seq_along(t_old) + 2),
  t_old[t_order]
)

# element-wise substitution (critical)
labels_new[t_idx] <- mapply(
  function(lbl, old) sub("^bl-T[0-9]+", t_map[as.character(old)], lbl),
  labels_new[t_idx],
  t_old,
  USE.NAMES = FALSE
)

# named lookup: old -> new
label_map <- setNames(labels_new, labels)

# apply to metadata
meta$annotation_new <- label_map[meta$annotation]

usages$annotation <- meta$annotation_new

# Create heatmap of each cluster by its TCAT usage
program_cols <- setdiff(
  colnames(usages),
  c("gene_symbol", "annotation")
)

# Compute mean program usage per annotation
program_means <- usages %>%
  filter(!grepl("^bl-NK", annotation)) %>%  # optional: remove NK
  group_by(annotation) %>%
  summarise(across(all_of(program_cols), mean, na.rm = TRUE)) %>%
  ungroup()

program_means_z <- program_means

program_means_z[program_cols] <- scale(program_means[program_cols])

mat <- t(as.matrix(program_means_z[, program_cols]))
colnames(mat) <- program_means_z$annotation

program_order <- c(
  # 1. T / NK cell subsets
  "CD4.Naive",
  "CD4.CM",
  "CD8.Naive",
  "CD8.EM",
  "CD8.Trm",
  "TEMRA",
  "Cytotoxic",
  "gdT",
  "MAIT",
  "Treg",

  # 2. Cell cycle / proliferation
  "ISG",
  "CellCycle.S",
  "CellCycle.Late.S",
  "CellCycle.G2M",
  "IEG",
  "IEG2",
  "IEG3",
  "Poor.Quality",

  # 3. Metabolism / housekeeping
  "Translation",
  "Mito",
  "RGCC.MYADM",
  "BCL2.FAM13A",

  # 4. Stress / protein handling
  "Heatshock",
  "Metallothionein",
  "Cytoskeleton",

  # 5. Myeloid / antigen presentation / doublets
  "HLA",
  "CD172a.MERTK",
  "Doublet.Myeloid",
  "Doublet.RBC",
  "Doublet.Platelet",
  "Doublet.Bcell",
  "Doublet.Plasmablast",
  "Doublet.Fibroblast",

  # 6. Cytokines / T helper programs
  "Th1.Like",
  "Th2.Activated",
  "Th2.Resting",
  "Th17.Activated",
  "Th17.Resting",
  "Th22",
  "Tfh.1",
  "Tfh.2",
  "Multi.Cytokine",
  "IL10.IL19",
  "OX40.EBI3",
  "ICOS.CD38",
  "CD40LG.TXNIP",
  "Tph",
  "Exhaustion",

  # 7. NK / other
  "CTLA4.CD38",
  "NME1.FABP5",
  "TIMD4.TIM3",
  "SOX4.TOX2"
)

cellOrder <- c(
    "bl-T8. CD8+ Central Memory/Naive",
    "bl-T4. CD8+ GZMK+ CD74+ HLA-DR+",
    "bl-T14. CD8+ GZMK+ CD74high HLA-DRhigh",
    "bl-T7. CD8+ GZMK+ TEMRA",
    "bl-T15. CD8+ GZMB+ DNMT1+ HELLS+ Proliferating",
    "bl-T18. CD8+ GZMB+ PCNAhigh Proliferating",
    "bl-T19. CD8+ GZMB+ CENPFhigh Proliferating",
    "bl-T6. CD8+ MT-high",
    "bl-T3. CD8+ GZMHhigh FGFBP2high",
    "bl-T5. CD4+ IL7Rhigh VIMhigh",
    "bl-T16. CD4+ T-reg",
    "bl-T10. CD4+ Central Memory/Naive",
    "bl-T13. CD4+ MAF+ IT2MA+ Effector Memory",
    "bl-T9. CD4+ Effector Memory",
    "bl-T11. TRDC+ Gamma/Delta",
    "bl-T12. TRGC1+ Gamma/Delta",
    "bl-T17. ISGhigh"
)

mat <- mat[rev(program_order), cellOrder]
pdf("TCAT_program_enrichment_heatmap.pdf", width = 10, height = 15)
    pheatmap(mat,
            color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
            cluster_rows = FALSE,      # cluster programs
            cluster_cols = FALSE,      # cluster annotations
             treeheight_row = 0,           # hide row dendrogram
            treeheight_col = 0,           # hide column dendrogram
            scale = "none",
            fontsize_row = 18,
            fontsize_col = 18,
            angle_col = 315,           # rotate column labels 45°
            main = "Program Enrichment by Cell Type (Z-scores)")
dev.off()

# Odds ratio with multinomial classifier
toAndFro <- data.frame(originalAnno = meta$annotation_new,
                       newAnno = scores$Multinomial_Label) %>% filter(!grepl("^bl-NK", originalAnno))

# 2️⃣ Wide frequency table: rows = clusters, columns = annotations
freq_table <- as.data.frame.matrix(table(toAndFro$newAnno, toAndFro$originalAnno))

# # Optional: add small pseudocount to avoid zero counts
# freq_table <- freq_table + 0.5

# 3️⃣ Long table for reference (optional)
freq_table_long <- freq_table %>%
  tibble::rownames_to_column("cell_type_pred_knn") %>%
  pivot_longer(cols = -cell_type_pred_knn, names_to = "annotation", values_to = "Freq")

# 4️⃣ Compute ORs and chi-square p-values per cluster vs rest
cluster_names <- rownames(freq_table)
annotations <- colnames(freq_table)

or_stats <- do.call(rbind, lapply(cluster_names, function(cluster) {
  cluster_counts <- as.numeric(freq_table[cluster, ])
  other_counts <- colSums(freq_table) - cluster_counts
  
  do.call(rbind, lapply(1:length(annotations), function(i) {
    # Build 2x2 table
    contin <- matrix(c(cluster_counts[i], other_counts[i],
                       sum(cluster_counts[-i]), sum(other_counts[-i])),
                     nrow = 2)
    
    # Skip degenerate tables
    if (any(rowSums(contin) == 0) || any(colSums(contin) == 0)) {
      return(data.frame(
        cluster = cluster,
        annotation = annotations[i],
        OR = NA,
        pvalue = NA,
        log_OR = NA
      ))
    }
    
    # Chi-square test
    chisq_res <- chisq.test(contin, simulate.p.value = TRUE, B = 100000)
    
    # Odds ratio
    OR <- (contin[1,1]*contin[2,2]) / (contin[1,2]*contin[2,1])
    
    data.frame(
      cluster = cluster,
      annotation = annotations[i],
      OR = OR,
      pvalue = chisq_res$p.value,
      log_OR = log(OR)
    )
  }))
}))

# Optional: replace infinite log_OR with NA to avoid plotting issues
or_stats <- or_stats %>%
  mutate(log_OR = ifelse(is.infinite(log_OR), NA, log_OR))

# Heatmap
pdf("T_Cluster_Enrichment_Heatmap_truncated.pdf", width=12, height=6)
    print(ggplot() +
    geom_tile(data = or_stats %>% filter(cluster %in% c("CD4_CM", "CD4_Naive")), aes(x = annotation, y = cluster, fill = log_OR), color = "grey80") +                # grey borders for tiles
    geom_text(data = freq_table_long %>% filter(cell_type_pred_knn %in% c("CD4_CM", "CD4_Naive")), aes(x = annotation, y = cell_type_pred_knn, label = Freq), size = 3) +  # add counts
    scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", # blue = depletion, red = enrichment
        midpoint = 0,                              # white = OR ~ 1
        na.value = "grey90"                        # color for NA values
    ) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank()
    ) +
    labs(
        x = "Annotation",
        y = "Cluster",
        fill = "log(OR)",
        title = "T Cluster Enrichment Heatmap"
    ))
dev.off()
