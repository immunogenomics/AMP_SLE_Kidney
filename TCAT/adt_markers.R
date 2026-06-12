suppressPackageStartupMessages({
    library(singlecellmethods)
    library(tidyverse)
    library(wesanderson)
    library(viridis)
})

adt <- readRDS("/data/srlab1/Yu/tissue_blood_eQTL/AMPII_SLE/AMPII_blood/data/pheno/AllSample/adt_raw.rds")
meta <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects_v3/t_nk_pbmc_meta.Rds") 

meta$fullCodes <- paste0("channel", meta$Channel, "-", meta$BARCODE)

table(meta$fullCodes %in% colnames(adt))

meta <- meta[meta$fullCodes %in% colnames(adt), ]
adt <- adt[, colnames(adt) %in% meta$fullCodes]

meta <- meta[match(colnames(adt), meta$fullCodes), ]

adtNorm <- adt %>% singlecellmethods::normalizeData(method='cellCLR')

# Plot CD45RA and CD45RO
plot_df <- data.frame(
    hUMAP1 = meta$UMAP_1,
    hUMAP2 = meta$UMAP_2,
    CD45RA = adtNorm["PTPRC-prot-HI100", ],
    CD45RO = adtNorm["PTPRC-prot-UCHL1", ],
    CD62L = adtNorm["SELL-prot", ],
    cellType = meta$annotation
)

# Scatterplot of CD45RA vs CD45RO
pdf("CD45RA_vs_CD45RO_scatter_T10.pdf", width = 6, height = 5)
      print(ggplot(plot_df %>% filter(cellType == "T10. CD4+ Central Memory/Naive"), aes(x = CD45RA, y = CD45RO)) + 
            geom_bin_2d(bins = 100) + 
            theme_bw(base_size = 14) + 
            scale_fill_viridis(option = "turbo") +
            labs(title = "CD45RA vs CD45RO Expression", x = "CD45RA Expression", y = "CD45RO Expression"))
dev.off()

# Scatterplot of CD45RA vs CD62L
pdf("CD45RA_vs_CD62L_scatter_T10.pdf", width = 6, height = 5)
      print(ggplot(plot_df %>% filter(cellType == "T10. CD4+ Central Memory/Naive"), aes(x = CD45RA, y = CD62L)) + 
            geom_bin_2d(bins = 100) + 
            theme_bw(base_size = 14) + 
            scale_fill_viridis(option = "turbo") +
            labs(title = "CD45RA vs CD62L Expression", x = "CD45RA Expression", y = "CD62L Expression"))
dev.off()

# Scatterplot of CD45RO vs CD45RO
pdf("CD45RO_vs_CD62L_scatter_T10.pdf", width = 6, height = 5)
      print(ggplot(plot_df %>% filter(cellType == "T10. CD4+ Central Memory/Naive"), aes(x = CD45RO, y = CD62L)) + 
            geom_bin_2d(bins = 100) + 
            theme_bw(base_size = 14) + 
            scale_fill_viridis(option = "turbo") +
            labs(title = "CD45RO vs CD62L Expression", x = "CD45RO Expression", y = "CD62L Expression"))
dev.off()

# New Central Memory Cluster
# Scatterplot of CD45RA vs CD45RO
pdf("CD45RA_vs_CD45RO_scatter_T5.pdf", width = 6, height = 5)
      print(ggplot(plot_df %>% filter(cellType == "T3. CD4+ IL7Rhigh VIMhigh"), aes(x = CD45RA, y = CD45RO)) + 
            geom_bin_2d(bins = 100) + 
            theme_bw(base_size = 14) + 
            scale_fill_viridis(option = "turbo") +
            labs(title = "CD45RA vs CD45RO Expression", x = "CD45RA Expression", y = "CD45RO Expression"))
dev.off()

# Scatterplot of CD45RA vs CD62L
pdf("CD45RA_vs_CD62L_scatter_T5.pdf", width = 6, height = 5)
      print(ggplot(plot_df %>% filter(cellType == "T3. CD4+ IL7Rhigh VIMhigh"), aes(x = CD45RA, y = CD62L)) + 
            geom_bin_2d(bins = 100) + 
            theme_bw(base_size = 14) + 
            scale_fill_viridis(option = "turbo") +
            labs(title = "CD45RA vs CD62L Expression", x = "CD45RA Expression", y = "CD62L Expression"))
dev.off()

# Scatterplot of CD45RO vs CD45RO
pdf("CD45RO_vs_CD62L_scatter_T5.pdf", width = 6, height = 5)
      print(ggplot(plot_df %>% filter(cellType == "T3. CD4+ IL7Rhigh VIMhigh"), aes(x = CD45RO, y = CD62L)) + 
            geom_bin_2d(bins = 100) + 
            theme_bw(base_size = 14) + 
            scale_fill_viridis(option = "turbo") +
            labs(title = "CD45RO vs CD62L Expression", x = "CD45RO Expression", y = "CD62L Expression"))
dev.off()

# Desntiy Plot of Effector Memory Cells
pdf("CD45RA_vs_CD45RO_scatter_EM.pdf", width = 6, height = 5)
      print(ggplot(plot_df %>% filter(cellType %in% c("T9. CD4+ Effector Memory", "T13. CD4+ MAF+ IT2MA+ Effector Memory" )), aes(x = CD45RA, y = CD45RO)) + 
            geom_bin_2d(bins = 100) + 
            theme_bw(base_size = 14) + 
            scale_fill_viridis(option = "turbo") +
            labs(title = "CD45RA vs CD45RO Expression", x = "CD45RA Expression", y = "CD45RO Expression"))
dev.off()

# Scatterplot of CD45RA vs CD62L
pdf("CD45RA_vs_CD62L_scatter_EM.pdf", width = 6, height = 5)
      print(ggplot(plot_df %>% filter(cellType %in% c("T9. CD4+ Effector Memory", "T13. CD4+ MAF+ IT2MA+ Effector Memory" )), aes(x = CD45RA, y = CD62L)) + 
            geom_bin_2d(bins = 100) + 
            theme_bw(base_size = 14) + 
            scale_fill_viridis(option = "turbo") +
            labs(title = "CD45RA vs CD62L Expression", x = "CD45RA Expression", y = "CD62L Expression"))
dev.off()

# Scatterplot of CD45RO vs CD45RO
pdf("CD45RO_vs_CD62L_scatter_EM.pdf", width = 6, height = 5)
      print(ggplot(plot_df %>% filter(cellType %in% c("T9. CD4+ Effector Memory", "T13. CD4+ MAF+ IT2MA+ Effector Memory" )), aes(x = CD45RO, y = CD62L)) + 
            geom_bin_2d(bins = 100) + 
            theme_bw(base_size = 14) + 
            scale_fill_viridis(option = "turbo") +
            labs(title = "CD45RO vs CD62L Expression", x = "CD45RO Expression", y = "CD62L Expression"))
dev.off()