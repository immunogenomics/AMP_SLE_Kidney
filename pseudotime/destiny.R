library(destiny)
library(tidyverse)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
name <- as.character(args[1])

# Test with normalized counts

# Test with hPCs
hPCs <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/Myeloid/harmony_0/Myeloid_combined_hPCs.Rds")
umap <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/Myeloid/harmony_0/Myeloid_umap.Rds")
meta <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/harmonyIntegration_finalFinal/PCs/Myeloid_combinedMeta.csv")

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

hPCs <- as.matrix(hPCs[rownames(hPCs) %in% meta$barcode ,])

dm <- DiffusionMap(hPCs)

plot_df <- as.data.frame(dm@eigenvectors)
plot_df$final_annotation <- meta$annotationNumber

pdf("testDestiny.pdf", width = 8, height = 6)
    print(ggplot(plot_df, aes(DC1, DC2, colour = final_annotation)) +
        ggrastr::rasterise(geom_point()))
dev.off()

plot_df$type <- case_when(
    grepl("bl-", plot_df$final_annotation) ~ "Blood Monocyte",
    plot_df$final_annotation %in% c("M0", "M1", "M2", "M3") ~ "Tissue Monocyte",
    plot_df$final_annotation %in% c("M5", "M7", "M9", "M11") ~ "SAM",
    plot_df$final_annotation %in% c("M10", "M12") ~ "Resident Macrophage",
    TRUE ~ "Tissue Macrophage"
)

pdf("testDestiny_broad.pdf", width = 8, height = 6)
    print(ggplot(plot_df, aes(DC1, DC2, colour = type)) +
        ggrastr::rasterise(geom_point()))
dev.off()

pdf("testDestiny_broad_faceted.pdf", width = 12, height = 8)
  print(ggplot(plot_df, aes(DC1, DC2, colour = type)) +
    ggrastr::rasterise(geom_point()) +
    facet_wrap(~type))
dev.off()