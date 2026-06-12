suppressPackageStartupMessages({
	library(patchwork)
	library(purrr)
	library(dplyr) 
	library(tidyr)
	library(cowplot)
	library(singlecellmethods)
	library(ggplot2)
	library(ggbeeswarm)
	library(stringr)
	library(viridis)
    library(emmeans)
	source("/data/srlab/anathan/scripts/scseq_utils.R")
    library(progressr)
	library(Matrix)
	library(MASS)
	library(future)
	library(future.apply)
	library(tidyverse)
	library(pheatmap)
	library("RColorBrewer")
	library(googlesheets4)
	gs4_deauth()
	library(ggrepel)
})

# Load data
ptMeta <- readRDS("../../finalObjects_v3/pt_meta.Rds")
ptRaw <- readRDS('/data/srlab2/qxiao/AMP-SLE/sc_nuc_data/2022-08-09_ScNuc_cell_QCed_RawCounts.rds')

pseudobulk <- function(counts, groups) {
    colnames(counts) = groups
    
    # Convert group labels into a sparse indicator matrix
    group_levels <- unique(groups)
    group_matrix <- sparseMatrix(
      i = match(groups, group_levels), 
      j = seq_len(length(groups)), 
    )

    # Perform sparse matrix multiplication to sum within groups
    pseudo_raw <- group_matrix %*% (counts %>% t)

    # Assign row names to match groups
    rownames(pseudo_raw) <- group_levels

    pseudo_raw = pseudo_raw %>% t
    
    return(pseudo_raw)
    
}

ptRaw <- ptRaw[, match(ptMeta$cell, colnames(ptRaw))]

clinData <- readRDS("../../finalObjects_v3/clinicalData.rds")

name = 'pt'

pt0_meta <- ptMeta %>% filter(final_annotation == "PT0. Malapdative/Failed Repair")

# Add treatment response info
pt0_meta <- pt0_meta %>% 
    mutate(TreatmentResponse = clinData$Responder.Status[match(pt0_meta$sample, clinData$sample)]) %>%
    filter(TreatmentResponse %in% c("NR", "PR", "CR"))

pt0_raw <- ptRaw[, match(pt0_meta$cell, colnames(ptRaw))]

pt0_pseudo <- pseudobulk(pt0_raw, pt0_meta$sample)

# Create sample level metadata with counts
ncounts = colSums(pt0_pseudo)
sample_meta = pt0_meta[!duplicated(pt0_meta$sample), ]
rownames(sample_meta) = sample_meta$sample
sample_meta$ncounts = ncounts

ndonor_thresh = round(nrow(sample_meta)*.9)
print(ndonor_thresh)
genes = rownames(pt0_pseudo)[rowSums(pt0_pseudo >= 1) > ndonor_thresh]
print(length(genes))

pseudo_raw = pt0_pseudo[genes, rownames(sample_meta)]

gc()

system.time(
with_progress({
suppressWarnings({
    p <- progressor(along = genes)  # one step per gene
    
    # Run across all genes
    res_list <- future_lapply(genes, future.seed = TRUE, function(gene) {
  
    p()

    # Extract pseudo-bulk expression
    gene_expr <- pseudo_raw[gene, ]
    
    # Fit full model with TreatmentResponse
    full <- glm.nb(gene_expr ~ TreatmentResponse + offset(log(ncounts)), data = sample_meta)
    
    # Compute estimated marginal means
    emm <- emmeans(full, ~ TreatmentResponse)
    
    # Pairwise contrasts with explicit directions
    contr_CR_NR <- contrast(emm, method = list("CR vs NR" = c(1, -1, 0))) %>% as.data.frame()
    contr_CR_PR <- contrast(emm, method = list("CR vs PR" = c(1, 0, -1))) %>% as.data.frame()
    contr_PR_NR <- contrast(emm, method = list("PR vs NR" = c(0, -1, 1))) %>% as.data.frame()

    # Tidy each contrast
    tidy_contrast <- function(df, gene) {
      df %>%
        dplyr::select(estimate, SE, z.ratio, p.value) %>%
        dplyr::rename(beta = estimate, stderr = SE, z_score = z.ratio, p = p.value) %>%
        dplyr::mutate(gene = gene)
    }

    completeVsNone <- tidy_contrast(contr_CR_NR, gene)
    completeVsPartial <- tidy_contrast(contr_CR_PR, gene)
    partialVsNone <- tidy_contrast(contr_PR_NR, gene)

    # Return as a list of 3 dataframes
    return(list(
        completeVsNone = completeVsNone,
        completeVsPartial = completeVsPartial,
        partialVsNone = partialVsNone
    ))
    })}
)})
)

completeVsNone_all <- do.call(rbind, lapply(res_list, `[[`, "completeVsNone"))
completeVsPartial_all <- do.call(rbind, lapply(res_list, `[[`, "completeVsPartial"))
partialVsNone_all <- do.call(rbind, lapply(res_list, `[[`, "partialVsNone"))

completeVsNone_all$padj = p.adjust(completeVsNone_all$p, method = "bonferroni")
completeVsPartial_all$padj = p.adjust(completeVsPartial_all$p, method = "bonferroni")
partialVsNone_all$padj = p.adjust(partialVsNone_all$p, method = "bonferroni")

saveRDS(list(
    completeVsNone = completeVsNone_all,
    completeVsPartial = completeVsPartial_all,
    partialVsNone = partialVsNone_all
), 'pairwise_geneExpresion.rds')

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
  filename = "volcano_plot_completeVsNone.pdf",
  plot = p,
  width = 5,
  height = 5
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
  filename = "volcano_plot_completeVsPartial.pdf",
  plot = p,
  width = 5,
  height = 5
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
  filename = "volcano_plot_partialVsNone.pdf",
  plot = p,
  width = 5,
  height = 5
)
