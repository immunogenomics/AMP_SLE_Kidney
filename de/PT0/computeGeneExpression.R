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

genecount_thresh = 1

celltypes = c('T/NK', 'B/Plasma', 'Myeloid', 'Glomerular', 'Interstitial/Stromal', 'Loop of Henle',
             'Distal Nephron', 'Proximal Tubule')

ctype_list = list()

meta_fn = '/data/srlab/ssg34/SLE_kidney_v2/data/tissue/pt_meta_qcd_harmony_umap_clusternames_11302023.rds'
meta_cr_fn = '/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/PT/chronicity/sc_meta.csv'
raw_fn = '/data/srlab2/qxiao/AMP-SLE/sc_nuc_data/2022-08-09_ScNuc_cell_QCed_RawCounts.rds'
name = 'pt'
ctype_list[['Proximal Tubule']] = c(meta_fn, raw_fn, name, meta_cr_fn)

raw_tissue_fn = '/data/srlab2/qxiao/AMP-SLE/sc_nuc_data/2022-08-09_ScNuc_cell_QCed_RawCounts.rds'
qcd_raw_tissue = readRDS(raw_tissue_fn)

parallel::detectCores()
plan(multisession, workers = 20)

# ndonor_thresh = round(nrow(sample_meta)*.9) #round(nrow(sample_meta)/2)
genecount_thresh = 1
pseudocount = 0

ctype = 'Proximal Tubule'
print(ctype)
qcd_meta = read_delim(ctype_list[[ctype]][4]) %>% as.data.frame
colnames(qcd_meta) = qcd_meta %>% colnames() %>% str_replace_all(' ', '_')
site_dummies = qcd_meta %>% dplyr::select(starts_with('Final_Site_')) %>% colnames

qcd_meta = qcd_meta %>% filter(final_annotation == "PT0. Late Injuryhigh")

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

qcd_raw = qcd_raw_tissue[, qcd_meta$cell]
pseudo_raw = pseudobulk(qcd_raw, qcd_meta$sample)

rm(qcd_raw)
gc()

ncounts = colSums(pseudo_raw)
sample_meta = qcd_meta[!duplicated(qcd_meta$sample), ]
rownames(sample_meta) = sample_meta$sample
sample_meta$ncounts = ncounts

ndonor_thresh = round(nrow(sample_meta)*.9)
print(ndonor_thresh)
genes = rownames(pseudo_raw)[rowSums(pseudo_raw >= genecount_thresh) > ndonor_thresh]
print(length(genes))

sample_meta$TreatmentResponse <- sample_meta$Responder.Status
sample_meta <- sample_meta %>% filter(TreatmentResponse %in% c("NR", "PR", "CR"))

pseudo_raw = pseudo_raw[genes, rownames(sample_meta)]

system.time(
with_progress({
suppressWarnings({
	p <- progressor(along = genes)  # one step per gene
    res_list = future_lapply(genes, future.seed = TRUE, function(gene){

        tryCatch({
    res = data.frame()

	p()

    gene_expr = pseudo_raw[gene, ]

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

    }, error = function(e) {
        message(paste("Error for gene:", gene))
        message("Error message:", e$message)
        return(NULL)  # Return NULL if there's an error to avoid stopping the loop
      })

    })})})
)

completeVsNone_all <- do.call(rbind, lapply(res_list, `[[`, "completeVsNone"))
completeVsPartial_all <- do.call(rbind, lapply(res_list, `[[`, "completeVsPartial"))
partialVsNone_all <- do.call(rbind, lapply(res_list, `[[`, "partialVsNone"))

completeVsNone_all$padj = p.adjust(completeVsNone_all$p, method = "BH")
completeVsPartial_all$padj = p.adjust(completeVsPartial_all$p, method = "BH")
partialVsNone_all$padj = p.adjust(partialVsNone_all$p, method = "BH")

saveRDS(list(
    completeVsNone = completeVsNone_all,
    completeVsPartial = completeVsPartial_all,
    partialVsNone = partialVsNone_all
), 'pairwise_geneExpresion_v2.rds')


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
  filename = "volcano_plot_completeVsNone_2.pdf",
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
  filename = "volcano_plot_completeVsPartial_2.pdf",
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
  filename = "volcano_plot_partialVsNone_2.pdf",
  plot = p,
  width = 3,
  height = 3
)