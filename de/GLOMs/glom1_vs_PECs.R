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

# genecount_thresh = 1

# celltypes = c('T/NK', 'B/Plasma', 'Myeloid', 'Glomerular', 'Interstitial/Stromal', 'Loop of Henle',
#              'Distal Nephron', 'Proximal Tubule')

# ctype_list = list()

# meta_fn = '/data/srlab/ssg34/SLE_kidney_v2/data/tissue/glom_meta_qcd_harmony_umap_clusternames_11302023.rds'
# meta_cr_fn = '/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/GLOM/chronicity/sc_meta.csv'
# raw_fn = '/data/srlab2/qxiao/AMP-SLE/sc_nuc_data/2022-08-09_ScNuc_cell_QCed_RawCounts.rds'
# name = 'glom'
# ctype_list[['Glomerular']] = c(meta_fn, raw_fn, name, meta_cr_fn)

# raw_tissue_fn = '/data/srlab2/qxiao/AMP-SLE/sc_nuc_data/2022-08-09_ScNuc_cell_QCed_RawCounts.rds'
# qcd_raw_tissue = readRDS(raw_tissue_fn)

# parallel::detectCores()
# plan(multisession, workers = 20)

# # ndonor_thresh = round(nrow(sample_meta)*.9) #round(nrow(sample_meta)/2)
# genecount_thresh = 1
# pseudocount = 0

# ctype = 'Glomerular'
# print(ctype)
# qcd_meta = read_delim(ctype_list[[ctype]][4]) %>% as.data.frame
# colnames(qcd_meta) = qcd_meta %>% colnames() %>% str_replace_all(' ', '_')
# site_dummies = qcd_meta %>% dplyr::select(starts_with('Final_Site_')) %>% colnames

# pec_ctypes = c(
#     "GLOM1. SDK1high HDAC9high PEC",
#     "GLOM0. BCAMhigh AKAP12high PEC",
#     "GLOM3. CXCL14high MMP7high PEC"
# )

# qcd_meta = qcd_meta %>% filter(final_annotation %in% pec_ctypes)

# qcd_raw = qcd_raw_tissue[, qcd_meta$cell]

# qcd_meta = qcd_meta %>% mutate(g1_ctype = ifelse(final_annotation=='GLOM1. SDK1high HDAC9high PEC', 
#                                                  'GLOM1', 'PEC'))
# qcd_meta = qcd_meta %>% mutate(sample_g1 = paste0(sample, '_', g1_ctype))

# pseudobulk <- function(counts, groups) {
#     colnames(counts) = groups
    
#     # Convert group labels into a sparse indicator matrix
#     group_levels <- unique(groups)
#     group_matrix <- sparseMatrix(
#       i = match(groups, group_levels), 
#       j = seq_len(length(groups)), 
#     )

#     # Perform sparse matrix multiplication to sum within groups
#     pseudo_raw <- group_matrix %*% (counts %>% t)

#     # Assign row names to match groups
#     rownames(pseudo_raw) <- group_levels

#     pseudo_raw = pseudo_raw %>% t
    
#     return(pseudo_raw)
    
# }

# pseudo_raw = pseudobulk(qcd_raw, qcd_meta$sample_g1)

# rm(qcd_raw)
# gc()

# # log TPM norm         
# a = 1e6
# pseudo_tpm = t(t(pseudo_raw) / colSums(pseudo_raw)) * a
# pseudo_log = log10(pseudo_tpm + 1)

# ncounts = colSums(pseudo_raw)
# sample_meta = qcd_meta[!duplicated(qcd_meta$sample_g1), ]
# rownames(sample_meta) = sample_meta$sample_g1
# sample_meta$ncounts = ncounts

# ndonor_thresh = round(nrow(sample_meta)*.9)
# print(ndonor_thresh)
# genes = rownames(pseudo_raw)[rowSums(pseudo_raw >= genecount_thresh) > ndonor_thresh]
# print(length(genes))

# pseudo_raw = pseudo_raw[genes, rownames(sample_meta)]
# pseudo_log = pseudo_log[genes, rownames(sample_meta)]

# sample_meta$g1_ctype = factor(sample_meta$g1_ctype, levels = c("PEC", "GLOM1"))

# system.time(
# with_progress({
# suppressWarnings({
# 	p <- progressor(along = genes)  # one step per gene
#     res_list = future_lapply(genes, future.seed = TRUE, function(gene){

#         tryCatch({
#     res = data.frame()

# 	p()

#     gene_expr = pseudo_log[gene, ]

#     null_fm = as.formula(paste0(c('gene_expr ~ First_biop + Final_Chronicity + Final_Activity + sample +', 
#                                   paste0(site_dummies, collapse = ' + ')), collapse = '') )
#     full_fm = as.formula(paste0(c('gene_expr ~ First_biop + Final_Chronicity + Final_Activity + sample + g1_ctype +', 
#                                   paste0(site_dummies, collapse = ' + ')), collapse = '') )

#     null = lm(null_fm, data = sample_meta)
#     full = lm(full_fm, data = sample_meta)
            
#     t =summary(full)$coefficients['g1_ctypeGLOM1', "t value"]
#     se = summary(full)$coefficients['g1_ctypeGLOM1', "Std. Error"]
#     b = summary(full)$coefficients['g1_ctypeGLOM1', "Estimate"] #full$coefficients[['TypeLN']]
#     pval = anova(full, null)[2, 'Pr(>F)']
#     res = rbind(res, data.frame('gene' = gene, 'beta' = b, 'stderr' = se, 't' = t, 'p' = pval))       
#     return(res)

#     }, error = function(e) {
#         message(paste("Error for gene:", gene))
#         message("Error message:", e$message)
#         return(NULL)  # Return NULL if there's an error to avoid stopping the loop
#       })

#     })})})
# )

# res = do.call("rbind", res_list)
# res['celltype'] = ctype
# res$padj = p.adjust(res$p, method = "bonferroni")

# saveRDS(res, "glom1_vs_pec_DE_results.rds")

# pdf("glom1_vs_pec_volcano.pdf", width = 5, height = 5)
# ggplot(res) +
# 	ggrastr::rasterise(geom_point(aes(x = beta, y = -log10(p)), color = 'grey'), dpi = 300) +	
# 	theme_classic(base_size = 14) +
# 	theme(
# 		legend.position = "none",
# 		panel.grid = element_blank(),
# 		axis.text.x = element_text(color = "black", size = 11, hjust=0.5),
# 		axis.text.y = element_text(color = "black", size = 11, hjust=1),
# 		axis.title = element_text(size=14, hjust = 0.5),
# 		plot.title = element_text(size=14, hjust = 0.5),
# 	) +
# 	geom_hline(yintercept = -log10(0.05/nrow(res)), color = 'gray55', linewidth = 1) +
# 	labs(x = 'Beta GLOM1 vs PEC', y = '-Log10(P-value)', title = "GLOM1 vs PEC")
# dev.off()	

# Run GSEA
library(msigdbr)
library(fgsea)

category = 'H'
hsig_df = msigdbr(species = "human", category = category)

category = 'C5'
gosig_df = msigdbr(species = "human", category = category)

genesig_df = rbind(hsig_df, gosig_df)
pathways = split(x = genesig_df$gene_symbol, f = genesig_df$gs_name)

betas = res %>% dplyr::select(starts_with('beta'))
ses = res %>% dplyr::select(starts_with('stderr'))
vars = ses^2
weights = 1/ vars
weighted_average <- rowSums(weights * betas, na.rm = TRUE) / rowSums(weights, na.rm = TRUE) 
weighted_se <- sqrt(1 / rowSums(weights, na.rm = TRUE))
z_scores <- weighted_average / weighted_se

names(z_scores) <- res$gene
ranks = sort(z_scores)
ranks %>% length

ranks <- sort(ranks, decreasing = TRUE)

gsea_res <- fgsea(
    pathways = pathways,
    stats    = ranks,
    minSize  = 15,
    maxSize  = 500
)
gsea_res$padj = p.adjust(gsea_res$pval, method = "bonferroni")

gsea_res_filt = gsea_res %>% filter(padj < 0.05) %>% dplyr::select(pathway, pval, padj, log2err, ES, NES) %>% arrange(desc(NES))

# Save GSEA results
saveRDS(gsea_res, "glom1_vs_pec_gsea_results.rds")

res <- readRDS("glom1_vs_pec_DE_results.rds")
# gsea_res <- readRDS("glom1_vs_pec_gsea_results.rds")

# Volcano plot with MT genes labelled
res = res %>% mutate(label_color = ifelse(gene %in% (gsea_res %>% 
    filter(pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% 
    pull(leadingEdge))[[1]], 
    'maroon', 'gray'))

# Set plotting order
res = res %>%
  mutate(plot_order = case_when(
    label_color == "gray" ~ 1,
    label_color == "maroon"  ~ 2
  )) %>%
  arrange(plot_order)
  
pdf("glom1_vs_pec_volcano_SPP1pathways.pdf", width = 5, height = 5)
ggplot(res) +
	ggrastr::rasterise(geom_point(aes(x = beta, y = -log10(p), color = label_color)), dpi = 300) +	
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
    ggrepel::geom_label_repel(data = res %>% 
        filter(label_color == 'maroon'), 
        aes(x = beta, y = -log10(p), label = gene), 
        size = 3, color = 'maroon', max.overlaps = 15, min.segment.length = 0.001, force = 10) +
    scale_color_identity() + 
	labs(x = 'Beta GLOM1 vs PEC', y = '-Log10(P-value)', title = "GLOM1 vs PEC")
dev.off()

# GSEA Select Pathway Heatmaps
gsea_pathways = c("HP_ABNORMAL_RENAL_CORPUSCLE_MORPHOLOGY", "HP_ABNORMAL_RENAL_CORTEX_MORPHOLOGY", "HP_ABNORMAL_RENAL_TUBULE_MORPHOLOGY",
                 "HP_FOCAL_SEGMENTAL_GLOMERULOSCLEROSIS", "HP_GLOMERULAR_SCLEROSIS", "HP_NEPHROPATHY", "HP_PROXIMAL_TUBULOPATHY",
                 "HP_RENAL_TUBULAR_DYSFUNCTION")

pdf("glom1_vs_pec_gsea_barplot.pdf", width = 12, height = 4)
    ggplot(gsea_res_filt %>% filter(NES > 0), aes(x = reorder(pathway, NES), y = NES, fill = -log10(padj))) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_classic(base_size = 14) +
        labs(x = "Pathway", y = "Normalized Enrichment Score (NES)", title = "GSEA Results for GLOM1 vs PEC") +
        scale_fill_viridis_c()
dev.off()