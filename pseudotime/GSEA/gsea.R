library(fgsea)
library(msigdbr)
library(tidyverse)

category = 'H'
genesig_df = msigdbr(species = "human", category = category)
pathways = split(x = genesig_df$gene_symbol, f = genesig_df$gs_name)

# Make and plot GSEA
gseaPlotMaker <- function(base_dir, name = "slingshot") {
  # Path to slingshot directory
  slingshot_dir <- file.path(base_dir, name)
	print(slingshot_dir)
  # Get lineage folders
  lineage_dirs <- list.dirs(
    path = slingshot_dir,
    full.names = TRUE,
    recursive = FALSE
  )
  
  # Keep only folders that match lineage_*
  lineage_dirs <- lineage_dirs[grepl("lineage_", basename(lineage_dirs))]
  
  for (ld in lineage_dirs) {
    
    lineage_name <- basename(ld)
    message("Running GSEA for ", lineage_name)
    
    # Example: read in correlation results
    res <- read.csv(file.path(ld, "correlatedGenes.csv"))
    
    # Create ranked gene list
    gene_list <- res$cor
    names(gene_list) <- res$gene
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    # Run preranked GSEA (example using fgsea)
    gsea_res <- fgsea::fgsea(
      pathways = pathways,
      stats = gene_list,
      nperm = 10000
    )
    
	gsea_res$padj <- p.adjust(gsea_res$pval, method = "BH")
	gsea_res_filt <- gsea_res %>%
	  filter(padj < 0.05) %>%
	  arrange(desc(NES))

    # Save results
    saveRDS(gsea_res, file.path(ld, "gsea_results.rds"))
    saveRDS(gsea_res_filt, file.path(ld, "gsea_results_filtered.rds"))
    

	pdf(file.path(ld, "gsea_results_barplot.pdf"), width = 12, height = 4)
    print(ggplot(gsea_res_filt %>% filter(NES > 0), aes(x = reorder(pathway, NES), y = NES, fill = -log10(padj))) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_classic(base_size = 14) +
        labs(x = "Pathway", y = "Normalized Enrichment Score (NES)", title = paste("GSEA Results for", lineage_name)) +
        scale_fill_viridis_c())
	dev.off()
  }
}

# gseaPlotMaker("../M0_PT_correlations/")
# gseaPlotMaker("../B_cells/")
gseaPlotMaker("../T_cells_noProlif/splitCD4and8/", "slingshot_CD4")
# gseaPlotMaker("../T_cells_v2/splitCD4and8/", "slingshot_CD8")