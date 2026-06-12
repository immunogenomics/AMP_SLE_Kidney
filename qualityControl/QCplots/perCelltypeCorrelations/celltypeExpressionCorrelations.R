# Read in single cell data
suppressPackageStartupMessages({
    library(tidyverse)
    library(Matrix)
    library(ggplot2)
    library(scales)
    library(ggrepel)
})

# Read in single cell
myeloid_meta <- readRDS("../../finalObjects_v3/myeloid_tissue_meta.Rds")  %>% select(cell, sample, dataset, nCount_RNA, nFeature_RNA) %>% mutate(label = "Myeloid")
b_meta <- readRDS("../../finalObjects_v3/b_tissue_meta.Rds") %>% select(cell, sample, dataset, nCount_RNA, nFeature_RNA) %>% mutate(label = "Bcell")
tnk_meta <- readRDS("../../finalObjects_v3/t_nk_tissue_meta.Rds")  %>% select(cell, sample, dataset, nCount_RNA, nFeature_RNA) %>% mutate(label = "TNK")
glom_meta <- readRDS("../../finalObjects_v3/glom_meta.Rds") %>% select(cell, sample, dataset, nCount_RNA, nFeature_RNA) %>% mutate(label = "Glom")
loh_meta <- readRDS("../../finalObjects_v3/loh_meta.Rds") %>% select(cell, sample, dataset, nCount_RNA, nFeature_RNA) %>% mutate(label = "LOH")
pt_meta <- readRDS("../../finalObjects_v3/pt_meta.Rds") %>% select(cell, sample, dataset, nCount_RNA, nFeature_RNA) %>% mutate(label = "PT")
intl_meta <- readRDS("../../finalObjects_v3/intl_meta.Rds") %>% select(cell, sample, dataset, nCount_RNA, nFeature_RNA) %>% mutate(label = "INTL")
dn_meta <- readRDS("../../finalObjects_v3/dn_meta.Rds") %>% select(cell, sample, dataset, nCount_RNA, nFeature_RNA) %>% mutate(label = "DN")

# Read in single cell
raw_tissue_fn = '/data/srlab2/qxiao/AMP-SLE/sc_nuc_data/2022-08-09_ScNuc_cell_QCed_RawCounts.rds'
raw_tissue = readRDS(raw_tissue_fn)

dir.create("figures/")

globalCorrs <- function(meta) {
	print(paste0("Processing ", meta$label[1]))
	sc_cells <- meta$cell[meta$dataset == "scRNAseq"]
	sn_cells <- meta$cell[meta$dataset == "snRNAseq"]

	sc_expr <- raw_tissue[, sc_cells]
	sn_expr <- raw_tissue[, sn_cells]

	sc_cpm <- Matrix::rowSums(sc_expr) / sum(sc_expr) * 1e6
	sn_cpm <- Matrix::rowSums(sn_expr) / sum(sn_expr) * 1e6

	print(cor(sc_cpm, sn_cpm, method = "spearman"))
	pdf(paste0("figures/", meta$label[1], "_sc_sn_correlation_density.pdf"), width = 4, height = 4)
	print(
		ggplot(data.frame(sc = log2(sc_cpm), sn = log2(sn_cpm)), aes(x = sc, y = sn)) +
			geom_bin2d(bins = 100) +
			scale_fill_viridis_c(trans = "log10") +
			geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
			labs(
			title = paste0(meta$label[1], "; Spearman R = ", round(cor(sc_cpm, sn_cpm, method = "spearman"), 2)),
			x = "log2(scRNA-seq CPM)",
			y = "log2(snRNA-seq CPM)"
			) +
			theme_minimal(base_size = 14)
		)
	dev.off()
}

lapply(list(myeloid_meta, b_meta, tnk_meta, glom_meta, loh_meta, pt_meta, intl_meta, dn_meta), globalCorrs)