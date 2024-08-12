#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")
path <- "/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis"

#-------------------------------------------------------

## FILTERING PARAMETERS
### NGENE : 500
### NUMI : 1000
### MT : 3% non-targeted

nGenepar <- 500
nUMIpar <- 1000
pctmt <- 3

#-------------------------------------------------------

#meta <- readRDS("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/full_metadata/meta_all_filtered_500nGene.rds")

#exprs_norm <- readRDS("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/seurat_output/normalized_data_all_filtered_500nGene.rds")

meta <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                            "meta_all_filtered_",
			    nGenepar,
                                 "nGene_",
                                 nUMIpar,
                                 "nUMI_",
                                 pctmt,
                                 "pctnontargetMTwdoubletandsample.rds",
                            sep = "")
)

exprs_norm <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/seurat_output/",
                            "normalized_data_all_filtered_",
                            nGenepar,
                                 "nGene_",
                                 nUMIpar,
                                 "nUMI_",
                                 pctmt,
                                 "pctnontargetMTwdoubletandsample.rds",
                            sep = "")
)


#markers <- c("DCN", "COL1A2", "LUM", "APOD", "COL1A1", "CCDC80", "FBLN1", "MFAP5", "IGFBP6", "CFD", "PTGDS", "PODXL", "HTRA1", "NPHS2", "DCN", "TPPP3",
#             "FGF1", "NPNT", "SPOCK2", "CTGF", "LST1", "TYROBP", "AIF1", "FCER1G", "LYZ", "S100A9", "CTSS", "S100A8", "FCN1", "FCGR3A", "RGS1", "HLA-DQA1",
#             "LYZ", "C1QA", "C1QB", "MS4A6A", "AIF1", "C1QC", "FCER1A", "FCER1G", "GNLY", "NKG7", "CCL4", "GZMB", "CCL3", "PRF1", "CST7", "CTSW", "KLRD1",
#             "GZMA", "CLDN5", "SERPINE2", "C10orf10", "FN1", "IGF2", "LTBP4", "PODXL", "SLC14A1", "PPP1R14A", "NPDC1", "RNASE1", "SDPR", "RAMP3", "DNASE1L3",
#             "PLVAP", "ENG", "EGFL7", "RGCC", "TMEM88", "RAMP2", "EMCN", "CRHBP", "PLAT", "RAMP2", "SOST", "FCN3", "EGFL7", "PTPRB", "ENG", "EHD3", "TAGLN",
#             "RGS5", "MYL9", "ACTA2", "TPM2", "C11orf96", "MYH11", "PLN", "RERGL", "SPARCL1", "MT1HL1", "SLC22A12", "RP11-119D9.1", "C9orf66", "PIPOX",
#             "DAO", "RNF186", "MT1A", "HAAO", "AFP", "SLC13A3", "SLC6A19", "FMO1", "PLG", "SLC34A1", "C14orf164", "ACP5", "SLC27A2", "ENPEP", "CALML3",
#             "KRT7", "SLC4A9", "ATP6V1G3", "CTB-27N1.1", "TSPAN8", "SLC26A4", "ATP6V0D2", "HEPACAM2", "CDA", "IL18", "SPINK1", "ATP6V1G3", "ATP6V0D2",
#             "SLC26A7", "TMEM101", "FAM24B", "RHCG", "SMIM6", "CTB-27N1.1", "DMRT2", "SPINK1", "ATP6V0D2", "ATP6V1G3", "DMRT2", "CTB-27N1.1", "SLC26A7",
#             "FAM24B", "SLC8A1", "KLK1", "STC1", "SLC12A1", "SFRP1", "CLDN16", "EGF", "PPP1R1A", "DUSP9", "CXCL12", "CTD-2228K2.5", "MFSD4", "PLAU", "SLC12A3",
#             "EGF", "SERPINA5", "DUSP9", "WNK4", "PVALB", "SFRP1", "RHCG", "KLK1", "SLC8A1", "RHCG", "PVALB", "SERPINA5", "TMEM178A", "TEX41", "ADM",
#             "TMPRSS2", "HMGCS2", "CD79A", "MS4A1", "IGJ", "CD79B", "VPREB3", "IRF8", "SELL", "GPR183", "BIRC3", "PLAC8", "CCL5", "CCL4", "NKG7", "GZMA",
#             "GZMH", "CST7", "GZMK", "CD7", "CTSW", "RGS1", "IL7R", "CD3E", "CD3G", "CD40LG", "LCK", "IL2RG", "TRAF3IP3", "IL7R", "CD2", "CD3E", "CD40LG",
#             "GPR171", "CD3G", "AC092580.4", "IL2RG", "IL2", "RGS1", "FXYD4", "AQP2", "STC1", "PTGER1", "GDF15", "ST6GAL1", "SMIM22", "GATA3", "SCNN1G",
#             "CALB1", "STC1", "HMGCS2", "SLC8A1", "SCNN1B", "SCNN1G", "GATA3-AS1", "RHCG", "TMEM178A", "PTGER1", "MMP7", "SLPI", "MMP7", "KRT7", "CLDN3",
#             "PAPPA2", "TSPAN8", "S100A14", "FIBIN", "CLDN16", "CTGF", "CFH", "TNNT2", "VCAM1", "BGN", "CYP1B1", "PTGDS", "SPOCK2", "SLPI", "MMP7", "VCAM1",
#             "CLDN3", "C1orf186", "ITGB8", "PRUNE2", "PROM1", "MDK", "MMP7", "KRT7", "RASSF4", "MZB1", "CLDN3", "C1orf186", "ITGB8", "PRUNE2", "PROM1",
#             "MDK", "MMP7", "KRT7", "RASSF4")

markers <- c("DCN", "COL1A2", "LUM", "APOD", "COL1A1", "CCDC80", "FBLN1", "MFAP5", "IGFBP6", "CFD", "PTGDS", "PODXL", "HTRA1", "NPHS2", "DCN", "TPPP3", "FGF1", "NPNT", "SPOCK2", "CTGF", "LST1", "TYROBP", "AIF1", "FCER1G", "LYZ", "S100A9", "CTSS", "S100A8", "FCN1", "FCGR3A", "RGS1", "HLA-DQA1", "LYZ", "C1QA", "C1QB", "MS4A6A", "AIF1", "C1QC", "FCER1A", "FCER1G", "GNLY", "NKG7", "CCL4", "GZMB", "CCL3", "PRF1", "CST7", "CTSW", "KLRD1", "GZMA", "CLDN5", "SERPINE2", "C10orf10", "FN1", "IGF2", "LTBP4", "PODXL", "SLC14A1", "PPP1R14A", "NPDC1", "RNASE1", "SDPR", "RAMP3", "DNASE1L3", "PLVAP", "ENG", "EGFL7", "RGCC", "TMEM88", "RAMP2", "EMCN", "CRHBP", "PLAT", "RAMP2", "SOST", "FCN3", "EGFL7", "PTPRB", "ENG", "EHD3", "TAGLN", "RGS5", "MYL9", "ACTA2", "TPM2", "C11orf96", "MYH11", "PLN", "RERGL", "SPARCL1", "MT1HL1", "SLC22A12", "RP11-119D9.1", "C9orf66", "PIPOX", "DAO", "RNF186", "MT1A", "HAAO", "AFP", "MT1HL1", "MT1A", "SLC13A3", "ALC6A19", "FMO1", "PLG", "SLC34A1", "C14orf164", "ACP5", "SLC27A2", "ENPEP", "CALML3", "KRT7", "SLC4A9", "ATP6V1G3", "CTB-27N1.1", "TSPAN8", "SLC26A4", "ATP6V0D2", "HEPACAM2", "CDA", "IL18", "SPINK1", "ATPV1G3", "ATP6V0D2", "SLC26A7", "TMEM101", "FAM24B", "RHCG", "SMIM6", "CTB-27N1.1", "SMRT2", "SPINK1", "ATP6V0D2", "ATP6V1G3", "DMRT2", "CTB-27N1.1", "SLC26A7", "FAM24B", "SLC8A1", "KLK1", "STC1", "SLC12A1", "SFRP1", "CLDN16", "EGF", "PPP1R1A", "DUSP9", "CXCL12", "CTD-2228K2.5", "MFSD4", "PLAU", "SLC12A3", "EGF", "SERPINA5", "DUSP9", "WNK4", "PVALB", "SFRP1", "RHCG", "KLK1", "SLC8A1", "RHCG", "PVALB", "SERPINA5", "TMEM178A", "TEX41", "ADM", "TMPRSS2", "HMGCS2", "CD79A", "MS4A1", "IGJ", "CD79B", "VPREB3", "IRF8", "SELL", "GPR183", "BIRC3", "PLAC8", "CCL5", "CCL4", "NKG7", "GZMA", "GZMH", "CST7", "GZMK", "CD7", "CTSW", "RGS1", "IL7R", "CD3E", "DC3G", "CD40LG", "LCK", "IL2RG", "TRAF3IP3", "IL7R", "CD2", "CD3E", "CD40LG", "GPR171", "CD3G", "AC092580.4", "IL2RG", "IL2", "RGS1", "FXYD4", "AQP2", "STC1", "PTGER1", "GDF15", "ST6GAL1", "SMIM22", "GATA3", "SCNN1G", "CALB1", "STC1", "HMGCS2", "SLC8A1", "SCNN1B", "SCNN1G", "GATA3-AS1", "RHCG", "TMEM178A", "PTGER1", "MMP7", "SLPI", "MMP7", "KRT7", "CLDN2", "PAPP2", "TSPAN8", "S100A14", "FIBIN", "CLDN16", "CTGF", "CFH", "TNNT2", "VCAM1", "BGN", "CYP1B1", "PTGDS", "SPOCK2", "SLPI", "MMP7", "VCAM1", "CLDN3", "C1orf186", "ITGB8", "PRUNE2", "PROM1", "MDK", "MMP7", "KRT7", "RASSF4")

#markers <- c("NPHS2", "CD3D", "DCN", "FCER1G", "IL7R")

#markers <- c("CD1C", "IRF7", "CLEC4C", "CLEC9A", "CADM1")

markers <- c("MZB1", "XBP1", "CD19")
markers <- c("GATA3", "PDGFRA", "ITGA8", "SEPT4", "COL6A1")
markers <- c("PDGFRB")
markers <- c("GATM", "GPX3", "ALDOB")
markers <- c("IGHG1", "IGHG3", "IGHG4", "IGHA")
#markers <- c("TUBA1B", "STMN1", "PTMA")

markers <- unique(markers)

markers <- markers[markers %in% rownames(exprs_norm)]

for (i in 1:length(markers)) {

  gene <- markers[i]

    max.cutoff = quantile(exprs_norm[gene,], .99)
    min.cutoff = quantile(exprs_norm[gene,], .01)

    tmp <- sapply(X = exprs_norm[gene,], FUN = function(x) {
        return(ifelse(test = x > max.cutoff, yes = max.cutoff,
            no = x))
    })
    tmp <- sapply(X = tmp, FUN = function(x) {
        return(ifelse(test = x < min.cutoff, yes = min.cutoff,
            no = x))
    })
  
    meta$gene <- as.numeric(tmp)


    ind <- ggplot(meta) +
      stat_summary_hex(aes(x = huwotUMAP1, y = huwotUMAP2, z = gene
#                          fill = cut(..value..,
#                                     c(0, 0.1, 0.5, 1, 3, Inf))
                        ),
                     fun = mean,
#                      colour = NA,
                     bins = 200, alpha = 0.75,
                     data = meta) +
#   geom_density2d(aes(x = harmonized_UMAP1, y = harmonized_UMAP2)) +
#   scale_fill_scico(palette = "berlin") +
    #scale_fill_viridis(option = "magma") +
#     scale_fill_viridis(option = "viridis", end = .9) +
    scale_fill_scico(palette = "imola") +
    guides(
      # fill = guide_colorbar(barwidth = 1, barheight = 10),
      fill = FALSE,
      color = FALSE,
      alpha = "none"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = gene
    ) +
    theme_classic(base_size = 5) +
    theme(
      plot.title = element_text(color="black", size=40, face = "italic"), # face="bold.italic"
      axis.text = element_blank(), 
      axis.ticks = element_blank(), 
      panel.grid = element_blank()
    ) 


    pdf(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/expression_plots/",
        gene,
        "_exp_all_doublet_samples.pdf",
        sep = ""),
        height = 4,
        width = 4)

    show(ind)
    dev.off()
}


rm(exprs_norm)
rm(meta)
gc()


