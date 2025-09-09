suppressPackageStartupMessages({
  library(dplyr)
  library(scales)
  library(pals)
  library(ggplot2)
  library(cowplot)
  library(tidyverse)
  library(extrafont)
})

labelfontsize = 20
tickfontsize = 16

suppressMessages({
  font_import(prompt = FALSE)
  loadfonts()
})


figdir = 'cna_plots/'

args = commandArgs(trailingOnly=TRUE)

name <- as.character(args[1])
type <- as.character(args[2])
file <- as.character(args[3])

nclus <- 4
xPos <- -Inf
yPos <- -Inf

pValues <- read.csv(paste0(name, "_", type, "_pValues.csv"))

dir.create(figdir)
dir.create(paste0(figdir, name))

if (name == "t_nk") {
  cellStates <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects/t_nk_pbmc_meta.Rds") %>% select(Cell, annotation)
  newAnnotations <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/newClusterConversion/TNK_PBMC.csv")
  cellStates <- cellStates %>%
    left_join(newAnnotations, by = c("annotation" = "oldAnnotation")) %>%
    mutate(annotation = ifelse(!is.na(newAnnotation), newAnnotation, annotation)) %>%
    select(-newAnnotation)
  xPos <- 5
  yPos <- -8
} else if (name == "mono_dc") {
  cellStates <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects/myeloid_pbmc_meta.Rds") %>% select(Cell, annotation)
  newAnnotations <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/newClusterConversion/Myeloid_PBMC.csv")
  cellStates <- cellStates %>%
    left_join(newAnnotations, by = c("annotation" = "oldAnnotation")) %>%
    mutate(annotation = ifelse(!is.na(newAnnotation), newAnnotation, annotation)) %>%
    select(-newAnnotation)
  xPos <- 9
  yPos <- -13
} else if (name == "b_cell") {
  cellStates <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects/b_pbmc_meta.Rds") %>% select(Cell, annotation)
  newAnnotations <- read.csv("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/newClusterConversion/B_PBMC.csv")
  cellStates <- cellStates %>%
    left_join(newAnnotations, by = c("annotation" = "oldAnnotation")) %>%
    mutate(annotation = ifelse(!is.na(newAnnotation), newAnnotation, annotation)) %>%
    select(-newAnnotation)
  xPos <- -2
  yPos <- -4.5
}

colnames(cellStates) <- c("Cell", "final_annotation")

plotType <- function(coVars, file) {
  meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/CNA_bloodOverlap/cna_inputs/', name, '/', type,'/meta.csv'))
  umap <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/CNA_bloodOverlap/cna_inputs/', name, '/', type,'/umap.csv'))
  ncorr <- read.csv(paste0("cna_results/", name, "/Final_", file, "_", coVars, "_ncorr.csv"), header = FALSE)
  fdrs <- read.csv(paste0("cna_results/", name, "/Final_", file, "_", coVars, "_fdrs.csv"), header = FALSE)

  colnames(fdrs) <- c('threshold', 'fdr', 'ncells')
  meta$ncorr <- ncorr$V1

  meta <- cbind(meta, umap)
  meta <- merge(meta, cellStates, by = "Cell", all.x = TRUE)


  if (all(fdrs$fdr >= 0.1)) {
    print("No neighborhoods pass significance")
    return(NULL)
  } else {
    fdr <- fdrs %>% filter(fdr < 0.1) %>% 
              mutate(fdr = round(fdr, 4)) %>% 
              filter(fdr == max(fdr)) %>% 
              pull(threshold)
  }

  meta <- meta %>% mutate(ncorr_thresh = ifelse(abs(meta$ncorr) > fdr, ncorr, NA))

  if (all(c("UMAP_1", "UMAP_2") %in% colnames(meta))) {
    meta = meta %>% rename(hUMAP1 = UMAP_1, hUMAP2 = UMAP_2)
  }

  tmp_meta = meta

  globalp = as.numeric(pValues[pValues$name == coVars, ]$pValue)
  pheno = type  

  title = name

  subset = c((aggregate(meta['ncorr'], meta['final_annotation'], FUN=mean) %>% arrange(ncorr) %>% tail(nclus))[[1]],
           (aggregate(meta['ncorr'], meta['final_annotation'], FUN=mean) %>% arrange(ncorr) %>% head(nclus))[[1]])
  subset = subset %>% unique

  pos_fdr_thresh <- fdr
  neg_fdr_thresh <- -1 * fdr


  vmax <- max(abs(tmp_meta$ncorr))

  vmin = -vmax
  custom_palette <- c("blue", "lightgray", "lightgray","lightgray",  "red") 
  break_points <- c(vmin, neg_fdr_thresh, 0, pos_fdr_thresh, vmax)  

  sig_cmap <- scale_color_gradientn(
    colours = custom_palette,
    values = rescale(break_points, to = c(0, 1)),  # Normalize thresholds
    limits = c(vmin, vmax),  # Ensure full scale
    na.value = "lightgray",
    guide = guide_colorbar(direction = "horizontal"),  # Horizontal legend
      breaks = c(-floor(abs(vmin)*100)/100, floor(vmax*100)/100),
  oob = scales::squish
  )


  tmp_meta <- tmp_meta %>% mutate(status = case_when(
    ncorr > pos_fdr_thresh ~ "Expanded",
    ncorr < neg_fdr_thresh ~ "Depleted",
    TRUE ~ "Unchanged"
  ))

  tmp_meta$status <- factor(tmp_meta$status, levels = c("Depleted", "Unchanged", "Expanded"))

  cna_umap <- ggplot() + 
        ggrastr::rasterise(geom_point(
        data = tmp_meta %>% rename(correlation = ncorr) %>% arrange(abs(correlation)), 
            aes(x = hUMAP1, y = hUMAP2, color = correlation), size = 0.75, stroke = 1e-10, shape = 20), dpi = 500) + 
        sig_cmap + 
      theme_classic(base_size = tickfontsize) +
        theme(
          legend.position = "right",
          legend.justification = c(1, 1),
          legend.text = element_text(size = labelfontsize-4),
          legend.title = element_text(size = labelfontsize),
            plot.title = element_text(hjust = 0.5, size = labelfontsize),
            axis.title = element_text(size = labelfontsize), 
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1)
      ) +
      annotate("text", x = xPos, y = yPos, 
          label = bquote(P[.(pheno)] < .(formatC(globalp, format = "e", digits = 2))), 
          size = 4) +
      labs(title = title, x = "UMAP1", y = "UMAP2", color = 'Correlation') + 
      theme(text=element_text(family="Arial Narrow")) 

  umap_legend <- cowplot::get_legend(cna_umap)
  cna_umap <- cna_umap + theme(legend.position = "none")

  ggsave(paste0(figdir, name, "/", type, "_", coVars, '_UMAP_sorted.pdf'), plot = cna_umap,  
        height = 4, width = 3.5, dpi = 500, device = cairo_pdf)

  cowplot::save_plot(paste0(figdir, name, '/', type, "_", coVars, '_UMAP_legend.pdf'), plot = umap_legend,  
        base_height = 1, base_width = 3.5, dpi = 500)

  meta = meta %>% mutate(final_annotation_num = final_annotation %>% str_split('. ') %>% map(1) %>% unlist)
  pos_fdr_thresh <- fdr
  neg_fdr_thresh <- -1 * fdr

  interval <- c(-max(abs(meta$ncorr)), max(abs(meta$ncorr)))

  medians <- meta %>%
    group_by(final_annotation) %>%
    summarise(
      median_nc = median(ncorr),
      .groups = "drop"
    )
  medians = medians %>% mutate(final_annotation_num = final_annotation %>% str_split('. ') %>% map(1) %>% unlist)

  breaks <- seq(round(-max(abs(meta$ncorr)), 2), round(max(abs(meta$ncorr)), 2), round(max(abs(meta$ncorr)), 2))

  violin_plot <- ggplot(meta, aes(y = reorder(final_annotation_num, ncorr), x = ncorr)) +
                      ggrastr::rasterise(ggbeeswarm::geom_quasirandom(aes(color = ncorr), width = 0.25, size = 0.5), dpi = 400) +
                      geom_point(data = medians, aes(x = median_nc, y = final_annotation_num)) + 
                      geom_vline(xintercept = neg_fdr_thresh, linetype = "dashed", color = "darkgrey") +
                      geom_vline(xintercept = pos_fdr_thresh, linetype = "dashed", color = "darkgrey") +
                      geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
                      sig_cmap + 
                      scale_x_continuous(breaks = breaks) +
                      labs( x= "Neighborhood Correlation", y = "", title = "") +
                      theme_classic(base_size = tickfontsize) +
                      theme(
                          text=element_text(family="Arial Narrow"),
                          legend.position = "none",
                          panel.grid = element_blank(),
                          axis.text.x = element_text(color = "black", size = labelfontsize, hjust=0.5),
                          axis.text.y = element_text(color = "black", size = labelfontsize, hjust=1),
                          axis.title = element_text(size=labelfontsize, hjust = 0.5),
                          axis.line.x.bottom = element_line(color = 'black'),
                          axis.line.y.left   = element_line(color = 'black'),
                          plot.margin = margin(10, 40, 10, 10)
                      )

  height = max(round(length(unique(meta$final_annotation_num)) / 3), 4)
  width = 4

  pdf(paste0(figdir, name, '/', coVars, "_full.pdf"), width = width, height = height)
      print(violin_plot)
  dev.off()

  # meta = meta %>% mutate(final_annotation_num = final_annotation %>% str_split('. ') %>% map(1) %>% unlist)
  # pos_fdr_thresh <- fdr
  # neg_fdr_thresh <- -1 * fdr

  # interval <- c(-max(abs(meta$ncorr)), max(abs(meta$ncorr)))

  # violin_plot <- ggplot(meta %>% filter(final_annotation %in% subset), aes(y = reorder(final_annotation_num, ncorr), x = ncorr)) +
  #                     ggrastr::rasterise(ggbeeswarm::geom_quasirandom(aes(color = ncorr), width = 0.25, size = 0.5), dpi = 400) +
  #                     geom_point(data = medians %>% filter(final_annotation %in% subset), aes(x = median_nc, y = final_annotation_num)) + 
  #                     geom_vline(xintercept = neg_fdr_thresh, linetype = "dashed", color = "darkgrey") +
  #                     geom_vline(xintercept = pos_fdr_thresh, linetype = "dashed", color = "darkgrey") +
  #                     geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  #                     sig_cmap + 
  #                     scale_x_continuous(breaks = breaks) +
  #                     labs( x= "Neighborhood Correlation", y = "", title = "") +
  #                     theme_classic(base_size = tickfontsize) +
  #                     theme(
  #                         text=element_text(family="Arial Narrow"),
  #                         legend.position = "none",
  #                         panel.grid = element_blank(),
  #                         axis.text.x = element_text(color = "black", size = labelfontsize, hjust=0.5),
  #                         axis.text.y = element_text(color = "black", size = labelfontsize, hjust=1),
  #                         axis.title = element_text(size=labelfontsize, hjust = 0.5),
  #                         axis.line.x.bottom = element_line(color = 'black'),
  #                         axis.line.y.left   = element_line(color = 'black'),
  #                         plot.margin = margin(10, 40, 10, 10)
                          
  #                         )

  # height = 4
  # width = 4
  # outplot = violin_plot #+ plot_layout(heights = c(6, 1.5))

  # ggsave(paste0(figdir, name, '/', coVars, '_vlnPlot_medians.pdf'), plot = outplot,  
  #       height = height, width = width, dpi = 500, device = cairo_pdf)
}

plotType("none", file)
plotType("tissueCorrection", file)
plotType("sexAge", file)