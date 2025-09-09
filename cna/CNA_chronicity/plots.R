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

font_import(prompt = FALSE)
loadfonts()

figdir = 'cna_plots/'

args = commandArgs(trailingOnly=TRUE)

name <- as.character(args[1])
nclus <- as.integer(args[2])

xPos <- as.integer(args[3])
yPos <- as.integer(args[4])

pValues <- read.csv(paste0(name, "_pValues.csv"))

plotType <- function(coVars) {
  meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', name, '/chronicity/sc_meta.csv'))
  umap <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', name, '/chronicity/sc_umap.csv'))
  ncorr <- read.csv(paste0("cna_results/", name, "/Final_Chronicity_", coVars, "_ncorr.csv"), header = FALSE)
  fdrs <- read.csv(paste0("cna_results/", name, "/Final_Chronicity_", coVars, "_fdrs.csv"), header = FALSE)

  colnames(fdrs) <- c('threshold', 'fdr', 'ncells')
  meta$ncorr <- ncorr$V1

  fdr <- fdrs %>% filter(fdr < 0.1) %>% 
              mutate(fdr = round(fdr, 4)) %>% 
              filter(threshold == min(threshold)) %>% 
              pull(threshold)

  meta <- meta %>% mutate(ncorr_thresh = ifelse(abs(meta$ncorr) > fdr, ncorr, NA))
  meta <- cbind(meta, umap)

  if (all(c("huwotUMAP1", "huwotUMAP2") %in% colnames(meta))) {
    meta = meta %>% rename(hUMAP1 = huwotUMAP1, hUMAP2 = huwotUMAP2)
  }

  tmp_meta = meta

  globalp = as.numeric(pValues[pValues$name == coVars, ]$pValue)
  pheno = 'Chronicity'  

  title = name

  subset = c((aggregate(meta['ncorr'], meta['final_annotation'], FUN=mean) %>% arrange(ncorr) %>% tail(nclus))[[1]],
           (aggregate(meta['ncorr'], meta['final_annotation'], FUN=mean) %>% arrange(ncorr) %>% head(nclus))[[1]])
  subset = subset %>% unique

  pos_fdr_thresh <- fdr
  neg_fdr_thresh <- -1 * fdr

  vmax = quantile(abs(tmp_meta$ncorr), .95)[[1]]
  if (vmax < pos_fdr_thresh) {
    vmax <- pos_fdr_thresh + 0.03
  }
  vmin = -vmax
  custom_palette <- c("blue", "lightgray", "lightgray","lightgray",  "red") 
  break_points <- c(vmin, neg_fdr_thresh, 0, pos_fdr_thresh, vmax)  

  sig_cmap <- scale_color_gradientn(
    colours = custom_palette,
    values = rescale(break_points, to = c(0, 1)),  # Normalize thresholds
    limits = c(vmin, vmax),  # Ensure full scale
    na.value = "pink",
    guide = guide_colorbar(direction = "horizontal"),  # Horizontal legend
      breaks = c(-floor(abs(vmin)*100)/100, floor(vmax*100)/100),
  oob = scales::squish
  )

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

  ggsave(paste0(figdir, name, "/", coVars, '_UMAP_sorted.pdf'), plot = cna_umap,  
        height = 4, width = 3.5, dpi = 500, device = cairo_pdf)

  cowplot::save_plot(paste0(figdir, name, '/', coVars, '_UMAP_legend.pdf'), plot = umap_legend,  
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

  height = max(round(length(unique(meta$final_annotation_num)) / 4), 4)
  width = 4

  pdf(paste0(figdir, name, '/', coVars, "_full.pdf"), width = width, height = height)
      print(violin_plot)
  dev.off()

  meta = meta %>% mutate(final_annotation_num = final_annotation %>% str_split('. ') %>% map(1) %>% unlist)
  pos_fdr_thresh <- fdr
  neg_fdr_thresh <- -1 * fdr

  interval <- c(-max(abs(meta$ncorr)), max(abs(meta$ncorr)))

  violin_plot <- ggplot(meta %>% filter(final_annotation %in% subset), aes(y = reorder(final_annotation_num, ncorr), x = ncorr)) +
                      ggrastr::rasterise(ggbeeswarm::geom_quasirandom(aes(color = ncorr), width = 0.25, size = 0.5), dpi = 400) +
                      geom_point(data = medians %>% filter(final_annotation %in% subset), aes(x = median_nc, y = final_annotation_num)) + 
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

  height = 4
  width = 4
  outplot = violin_plot #+ plot_layout(heights = c(6, 1.5))

  ggsave(paste0(figdir, name, '/', coVars, '_vlnPlot_medians.pdf'), plot = outplot,  
        height = height, width = width, dpi = 500, device = cairo_pdf)
}

plotType("None")
plotType("SiteFirstBiopsy")
plotType("FirstBiopResponse")
plotType("SiteFirstBiopsyResponse")