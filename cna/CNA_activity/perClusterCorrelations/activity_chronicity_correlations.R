suppressPackageStartupMessages({
    library(dplyr)
    library(scales)
    library(pals)
    library(ggplot2)
    library(cowplot)
    library(tidyverse)
    library(extrafont)
    library(ggrastr)
    library(janitor)
})

args = commandArgs(trailingOnly=TRUE)

name <- args[1]

dir.create(name)

c_meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', name, '/chronicity/sc_meta.csv'))
c_ncorr <- read.csv(paste0("../../CNA_chronicityCorrections_siteModified/cna_results/", name, "/Final_Chronicity_SiteFirstBiopsy_ncorr.csv"), header = FALSE)
c_fdrs <- read.csv(paste0("../../CNA_chronicityCorrections_siteModified/cna_results/", name, "/Final_Chronicity_SiteFirstBiopsy_fdrs.csv"), header = FALSE)
colnames(c_fdrs) <- c('threshold', 'fdr', 'ncells')

c_fullCorr <- data.frame(cell = c_meta$cell, Chronicity = c_ncorr$V1, final_annotation = c_meta$final_annotation)

a_meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', name, '/activity/sc_meta.csv'))
a_ncorr <- read.csv(paste0("../cna_results/", name, "/Final_Activity_SiteFBChron_ncorr.csv"), header = FALSE)
a_fdrs <- read.csv(paste0("../cna_results/", name, "/Final_Activity_SiteFBChron_fdrs.csv"), header = FALSE)
colnames(a_fdrs) <- c('threshold', 'fdr', 'ncells')

a_fullCorr <- data.frame(cell = a_meta$cell, Activity = a_ncorr$V1)

plot_df <- merge(c_fullCorr, a_fullCorr, by = "cell")

a_threshold <- a_fdrs %>% filter(fdr < 0.1) %>% 
          mutate(fdr = round(fdr, 4)) %>% 
          filter(threshold == min(threshold)) %>% 
          pull(threshold)

c_threshold <- c_fdrs %>% filter(fdr < 0.1) %>% 
          mutate(fdr = round(fdr, 4)) %>% 
          filter(threshold == min(threshold)) %>% 
          pull(threshold)

plot_df <- plot_df %>% mutate(c_status = case_when(
    Chronicity > c_threshold ~ "Expanded",
    Chronicity < -1 * c_threshold ~ "Depleted",
    TRUE ~ "Unchanged"
))

plot_df <- plot_df %>% mutate(a_status = case_when(
    Activity > a_threshold ~ "Expanded",
    Activity < -1 * a_threshold ~ "Depleted",
    TRUE ~ "Unchanged"
))

plot_df <- plot_df %>%
mutate(color = case_when(
    a_status == "Depleted" & c_status == "Unchanged" ~ "steelblue1",
    c_status == "Depleted" & a_status == "Unchanged" ~ "steelblue1",
    a_status == "Depleted" & c_status == "Depleted"  ~ "darkblue",
    a_status == "Unchanged" & c_status == "Unchanged" ~ "lightgrey",
    a_status == "Expanded" & c_status == "Unchanged" ~ "tomato",
    c_status == "Expanded" & a_status == "Unchanged" ~ "tomato",
    a_status == "Expanded" & c_status == "Expanded"  ~ "darkred",
    TRUE                                             ~ "darkorchid1"
))

plot_df$color <- factor(plot_df$color, levels = c("darkblue", "steelblue1", "lightgrey", "tomato", "darkred", "darkorchid1"))

p <- ggplot(data = plot_df, aes(x = Chronicity, y = Activity, color = color)) + 
    ggrastr::rasterise(geom_point(size = 1), dpi = 400) + 
    scale_color_manual(values = c("darkblue", "steelblue1", "lightgrey", "tomato", "darkred", "purple")) + 
    geom_vline(xintercept = c_threshold, linetype="dotted", color = "grey") + 
    geom_vline(xintercept = -1 * c_threshold, linetype="dotted", color = "grey") + 
    geom_hline(yintercept = a_threshold, linetype="dotted", color = "grey") + 
    geom_hline(yintercept = -1 * a_threshold, linetype="dotted", color = "grey") + 
    annotate("text", x = Inf, y = -Inf, label = paste0("R = ", round(cor(plot_df$Chronicity, plot_df$Activity), 3)), hjust = 1.1, vjust = -0.5, size = 5) + 
    theme_classic() + ggtitle(name) + 
    theme(legend.position="none", plot.title = element_text(hjust = 0.5)) 

pdf(paste0(name, "/", name, "_corr.pdf"), width = 5, height = 5.2)
    print(p)
dev.off()

# Violin Plot
plot_df = plot_df %>% mutate(final_annotation_num = final_annotation %>% str_split('. ') %>% map(1) %>% unlist)
c_interval <- c(-max(abs(plot_df$Chronicity)), max(abs(plot_df$Chronicity)))
a_interval <- c(-max(abs(plot_df$Activity)), max(abs(plot_df$Activity)))

medians <- plot_df %>%
  group_by(final_annotation) %>%
  summarise(
    median_act = median(Activity),
    median_chr = median(Chronicity), 
    .groups = "drop"
  )

medians = medians %>% mutate(final_annotation_num = final_annotation %>% str_split('. ') %>% map(1) %>% unlist)

a_breaks <- seq(round(-max(abs(plot_df$Activity)), 2), round(max(abs(plot_df$Activity)), 2), round(max(abs(plot_df$Activity)), 2))
c_breaks <- seq(round(-max(abs(plot_df$Chronicity)), 2), round(max(abs(plot_df$Chronicity)), 2), round(max(abs(plot_df$Chronicity)), 2))

# Create color scales
# Activity
vmax = quantile(abs(plot_df$Activity), .95)[[1]]
if (vmax < a_threshold) {
  vmax <- a_threshold + 0.03
}
vmin = -vmax
custom_palette <- c("blue", "lightgray", "lightgray","lightgray",  "red") 
break_points <- c(vmin, -1 * a_threshold, 0, a_threshold, vmax)  

a_sig_cmap <- scale_color_gradientn(
  colours = custom_palette,
  values = rescale(break_points, to = c(0, 1)),  # Normalize thresholds
  limits = c(vmin, vmax),  # Ensure full scale
  na.value = "pink",
  guide = guide_colorbar(direction = "horizontal"),  # Horizontal legend
    breaks = c(-floor(abs(vmin)*100)/100, floor(vmax*100)/100),
oob = scales::squish
)

# Chronicity
vmax = quantile(abs(plot_df$Chronicity), .95)[[1]]
if (vmax < c_threshold) {
  vmax <- c_threshold + 0.03
}
vmin = -vmax
custom_palette <- c("blue", "lightgray", "lightgray","lightgray",  "red") 
break_points <- c(vmin, -1 * c_threshold, 0, c_threshold, vmax)  

c_sig_cmap <- scale_color_gradientn(
  colours = custom_palette,
  values = rescale(break_points, to = c(0, 1)),  # Normalize thresholds
  limits = c(vmin, vmax),  # Ensure full scale
  na.value = "pink",
  guide = guide_colorbar(direction = "horizontal"),  # Horizontal legend
    breaks = c(-floor(abs(vmin)*100)/100, floor(vmax*100)/100),
oob = scales::squish
)

labelfontsize = 20
tickfontsize = 16

a_vlnPlot <- ggplot(plot_df, aes(y = reorder(final_annotation_num, Activity), x = Activity)) +
                    ggrastr::rasterise(ggbeeswarm::geom_quasirandom(aes(color = Activity), width = 0.25, size = 0.5), dpi = 400) +
                    geom_point(data = medians, aes(x = median_act, y = final_annotation_num)) + 
                    geom_vline(xintercept = -1 * a_threshold, linetype = "dashed", color = "darkgrey") +
                    geom_vline(xintercept = a_threshold, linetype = "dashed", color = "darkgrey") +
                    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
                    a_sig_cmap + 
                    scale_x_continuous(breaks = a_breaks) +
                    labs( x= "Neighborhood Correlation", y = "", title = "") +
                    theme_classic(base_size = tickfontsize) +
                    ggtitle("Activity") + 
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

c_vlnPlot <- ggplot(plot_df, aes(y = reorder(final_annotation_num, Activity), x = Chronicity)) +
                    ggrastr::rasterise(ggbeeswarm::geom_quasirandom(aes(color = Chronicity), width = 0.25, size = 0.5), dpi = 400) +
                    geom_point(data = medians, aes(x = median_chr, y = final_annotation_num)) + 
                    geom_vline(xintercept = -1 * c_threshold, linetype = "dashed", color = "darkgrey") +
                    geom_vline(xintercept = c_threshold, linetype = "dashed", color = "darkgrey") +
                    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
                    a_sig_cmap + 
                    scale_x_continuous(breaks = c_breaks) +
                    labs( x= "Neighborhood Correlation", y = "", title = "") +
                    theme_classic(base_size = tickfontsize) +
                    ggtitle("Chronicity") + 
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

plots <- ggpubr::ggarrange(a_vlnPlot, c_vlnPlot, ncol = 2)

height = max(round(length(unique(plot_df$final_annotation_num)) / 3), 4)

pdf(paste0(name, "/", name, "_combinedVlns.pdf"), width = 10, height = height)
  print(plots)
dev.off()

smallPlot <- function(cluster) {
  small_df <- plot_df[plot_df$final_annotation == cluster, ]
  
  small_df$color <- factor(small_df$color, levels = unique(small_df$color))

  p <- ggplot(data = small_df, aes(x = Chronicity, y = Activity, color = color)) + 
    ggrastr::rasterise(geom_point(size = 1), dpi = 400) + 
    scale_color_manual(values = as.vector(levels(small_df$color))) + 
    geom_vline(xintercept = c_threshold, linetype="dotted", color = "grey") + 
    geom_vline(xintercept = -1 * c_threshold, linetype="dotted", color = "grey") + 
    geom_hline(yintercept = a_threshold, linetype="dotted", color = "grey") + 
    geom_hline(yintercept = -1 * a_threshold, linetype="dotted", color = "grey") + 
    annotate("text", x = Inf, y = -Inf, label = paste0("R = ", round(cor(small_df$Chronicity, small_df$Activity), 3)), hjust = 1.1, vjust = -0.5, size = 5) + 
    theme_classic() + ggtitle(cluster) + 
    theme(legend.position="none", plot.title = element_text(hjust = 0.5)) 

    pdf(paste0(name, "/", make_clean_names(cluster), "_corr.pdf"), height = 5, width = 5)
      print(p)
    dev.off()
}

for (cluster in unique(plot_df$final_annotation)) {
  smallPlot(cluster)
}