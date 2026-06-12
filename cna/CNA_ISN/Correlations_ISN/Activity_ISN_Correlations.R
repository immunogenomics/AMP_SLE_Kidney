suppressPackageStartupMessages({
    library(dplyr)
    library(scales)
    library(pals)
    library(ggplot2)
    library(cowplot)
    library(tidyverse)
    library(extrafont)
})

plotCellType <- function(cellType, xlim, ylim) {
    i_corr <- read.csv(paste0("../cna_results/", cellType, "/Final_ISN_[V]_SiteFBChron_Activity_ncorr.csv"), header = FALSE)
    i_fdrs <- read.csv(paste0("../cna_results/", cellType, "/Final_ISN_[V]_SiteFBChron_Activity_fdrs.csv"), header = FALSE)
    i_meta <- read.csv(paste0("/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/", cellType, "/ISN/sc_meta.csv"))
    colnames(i_fdrs) <- c('threshold', 'fdr', 'ncells')
    
    i_corr$cell <- i_meta$cell

    a_corr <- read.csv(paste0("../cna_results/", cellType, "/Final_ISN_[V]_SiteFBChron_ncorr.csv"), header = FALSE)
    a_fdrs <- read.csv(paste0("../cna_results/", cellType, "/Final_ISN_[V]_SiteFBChron_fdrs.csv"), header = FALSE)
    a_meta <- read.csv(paste0("/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/", cellType, "/ISN/sc_meta.csv"))
    colnames(a_fdrs) <- c('threshold', 'fdr', 'ncells')

    a_corr$cell <- a_meta$cell

    cellsToKeep <- i_corr$cell[i_corr$cell %in% a_corr$cell]

    plot_df <- data.frame(Activity = i_corr$V1[i_corr$cell %in% cellsToKeep], Standard = a_corr$V1[a_corr$cell %in% cellsToKeep])
 
    a_threshold <- a_fdrs %>% filter(fdr < 0.1) %>% 
              mutate(fdr = round(fdr, 4)) %>% 
              filter(threshold == min(threshold)) %>% 
              pull(threshold)

    if (all(i_fdrs$fdr >= 0.1)) {
        i_threshold <- Inf
    } else {
        i_threshold <- i_fdrs %>% filter(fdr < 0.1) %>% 
              mutate(fdr = round(fdr, 4)) %>% 
              filter(threshold == min(threshold)) %>% 
              pull(threshold)
    }

    plot_df <- plot_df %>% mutate(a_status = case_when(
        Standard > a_threshold ~ "Expanded",
        Standard < -1 * a_threshold ~ "Depleted",
        TRUE ~ "Unchanged"
    ))

    plot_df <- plot_df %>% mutate(i_status = case_when(
        Activity > i_threshold ~ "Expanded",
        Activity < -1 * i_threshold ~ "Depleted",
        TRUE ~ "Unchanged"
    ))

    plot_df <- plot_df %>%
    mutate(color = case_when(
        a_status == "Depleted" & i_status == "Unchanged" ~ "lightblue",
        i_status == "Depleted" & a_status == "Unchanged" ~ "lightblue",
        a_status == "Depleted" & i_status == "Depleted"  ~ "darkblue",
        a_status == "Unchanged" & i_status == "Unchanged" ~ "grey",
        a_status == "Expanded" & i_status == "Unchanged" ~ "lightred",
        i_status == "Expanded" & a_status == "Unchanged" ~ "lightred",
        a_status == "Expanded" & i_status == "Expanded"  ~ "darkred",
        TRUE                                             ~ "purple"
    ))

    plot_df$color <- factor(plot_df$color, levels = c("darkblue", "lightblue", "grey", "lightred", "darkred", "purple"))

    p <- ggplot(data = plot_df, aes(x = Standard, y = Activity, color = color)) + 
        ggrastr::rasterise(geom_point(size = 1), dpi = 400) + 
        scale_color_manual(values = c("steelblue1", "lightgrey", "tomato", "purple")) + 
        geom_vline(xintercept = a_threshold, linetype="dotted", color = "grey") + 
        geom_vline(xintercept = -1 * a_threshold, linetype="dotted", color = "grey") + 
        geom_hline(yintercept = i_threshold, linetype="dotted", color = "grey") + 
        geom_hline(yintercept = -1 * i_threshold, linetype="dotted", color = "grey") + 
        annotate("text", x = Inf, y = -Inf, label = paste0("R = ", round(cor(plot_df$Standard, plot_df$Activity), 3)), hjust = 1.1, vjust = -0.5, size = 5) + 
        theme_classic() + ggtitle(cellType) + 
        # xlim(xlim) + ylim(ylim) + 
        theme(legend.position="none", plot.title = element_text(hjust = 0.5)) 

    pdf(paste0(cellType, "_ActivityCorrected_corr.pdf"), width = 5, height = 5.2)
        print(p)
    dev.off()
}

plotCellType("GLOM", c(-0.3, 0.55), c(-0.4, 0.3))
plotCellType("myeloid", c(-0.45, 0.65), c(-0.5, 0.5))