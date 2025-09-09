suppressPackageStartupMessages({
    library(dplyr)
    library(scales)
    library(pals)
    library(ggplot2)
    library(cowplot)
    library(tidyverse)
    library(extrafont)
})

plotCellType <- function(cellType) {
    c_ncorr <- read.csv(paste0("../cna_results/", cellType, "/Final_Chronicity_none_ncorr.csv"), header = FALSE)
    c_fdrs <- read.csv(paste0("../cna_results/", cellType, "/Final_Chronicity_none_fdrs.csv"), header = FALSE)
    colnames(c_fdrs) <- c('threshold', 'fdr', 'ncells')

    p_ncorr <- read.csv(paste0("../cna_results/", cellType, "/injured_pt_prop_none_ncorr.csv"), header = FALSE)
    p_fdrs <- read.csv(paste0("../cna_results/", cellType, "/injured_pt_prop_none_fdrs.csv"), header = FALSE)
    colnames(p_fdrs) <- c('threshold', 'fdr', 'ncells')

    plot_df <- data.frame(Chronicity = c_ncorr$V1, PT = p_ncorr$V1)
 
    p_threshold <- p_fdrs %>% filter(fdr < 0.1) %>% 
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

    plot_df <- plot_df %>% mutate(p_status = case_when(
        PT > p_threshold ~ "Expanded",
        PT < -1 * p_threshold ~ "Depleted",
        TRUE ~ "Unchanged"
    ))

    plot_df <- plot_df %>%
    mutate(color = case_when(
        p_status == "Depleted" & c_status == "Unchanged" ~ "lightblue",
        c_status == "Depleted" & p_status == "Unchanged" ~ "lightblue",
        p_status == "Depleted" & c_status == "Depleted"  ~ "darkblue",
        p_status == "Unchanged" & c_status == "Unchanged" ~ "grey",
        p_status == "Expanded" & c_status == "Unchanged" ~ "lightred",
        c_status == "Expanded" & p_status == "Unchanged" ~ "lightred",
        p_status == "Expanded" & c_status == "Expanded"  ~ "darkred",
        TRUE                                             ~ "purple"
    ))

    plot_df$color <- factor(plot_df$color, levels = c("darkblue", "lightblue", "grey", "lightred", "darkred", "purple"))

    p <- ggplot(data = plot_df, aes(x = Chronicity, y = PT, color = color)) + 
        ggrastr::rasterise(geom_point(size = 1), dpi = 400) + 
        scale_color_manual(values = c("blue", "steelblue1", "lightgrey", "tomato", "red")) + 
        geom_vline(xintercept = c_threshold, linetype="dotted", color = "grey") + 
        geom_vline(xintercept = -1 * c_threshold, linetype="dotted", color = "grey") + 
        geom_hline(yintercept = p_threshold, linetype="dotted", color = "grey") + 
        geom_hline(yintercept = -1 * p_threshold, linetype="dotted", color = "grey") + 
        annotate("text", x = Inf, y = -Inf, label = paste0("R = ", round(cor(plot_df$Chronicity, plot_df$PT), 3)), hjust = 1.1, vjust = -0.5, size = 5) + 
        theme_classic() + ggtitle(cellType) + 
        theme(legend.position="none", plot.title = element_text(hjust = 0.5)) 

    pdf(paste0(cellType, "_corr.pdf"), width = 5, height = 5.2)
        print(p)
    dev.off()
}

plotCellType("DN")
plotCellType("GLOM")
plotCellType("INTL")
plotCellType("LOH")
plotCellType("myeloid")
plotCellType("t_nk")