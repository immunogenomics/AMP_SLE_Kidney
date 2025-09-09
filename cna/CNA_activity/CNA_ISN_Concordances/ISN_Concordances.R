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

a_meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', name, '/activity/sc_meta.csv'))
a_ncorr <- read.csv(paste0("../cna_results/", name, "/Final_Activity_SiteFBChron_ncorr.csv"), header = FALSE)
a_fdrs <- read.csv(paste0("../cna_results/", name, "/Final_Activity_SiteFBChron_fdrs.csv"), header = FALSE)
colnames(a_fdrs) <- c('threshold', 'fdr', 'ncells')

a_threshold <- a_fdrs %>% filter(fdr < 0.1) %>% 
          mutate(fdr = round(fdr, 4)) %>% 
          filter(threshold == min(threshold)) %>% 
          pull(threshold)

activityCorrs <- function(comparison) {
  ncorr_2 <- read.csv(paste0("../cna_results/", name, "/Final_Activity_SiteFBChron_", comparison,"_ncorr.csv"), header = FALSE)
  
  plot_df <- data.frame(Standard = a_ncorr$V1, ISN = ncorr_2$V1)

  plot_df <- plot_df %>% mutate(a_status = case_when(
    Standard > a_threshold ~ "Expanded",
    Standard < -1 * a_threshold ~ "Depleted",
    TRUE ~ "Unchanged"
  ))

  plot_df$a_status <- factor(plot_df$a_status, levels = c("Depleted", "Unchanged", "Expanded"))

  p <- ggplot(data = plot_df, aes(x = Standard, y = ISN, color = a_status)) + 
      ggrastr::rasterise(geom_point(size = 1), dpi = 400) + 
      scale_color_manual(values = c("steelblue1", "lightgrey", "tomato")) + 
      geom_vline(xintercept = a_threshold, linetype="dotted", color = "grey") + 
      geom_vline(xintercept = -1 * a_threshold, linetype="dotted", color = "grey") + 
      annotate("text", x = Inf, y = -Inf, label = paste0("R = ", round(cor(plot_df$Standard, plot_df$ISN), 3)), hjust = 1.1, vjust = -0.5, size = 5) + 
      theme_classic() + ggtitle(comparison) + 
      theme(legend.position="none", plot.title = element_text(hjust = 0.5)) 

  pdf(paste0(name, "/", name, "_", comparison, "_corr.pdf"), width = 5, height = 5.2)
      print(p)
  dev.off()
}

activityCorrs("ISN_IV")
activityCorrs("ISN_IV_V")
activityCorrs("ISN_V")