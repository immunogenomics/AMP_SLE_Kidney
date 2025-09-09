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

percents <- data.frame(
  type = character(),  # Empty character column
  percent = numeric()   # Empty numeric column
)

plotType <- function(coVars) {
  meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', name, '/chronicity/sc_meta.csv'))
  umap <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', name, '/chronicity/sc_umap.csv'))
  ncorr <- read.csv(paste0("cna_results/", name, "/Final_Chronicity_", coVars, "_ncorr.csv"), header = FALSE)
  fdrs <- read.csv(paste0("cna_results/", name, "/Final_Chronicity_", coVars, "_fdrs.csv"), header = FALSE)

  colnames(fdrs) <- c('threshold', 'fdr', 'ncells')
  meta$ncorr <- ncorr$V1

  if (all(fdrs$fdr >= 0.1)) {
    percents <<- rbind(percents, data.frame(type = coVars, percent = 0))
    print("No neighborhoods pass significance")
    return(NULL)
  }

  fdr <- fdrs %>% filter(fdr < 0.1) %>% 
              mutate(fdr = round(fdr, 4)) %>% 
              filter(threshold == min(threshold)) %>% 
              pull(threshold)

  percents <<- rbind(percents, data.frame(type = coVars, percent = sum(meta$ncorr > fdr) / nrow(meta)))
}

plotType("None")
plotType("SiteFirstBiopsy")
plotType("SiteFirstBiopsyResponse")

dir.create("cna_percents_broad")
write.csv(percents, paste0("cna_percents_broad/", name, "_percents.csv"), row.names = FALSE)
