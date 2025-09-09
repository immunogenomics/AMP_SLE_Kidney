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

figdir = 'cna_percents/'
args = commandArgs(trailingOnly=TRUE)

name <- as.character(args[1])

dir.create(figdir)

percentages <- function(coVars) {
  umap <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', name, '/chronicity/sc_umap.csv'))
  meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', name, '/chronicity/sc_meta.csv'))
  ncorr <- read.csv(paste0("cna_results/", name, "/", coVars, "_ncorr.csv"), header = FALSE)
  fdrs <- read.csv(paste0("cna_results/", name, "/", coVars, "_fdrs.csv"), header = FALSE)

  colnames(fdrs) <- c('threshold', 'fdr', 'ncells')
  meta$ncorr <- ncorr$V1

  if (all(fdrs$fdr >= 0.1)) {
    print("No neighborhoods pass significance")

    result <- meta %>%
    group_by(final_annotation) %>%
    summarise(
      total = n(),
      passing = 0,
      percent = 0,
      corr = median(ncorr)
  )

  result$cellType <- name

  result <- result[, c("cellType", "final_annotation", "passing", "total", "percent", "corr")]
  colnames(result) <- c("Cell Type", "Cell State", "Passing", "Total", "Percent", "Median Correlation")

  write.csv(result, paste0(figdir, name, "_", coVars, "_percents.csv"), row.names = FALSE)
  return(NULL)
  }

  fdr <- fdrs %>% filter(fdr < 0.1) %>% 
              mutate(fdr = round(fdr, 4)) %>% 
              filter(threshold == min(threshold)) %>% 
              pull(threshold)

  result <- meta %>%
  group_by(final_annotation) %>%
  summarise(
    total = n(),
    passing = sum(ncorr > fdr),
    percent = (passing / total) * 100,
    corr = median(ncorr)
  )

  result$cellType <- name

  result <- result[, c("cellType", "final_annotation", "passing", "total", "percent", "corr")]
  colnames(result) <- c("Cell Type", "Cell State", "Passing", "Total", "Percent", "Median Correlation")

  write.csv(result, paste0(figdir, name, "_", coVars, "_percents.csv"), row.names = FALSE)
}

percentages("Final_Chronicity_SiteFirstBiopsy")