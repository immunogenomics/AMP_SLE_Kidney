suppressPackageStartupMessages({
    library(dplyr)
    library(scales)
    library(pals)
    library(ggplot2)
    library(cowplot)
    library(tidyverse)
    library(extrafont)
})

csvReader <- function(cellType) {
    ps <- read.csv(paste0(cellType, "_pValues.csv"))
    ps$cellType <- cellType
    return(ps)
}

big_ps <- csvReader("DN")
big_ps <- rbind(big_ps, csvReader("GLOM"))
big_ps <- rbind(big_ps, csvReader("INTL"))
big_ps <- rbind(big_ps, csvReader("LOH"))
big_ps <- rbind(big_ps, csvReader("PT"))
big_ps <- rbind(big_ps, csvReader("myeloid"))
big_ps <- rbind(big_ps, csvReader("b_plasma"))
big_ps <- rbind(big_ps, csvReader("t_nk"))

big_ps$name <- factor(big_ps$name, 
    levels = c(
        "None",
        "SiteFirstBiopsy",
        "SiteFBChron",
        "SiteFBChron_ISN_IV",
        "SiteFBChron_ISN_IV_V",
        "SiteFBChron_ISN_V",
        "SiteFBChron_Pred",
        "SiteFBChron_Resp",
        "SiteFBChron_Soup"
))

p <- ggplot(big_ps, aes(x = cellType, y = -log10(pValue), fill = name)) + 
    geom_bar(stat = "identity", position=position_dodge()) + 
    # scale_fill_manual(values = c("#5A189A", "#9D4EDD", "#2D6A4F", "#95D5B2")) + 
    theme_classic() + ggtitle("-log10(p) Before and After Correction") + 
    theme(plot.title = element_text(hjust = 0.5)) 

pdf("pVals.pdf", width = 12, height = 4)
    print(p)
dev.off()