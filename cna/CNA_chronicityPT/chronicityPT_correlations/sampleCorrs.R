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
    ps <- read.csv(paste0("../ptChronicities/", cellType, "/sc_meta.csv"))
    ps <- ps %>% distinct(sample, .keep_all = TRUE)
    return(ps)
}

big_ps <- csvReader("DN")
big_ps <- rbind(big_ps, csvReader("GLOM"))
big_ps <- rbind(big_ps, csvReader("INTL"))
big_ps <- rbind(big_ps, csvReader("LOH"))
big_ps <- rbind(big_ps, csvReader("myeloid"))
big_ps <- rbind(big_ps, csvReader("t_nk"))

big_ps <- big_ps %>% distinct(sample, .keep_all = TRUE)
test <- cor.test(big_ps$injured_pt_prop, big_ps$Final_Chronicity)

p <- ggplot(data = big_ps, aes(x = Final_Chronicity, y = injured_pt_prop)) + 
    ggrastr::rasterise(geom_point(size = 1), dpi = 400) + 
    geom_smooth(method="lm", se=F) + 
    annotate("text", x = Inf, y = -Inf, 
        label = paste0("R = ", round(test$estimate, 4), "; p =", signif(test$p.value, 4)), hjust = 1.1, vjust = -0.5, size = 5) + 
    theme_classic() + xlab("Chronicity") + ylab("Late Injury High Proportion")

pdf("pt_chronicity_corr.pdf", width = 5, height = 5.2)
    print(p)
dev.off()