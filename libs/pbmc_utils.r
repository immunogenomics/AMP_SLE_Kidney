set.seed(0)
library(Seurat)
library(dplyr) 
library(tidyr)
library(viridis)
library(stringr)
library(pheatmap)
library(presto)
library(pals)
library(harmony)
library(cowplot)
library(singlecellmethods)  
library(ggplot2)
#library(Cairo)
source("/data/srlab/anathan/scripts/scseq_utils.R")
library(parallel)
.libPaths("/PHShome/ssg34/.conda/envs/plswork/lib/R/library")

plot_shuffled_features <- function(input_df, input_norm, feature, pct, pt_size, legend_break) {
    
max.cutoff = quantile(input_norm[feature, input_df$cell], pct)
min.cutoff = quantile(input_norm[feature, input_df$cell], 1-pct)

plot_df <- input_df 
    
plot_df$norm_expression <- input_norm[feature, input_df$Cell]
    
plot_df <- plot_df %>% 
                mutate(norm_expression = ifelse(norm_expression < min.cutoff, min.cutoff, norm_expression)) %>%
                mutate(norm_expression = ifelse(norm_expression > max.cutoff, max.cutoff, norm_expression)) 

legend_breaks <- c(seq(0, max(plot_df$norm_expression), legend_break), round(max(plot_df$norm_expression), 1))

p <- ggplot() +
        geom_point(
            data = plot_df[sample(nrow(plot_df)), ] %>% 
                      select(UMAP_1, UMAP_2, norm_expression), 
            aes(x = UMAP_1, y = UMAP_2, color = norm_expression),
            size = pt_size, stroke = 0.0001, shape = 20) +
      scale_color_viridis(name = "Expression", breaks = legend_breaks) +
      labs(x="", y="", title = feature) + 
    theme_classic(base_size = 30) + 
    theme(plot.title = element_text(size = 30, hjust = 0.5, face = "bold.italic"),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20),
          axis.title = element_text(hjust = 0.75, 
                                    size = 20, 
                                    face = "bold"), 
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    labs(x = "UMAP1", y = "UMAP2")
return(p)
}

pivot_wilcox_expr <- function(wilcox_obj) {
    df <- wilcox_obj %>% select(group, feature, avgExpr) %>% 
                    pivot_wider(names_from = 'feature', values_from = 'avgExpr') %>% data.frame()
    rownames(df) <- df$group
    df <- df %>% select(-group) %>% as.matrix()
    return(df)
}

pivot_wilcox_logFC <- function(wilcox_obj) {
    df <- wilcox_obj %>% select(group, feature, logFC) %>% 
                    pivot_wider(names_from = 'feature', values_from = 'logFC') %>% data.frame()
    rownames(df) <- df$group
    df <- df %>% select(-group) %>% as.matrix()
    return(df)
}

pseudobulk <- function(individual, meta, norm) {
    cells <- meta %>% filter(sample == individual) %>% pull(cell)
    if (length(cells) > 2) {
        pb <- rowMeans(norm[, cells])
        pb <- c(pb, individual)
        names(pb) <- c(rownames(norm), 'sample')
        return(pb)
    }
}



de <- function(feature, df) {
    model_df <- df %>% select(feature, Unified_Visit, Age, Sex, Type, avg_count, avg_mt) %>% rename(Exp = feature)
    if (sum(model_df$Exp > 0) > 0.10 * nrow(model_df)) { 
        m_0 <- lm(Exp ~ log(avg_count) + avg_mt + scale(Age) + Sex, data = model_df)
        m_1 <- lm(Exp ~ log(avg_count) + avg_mt + scale(Age) + Sex + Type, data = model_df)
        ANNO <- anova(m_0, m_1)
        LRP <- ANNO[2,6]
        F <- ANNO[2,5]
        Beta <- summary(m_1)$coefficients['Type', 'Estimate']
        SE <- summary(m_1)$coefficients['Type', 'Std. Error']
        res <- c(gene = feature, LRP = LRP, F = F, Beta = Beta, SE = SE)
    }
    else {
        res <- c(gene = feature, LRP = NA, F = NA, Beta = NA, SE = NA)
    }
    return(res)
}