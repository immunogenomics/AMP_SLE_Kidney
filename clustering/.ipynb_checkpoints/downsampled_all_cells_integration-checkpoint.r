#@Author: Sid Gurajala
#@Date: 11/3/2023
#@Script: Integrate 100K Downsampled Cells

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
library(lisi)
source("/data/srlab/anathan/scripts/scseq_utils.R")
library(parallel)
library('lme4', lib.loc = "/PHShome/ssg34/.conda/envs/jupy/lib/R/library")
library(glmmTMB, lib.loc ="/PHShome/ssg34/.conda/envs/jupy/lib/R/library")
.libPaths("/PHShome/ssg34/.conda/envs/plswork/lib/R/library")


## Read in data

BuildSNNSeurat <- function (data.use, k.param = 30, prune.SNN = 1/15, nn.eps = 0) {
    my.knn <- nn2(data = data.use, k = k.param, searchtype = "standard", eps = nn.eps)
    nn.ranked <- my.knn$nn.idx

    snn_res <- ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    rownames(snn_res) <- row.names(data.use)
    colnames(snn_res) <- row.names(data.use)
    return(snn_res)
}
#environment(BuildSNNSeurat) <- asNamespace("Seurat")

NormalizeDataSeurat <- function(A, scaling_factor = 1e4, do_ftt = FALSE) {
	A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
    A@x <- scaling_factor * A@x
    if (do_ftt) {
        A@x <- sqrt(A@x) + sqrt(1 + A@x)
    } else {
        A@x <- log(1 + A@x)
    }
	return(A)
}

ScaleDataSeurat <- function (data.use, margin = 1, scale.max = 10,
                                block.size = 1000) {

    if (margin == 2) data.use %<>% t
    max.block <- ceiling(nrow(data.use)/block.size)

    ## Define data and functions to use in sparse and dense cases
    if (class(data.use) == "dgCMatrix" | class(data.use) == "dgTMatrix") {
        scale_fxn <- function(x) {
            FastSparseRowScale(mat = x, scale = TRUE, center = TRUE,
                               scale_max = scale.max, display_progress = FALSE)
        }
    } else {
        scale_fxn <- function(x) {
            FastRowScale(mat = x, scale = TRUE, center = TRUE,
                               scale_max = scale.max, display_progress = FALSE)
       }
        data.use <- as.matrix(data.use)
    }

    ## Do scaling, at once or in chunks
    if (max.block == 1) {
        scaled.data <- scale_fxn(data.use)
    } else {
        scaled.data <- matrix(NA, nrow(data.use), ncol(data.use))
        for (i in 1:max.block) {
            idx.min <- (block.size * (i - 1))
            idx.max <- min(nrow(data.use), (block.size * i - 1) + 1)
            my.inds <- idx.min:idx.max
            scaled.data[my.inds, ] <- scale_fxn(data.use[my.inds, , drop = F])
        }
    }

    colnames(scaled.data) <- colnames(data.use)
    row.names(scaled.data) <- row.names(data.use)
    scaled.data[is.na(scaled.data)] <- 0
    if (margin == 2) scaled.data %<>% t
    return(scaled.data)
}
#environment(ScaleDataSeurat) <- asNamespace("Seurat")

FindVariableGenesBatch <- function(exprs_mat, meta_df, genes_exclude = NULL, ngenes_use = 1e3, expr_min = .1) {
    if (!is.null(genes_exclude)) {
        genes_use <- setdiff(row.names(exprs_mat), genes_exclude)
    }
    x_res <- split(meta_df$cell, meta_df$sample) %>% lapply(function(x) {
        FindVariableGenesSeurat(exprs_mat[genes_use, x]) %>% 
            subset(gene.mean >= expr_min) %>% 
            tibble::rownames_to_column("gene") %>% 
            dplyr::arrange(-gene.dispersion) %>%
            head(ngenes_use)
    })
    data.table(Reduce(rbind, x_res))[, .N, by = gene][order(-N)]    
}


FindVariableGenesSeurat <- function (data, x.low.cutoff = 0.1, x.high.cutoff = 8,
                                     y.cutoff = 1, y.high.cutoff = Inf, num.bin = 0,
                                     binning.method = "equal_width", sort.results = TRUE,
                                     display.progress = TRUE, ...)
{
    genes.use <- rownames(data)
    if (class(data) != "dgCMatrix") {
        data <- as(as.matrix(data), "dgCMatrix")
    }
    ## (1) get means and variances
    gene.mean <- FastExpMean(data, display.progress)
    names(gene.mean) <- genes.use
    gene.dispersion <- FastLogVMR(data, display.progress)
    names(gene.dispersion) <- genes.use

    gene.dispersion[is.na(x = gene.dispersion)] <- 0
    gene.mean[is.na(x = gene.mean)] <- 0

    mv.df <- data.frame(gene.mean, gene.dispersion)
    rownames(mv.df) <- rownames(data)

    ## (OPTIONAL) do the binning correction
    if (num.bin > 0) {
      if (binning.method == "equal_width") {
          data_x_bin <- cut(x = gene.mean, breaks = num.bin)
      }
      else if (binning.method == "equal_frequency") {
          data_x_bin <- cut(x = gene.mean, breaks = c(-1, quantile(gene.mean[gene.mean >
              0], probs = seq(0, 1, length.out = num.bin))))
     }
      else {
          stop(paste0("Invalid selection: '", binning.method,
              "' for 'binning.method'."))
      }
      names(x = data_x_bin) <- names(x = gene.mean)
      mean_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
          FUN = mean)
      sd_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
          FUN = sd)
      gene.dispersion.scaled <- (gene.dispersion - mean_y[as.numeric(x = data_x_bin)])/sd_y[as.numeric(x = data_x_bin)]
      gene.dispersion.scaled[is.na(x = gene.dispersion.scaled)] <- 0
      ##names(gene.dispersion.scaled) <- names(gene.mean)

      mv.df$gene.dispersion.scaled <- gene.dispersion.scaled
    }

    return(mv.df)
}

FastSparseRowScale <- function(mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = TRUE) {
    .Call('_Seurat_FastSparseRowScale', PACKAGE = 'Seurat', mat, scale, center, scale_max, display_progress)
}

# weighted PCA is in immunogenomics/singlecellmethods, function weighted_pca
#library(singlecellmethods)
weighted_pca <- function(X, weights, genes_use=NULL, npc=20, do_corr=TRUE, scale_thresh=10) {
    if (!identical(length(weights), ncol(X))) {
        stop('Columns in X must match length of weights')
    }
    
#     y <- factor(y)
#     weights <- as.numeric((1 / prop.table(table(y)))[y]) / nlevels(y)
    if (any(is.na(weights))) {
        idx_keep <- which(is.na(weights))
#         y <- y[idx_keep]
        weights <- weights[idx_keep]
        X <- X[, idx_keep]
    }
    if (is.null(genes_use)) {
        genes_use <- row.names(X)
    } else if (length(genes_use) < nrow(X)) {
        if (any(!genes_use %in% row.names(X))) {
            stop('genes_use not in rownames of X')
        }
        X <- X[genes_use, ]
    }
    
    ## weighted z-scores
#     mu <- X %>% apply(1, function(x) {SDMTools:::wt.mean(x, weights)})
#     sig <- X %>% apply(1, function(x) {SDMTools:::wt.sd(x, weights)})
    mu <- singlecellmethods::rowMeans(X, weights = weights)
    sig <- singlecellmethods::rowSDs(X, weights = weights)
    
    # Added 12/9/19: save weighted scaling means and std devs
    vargenes_means_sds <- tibble(
        symbol = genes_use,
        mean = mu
    )
    vargenes_means_sds$stddev <- sig
    # finish added 12/9/19
    
    X <- scaleDataWithStats(X, mu, sig) 
    X <- X[which(is.na(rowSums(X)) == 0), ]
    if (do_corr) {
        X <- X %>% scale() %>% pmin(scale_thresh) %>% pmax(-scale_thresh)
    }
    
    ## weighted SVD
#     pres <- rsvd::rsvd(X %*% Matrix::Diagonal(x = sqrt(weights)), k = npc)
    pres <- RSpectra::svds(X %*% Matrix::Diagonal(x = sqrt(weights)), npc)
    V <- (Matrix::Diagonal(x = 1 / sqrt(weights)) %*% pres$v) %*% diag(pres$d)
    V <- as.matrix(V)
    colnames(V) <- paste0('PC', 1:npc)
    row.names(V) <- colnames(X)
    colnames(pres$u) <- paste0('PC', 1:npc)
    row.names(pres$u) <- row.names(X)
    return(list(loadings = pres$u, embeddings = V, vargenes = vargenes_means_sds))
}

# Cosine normalize values
cosine_normalize <- function(X, MARGIN = 1, do_safe = TRUE) {
    if (do_safe) {
        X <- sweep(X, MARGIN, apply(X, MARGIN, max), "/")
    }
    sweep(X, MARGIN, apply(X, MARGIN, function(x) sqrt(sum(x^2))), "/")
}
                           
                                                    
do_pca <- function(X, weights, genes_use=NULL, npc=10, do_corr=TRUE) {
    if (is.null(genes_use)) {
        genes_use <- row.names(X)
    }
    mu <- X[genes_use, ] %>% apply(1, function(x) {SDMTools:::wt.mean(x, weights)})
    sig <- X[genes_use, ] %>% apply(1, function(x) {SDMTools:::wt.sd(x, weights)})
    
    X <- X[genes_use, ] %>% scaleDataWithStats(mu, sig) 
    X <- X[which(is.na(rowSums(X)) == 0), ]
    if (do_corr) {
        X <- scale(X)
    }
    pres <- rsvd::rsvd(X, k = npc)
    V <- pres$v %*% diag(pres$d)
    V <- data.table(V)
    colnames(V) <- paste0('PC', 1:npc)
    return(V)    
}
                           
get_stats <- function(X, weights, genes_use=NULL) {
    if (is.null(genes_use)) {
        genes_use <- row.names(X)
    }
    mu <- X[genes_use, ] %>% apply(1, function(x) {SDMTools:::wt.mean(x, weights)})
    sig <- X[genes_use, ] %>% apply(1, function(x) {SDMTools:::wt.sd(x, weights)})
    return(list(mu = mu, sig = sig))
}
                           

meta <- readRDS(
        "/data/srlab/ssg34/SLE_kidney_v2/data/downsampled_all/downsampled_meta.rds")
norm <- readRDS(
        "/data/srlab/ssg34/SLE_kidney_v2/data/downsampled_all/downsampled_norm.rds")

genes_exclude <- c(grep("^MT-|^RPL|^RPS|MALAT1|MIR-", row.names(norm), value = TRUE), Seurat::cc.genes$s.genes, Seurat::cc.genes$g2m.genes) 

sc_meta <- meta %>% filter(dataset == "scRNAseq")
sn_meta <- meta %>% filter(dataset == "snRNAseq")

sc_norm <- norm[, sc_meta$cell]
sn_norm <- norm[, sn_meta$cell]

sc_threshold_genes <- sc_norm[rownames(sc_norm[rowMeans(sc_norm) >= .01, ]), ]
sn_threshold_genes <- sn_norm[rownames(sn_norm[rowMeans(sn_norm) >= .01, ]), ]

sc_threshold_genes <- sc_threshold_genes[ -which( rownames(sc_threshold_genes) %in% genes_exclude ) ,]
sn_threshold_genes <- sn_threshold_genes[ -which( rownames(sn_threshold_genes) %in% genes_exclude ) ,]

sc_samples_25 <- sc_meta %>% group_by(sample) %>% tally() %>% filter(n > 25) %>% pull(sample)

sc_threshold_samples <- sc_threshold_genes[, sc_meta %>% filter(sample %in% sc_samples_25) %>% pull(cell)]


sn_samples_25 <- sn_meta %>% group_by(sample) %>% tally() %>% filter(n > 25) %>% pull(sample)

sn_threshold_samples <- sn_threshold_genes[, sn_meta %>% filter(sample %in% sn_samples_25) %>% pull(cell)]

sc_var_genes_1000_raw <- vargenes_vst(object = sc_threshold_samples, groups = sc_meta %>% filter(sample %in% sc_samples_25) %>% pull(sample), topn = 500)

sn_var_genes_1000_raw <- vargenes_vst(object = sn_threshold_samples, groups = sn_meta %>% filter(sample %in% sn_samples_25) %>% pull(sample), topn = 500)


saveRDS(sc_var_genes_1000_raw, '/data/srlab/ssg34/SLE_kidney_v2/data/downsampled_all/vargenes_sc_raw.rds')
saveRDS(sc_var_genes_1000_raw, '/data/srlab/ssg34/SLE_kidney_v2/data/downsampled_all/vargenes_sn_raw.rds')

inter_var <- intersect(sc_var_genes_1000_raw, sn_var_genes_1000_raw) 

vargenes_df <- inter_var

norm_scaled <- norm[ vargenes_df,] %>% ScaleDataSeurat() 
                           
table(meta$dataset)
y <- factor(meta$dataset)
weights <- as.numeric((1 / prop.table(table(y)))[y]) / nlevels(y)

temp <- Matrix::Matrix(norm_scaled, sparse = TRUE)

pca_res <- weighted_pca(temp, weights, rownames(temp), 20, TRUE, 10)#$embeddings
                           
saveRDS(pca_res, '/data/srlab/ssg34/SLE_kidney_v2/data/downsampled_all/pca.rds')

pca_res <- pca_res$embeddings 

umap_object <- uwot::umap(
    X = pca_res[,1:20],
    ret_extra = c('nn', 'fgraph', 'model')
) 
meta$uwotUMAP1 <- umap_object$embedding[, 1]
meta$uwotUMAP2 <- umap_object$embedding[, 2]                           

library(harmony)

harmony <- HarmonyMatrix(pca_res[,1:20], glom_meta, vars_use = c("sample", "processing.batch", "Site", 'dataset'), theta = c(1,0,0,1),
                         # lambda = 1, 
#                          tau = 0, 
                          epsilon.cluster = -Inf,
                          epsilon.harmony = -Inf,
                         sigma = 0.15,
                         max.iter.cluster = 50,
                         max.iter.harmony = 20,
                         plot_convergence = TRUE,
                             npcs = 20,
                            do_pca = F,
                            weights = weights)

colnames(harmony) <-  c("hPC-1", "hPC-2", "hPC-3", "hPC-4",
                      "hPC-5", "hPC-6", "hPC-7", "hPC-8",
                      "hPC-9", "hPC-10", "hPC-11", "hPC-12",
                      "hPC-13", "hPC-14", "hPC-15", "hPC-16",
                      "hPC-17", "hPC-18", "hPC-19", "hPC-20")
                           
humap_object <- uwot::umap(
    X = harmony[,1:20],
    ret_extra = c('nn', 'fgraph', 'model')
) 
                           
meta$huwotUMAP1 <- humap_object$embedding[, 1]
meta$huwotUMAP2 <- humap_object$embedding[, 2]
meta <- cbind(meta, harmony[, 1:20])
                           
saveRDS(harmony, '/data/srlab/ssg34/SLE_kidney_v2/data/downsampled_all/harmony.rds')  
saveRDS(meta, '/data/srlab/ssg34/SLE_kidney_v2/data/downsampled_all/final_meta.rds')                           