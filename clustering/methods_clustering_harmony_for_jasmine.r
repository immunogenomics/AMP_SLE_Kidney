#@Script: Normalize, ScaleData, Harmonize, UMAP, Cluster For Jasmine
#@Author: Sid Gurajala
#@Date: 11/22/2023

library(dplyr)
library(singlecellmethods)
library(harmony)
library(Seurat)

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
                           
vargenes_vst <- function(object, groups, topn, loess.span = 0.3) {
    clip.max <- sqrt(ncol(object))

    N <- ncol(object)
    if (missing(groups)) {
        groups <- rep('A', N)
    }
    
    res <- split(seq_len(N), groups) %>% lapply(function(idx) {
        object_group <- object[, idx]
        ## row means
        hvf.info <- data.frame(
          symbol = rownames(object_group), 
          mean = Matrix::rowMeans(object_group)
        )

        ## row vars
        hvf.info$variance <- rowVars(object_group, hvf.info$mean)

        ## initialize
        hvf.info$variance.expected <- 0
        hvf.info$variance.standardized <- 0

        not.const <- hvf.info$variance > 0

        ## loess curve fit 
        suppressWarnings({
            fit <- loess(formula = log10(variance) ~ log10(mean), 
                data = hvf.info[not.const, ], span = loess.span)            
        })

        ## extract fitted variance 
        hvf.info$variance.expected[not.const] <- 10^fit$fitted

        ## get row standard deviations after clipping
        hvf.info$variance.standardized <- rowVarsStd(
            object_group, 
            hvf.info$mean, 
            sqrt(hvf.info$variance.expected), 
            clip.max
        )

        hvf.info <- hvf.info %>% 
#             tibble::rownames_to_column('symbol') %>% 
            arrange(-variance.standardized) %>% 
            tibble::rowid_to_column('rank') %>% 
            transform(group = unique(groups[idx]))

        return(hvf.info)        
    })
    
    
    if (missing(topn)) {
        ## MODE 1: return table 
        res <- Reduce(rbind, res) %>% 
            dplyr::select(group, symbol, rank, everything())

        if (length(unique(res$group)) == 1) {
            res$group <- NULL
        }
    } else {
        ## MODE 2: return genes
        res <- lapply(res, function(x) head(x, topn)$symbol)
    }
    return(res)
}



#list of excluded genes from vargenes
genes_exclude <- c(grep("^MT-|^RPL|^RPS|MALAT1|MIR-", row.names(glom_norm), value = TRUE), Seurat::cc.genes$s.genes, Seurat::cc.genes$g2m.genes) 

#split metadata into single cell and single nuc
sc_meta <- meta %>% filter(dataset == "scRNAseq")
sn_meta <- meta %>% filter(dataset == "snRNAseq")

#split normalized data
sc_norm <- norm[, sc_meta$cell]
sn_norm <- norm[, sn_meta$cell]

#restrict to genes with measurable expression
sc_threshold_genes <- sc_norm[rownames(sc_norm[rowMeans(sc_norm) >= .01, ]), ]
sn_threshold_genes <- sn_norm[rownames(sn_norm[rowMeans(sn_norm) >= .01, ]), ]

#remove excluded genes
sc_threshold_genes <- sc_threshold_genes[ -which( rownames(sc_threshold_genes) %in% genes_exclude ) ,]
sn_threshold_genes <- sn_threshold_genes[ -which( rownames(sn_threshold_genes) %in% genes_exclude ) ,]


#restrict to samples with > 25 cells for var gene selection
sc_samples_25 <- sc_meta %>% group_by(sample) %>% tally() %>% filter(n > 25) %>% pull(sample)

sc_threshold_samples <- sc_threshold_genes[, sc_meta %>% filter(sample %in% sc_samples_25) %>% pull(cell)]


sn_samples_25 <- sn_meta %>% group_by(sample) %>% tally() %>% filter(n > 25) %>% pull(sample)

sn_threshold_samples <- sn_threshold_genes[, sn_meta %>% filter(sample %in% sn_samples_25) %>% pull(cell)]


#find top 500 variable genes per sample per modality
sc_var_genes_500_raw <- vargenes_vst(object = sc_threshold_samples, groups = sc_meta %>% filter(sample %in% sc_samples_25) %>% pull(sample), topn = 500)

sn_var_genes_500_raw <- vargenes_vst(object = sn_threshold_samples, groups = sn_meta %>% filter(sample %in% sn_samples_25) %>% pull(sample), topn = 500)
 
#find intersection of sc and sn var features
inter_var <- intersect(sc_var_genes_500_raw, sn_var_genes_500_raw)

#scale data with Seurat function
scaled <- norm[ vargenes_df,] %>% ScaleDataSeurat() 

#find sc/sn weights for PCA                    
table(meta$dataset)
y <- factor(meta$dataset)
weights <- as.numeric((1 / prop.table(table(y)))[y]) / nlevels(y)    
                    
#reformat as sparse matrix
temp <- Matrix::Matrix(scaled, sparse = TRUE)
#perform weighted pca                      
pca_res <- weighted_pca(temp, weights, rownames(temp), 20, TRUE, 10)$embeddings
                      

#Harmony from dev version on github, change function to change full object                      
harmony <- HarmonyMatrix(pca_res[,1:20], 
                         meta, 
                         vars_use = c("sample", "processing.batch", "Site", 'dataset'), 
                         theta = c(1,0,0,1),
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

#reassign harmony colnames
colnames(harmony) <-  c("hPC-1", "hPC-2", "hPC-3", "hPC-4",
                      "hPC-5", "hPC-6", "hPC-7", "hPC-8",
                      "hPC-9", "hPC-10", "hPC-11", "hPC-12",
                      "hPC-13", "hPC-14", "hPC-15", "hPC-16",
                      "hPC-17", "hPC-18", "hPC-19", "hPC-20")

#make umap, save nn fgraph and model                      
humap_object <- uwot::umap(
    X = harmony[,1:20],
    ret_extra = c('nn', 'fgraph', 'model')
) 

#Run modularity clustering                      
resolution_list <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0)
ids_ref_cca <- Reduce(cbind, mclapply(resolution_list, function(res_use) {
    Seurat:::RunModularityClustering(SNN = humap_object$fgraph, modularity = 1, 
        resolution = res_use, algorithm = 1, n.start = 20, 
        n.iter = 20, random.seed = 100, print.output = FALSE, 
        temp.file.location = NULL, edge.file.name = NULL)
}, mc.preschedule = FALSE, mc.cores = min(20, length(resolution_list))))
gc()
colnames(ids_ref_cca) <- sprintf("hres.%.2f", resolution_list)                      