source("/data/srlab/ssg34/SLE_kidney_v2/scripts/libs/kidney_utils.r")
set.seed(0)

args <- commandArgs(trailingOnly=TRUE)

norm_file_path <- as.character(args[1])
meta_file_path <- as.character(args[2])
cell_type <- as.character(args[3])


norm <- readRDS(input_file_path)[['RNA']]@data
meta <- readRDS(meta_file_path)

norm <- norm[, meta$Cell]

gene_list <-  mclapply(rownames(norm), check_5p_exp, norm, mc.cores = 20)

genes <- data.frame(genes = do.call(rbind, gene_list))

output_file_path <- paste0('/data/srlab/ssg34/SLE_kidney_v2/data/pseudobulk_de/', cell_type, '_5pgenes_de.rds')

saveRDS(genes, output_file_path)
