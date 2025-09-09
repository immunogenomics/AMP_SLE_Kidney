#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")

#-------------------------------------------------------

samples <- read.delim("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/list_of_cell_samples.txt",
                            sep = "\n", header = F)
samples <- as.character(samples[,1])

samples <- samples[!(samples %in% c("AMPSLEkid_cells_1258", "AMPSLEkid_cells_0135"))]

#-------------------------------------------------------

## FILTERING PARAMETERS
### NGENE : 500
### NUMI : 1000
### MT : 3% non-targeted

nGenepar <- 500
nUMIpar <- 1000
pctmt <- 3

#-------------------------------------------------------


meta_full <- c()
#data_full <- c()


for (i in 1:length(samples)) {

## Read in filtered data
    data <- readRDS(paste("../cell-seq-output/raw_data/corrected_counts_",
                                 samples[i],
				 ".rds",
                      		 sep = "")
                   )

   colnames(data) <- paste(colnames(data), samples[i], sep = "-")
   print(samples[i])
   sample <- samples[i]
## Make Raw Meta
    meta_raw <- data.frame(cells = colnames(data),
                       nGene = colSums(data > 0),
                       nUMI = colSums(data),
                       sample = sample
                       )
## Identify mitochondrial genes
    mito_genes <- rownames(data)[grep("^Mt-", rownames(data),
                                                        value = FALSE, ignore.case = TRUE)]
## Calculate % Mitochondrial UMI
    meta_raw$percent.mito <- Matrix::colSums(data[mito_genes, ])/Matrix::colSums(data)

## Calculate $ Mitochondrial excluding DASHed MT-genes
    mito_genes_include <- c("MT-ND5", "MT-ND6")
    mito_genes_exclude <- c('MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP6','MT-ATP8','MT-CO3','MT-ND3','MT-ND4','MT-ND4L','MT-CYB')
    tmp <- data[!(rownames(data) %in% mito_genes_exclude),]
    meta_raw$percent.mito.sub <- Matrix::colSums(tmp[mito_genes_include, ])/Matrix::colSums(tmp)
    tmp <- data[!(rownames(data) %in% mito_genes),]
    meta_raw$nUMInonMT <- colSums(tmp)

    rm(tmp)
    rm(data)
    gc()

## Read in metadata
#    meta <- readRDS(paste("../cell-seq-output/filtered_metadata/",
#                                 samples[i],
#                                 "_meta_filtered_",
#                                 nGenepar,
#                                 "nGene_",
#                                 nUMIpar,
#                                 "nUMI_",
#                                 pctmt,
#                                 "pctnontargetMT.rds",
#                                 sep = "")
#                   )

## Add to full meta and data
    
#    data_full <- cbind(data_full, data)
    meta_full <- rbind(meta_full, meta_raw)    
}

print("Dimensions of full meta and data...")
print(dim(meta_full))
#print(dim(data_full))

## Merge meta with map
#cell_map <- readRDS("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/outside_meta/cell_iso_map.rds")
#meta_full <- merge(meta_full, cell_map, by = "sample", all.x = T)
#print("Have we added columns?")
#print(dim(meta_full))


#saveRDS(meta_full, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
#                                 "meta_filtered_all",
#                                 nGenepar,
#                                 "nGene_",
#                                 nUMIpar,
#                                 "nUMI_",
#                                 pctmt,
#                                 "pctnontargetMT.rds",
#                                 sep = "")
#       )

saveRDS(meta_full, paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/full_metadata/",
                           "meta_corrected_all.rds",
                            sep = "")
       )

rm(meta_full)
#rm(data_full)
#rm(meta)
rm(data)
rm(samples)
rm(nGenepar)

gc()


