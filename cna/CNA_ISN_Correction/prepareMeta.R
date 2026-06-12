suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(tidyverse)
    library(caret)
})

set.seed(0)

args <- commandArgs(trailingOnly=TRUE)

tissue_meta <- readRDS(as.character(args[1]))  # Tissue meta
pbmc_meta <- readRDS(as.character(args[2]))  # PBMC meta
pbmc_hpcs <- readRDS(as.character(args[3]))  # PBMC hPCs

name <- as.character(args[4])

# Clinical information, and whether it comes from the first or second biopsy
clinical <- readRDS('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/clinical_data_05042023.rds')
first_biop_pred <- readRDS("/data/srlab2/qxiao/AMP-SLE/data/clinical/df_pred_biop.rds")

# Add clinical data and biopsy information
clinical <- merge(clinical, first_biop_pred, by.x = "AMP.Subject_ID", by.y = "Subject_ID")

# Overall, clean meta file
sampleToDonor <- data.frame(
    sample = clinical$sample, Unified_Visit = clinical$AMP.ID,
    ISN = clinical$Final_ISN,
    Activity = clinical$Final_Activity,
    Chronicity =  clinical$Final_Chronicity,
    Responder = clinical$Responder.Status, Response12 = clinical$Clinical.Response.at.12.weeks,
    Response26 = clinical$Clinical.Response.at.26.weeks, Response52 = clinical$Clinical.Response.at.52.weeks,
    First_biop = clinical$First_biop
)

# 
tissue_meta <- merge(tissue_meta, sampleToDonor, by = "sample")
pbmc_meta <- merge(pbmc_meta, sampleToDonor, by = "Unified_Visit")

tissue_cna_meta <- data.frame(
    barcode = tissue_meta$cell,
    sample = tissue_meta$sample,
    cellType = tissue_meta$final_annotation,
    batch = tissue_meta$processing.batch,
    unified_visit = tissue_meta$Unified_Visit,
    chronicity = tissue_meta$Chronicity,
    ISN = tissue_meta$ISN,
    activity = tissue_meta$Activity,
    responder = tissue_meta$Responder, 
    responder12 = tissue_meta$Response12,
    responder26 = tissue_meta$Response26, 
    responder52 = tissue_meta$Response52,
    first_biopsy = tissue_meta$First_biop
)

tissue_cna_meta <- tissue_cna_meta[!(is.na(tissue_cna_meta$activity)), ]

# One Hot Encoder
oneHotify <- function(meta, variable) {
    dummy <- dummyVars(" ~ .", data = meta[variable])
    oneHotified <- predict(dummy, newdata = meta[variable])
    return(cbind(meta, oneHotified))
}

tissue_cna_meta <- oneHotify(tissue_cna_meta, "ISN")
tissue_cna_meta <- oneHotify(tissue_cna_meta, "responder")
tissue_cna_meta <- oneHotify(tissue_cna_meta, "responder12")
tissue_cna_meta <- oneHotify(tissue_cna_meta, "responder26")
tissue_cna_meta <- oneHotify(tissue_cna_meta, "responder52")

pbmc_cna_meta <- data.frame(
    barcode = pbmc_meta$Cell,
    cellType = pbmc_meta$annotation,
    batch = pbmc_meta$Batch,
    unified_visit = pbmc_meta$Unified_Visit,
    chronicity = pbmc_meta$Chronicity,
    ISN = pbmc_meta$ISN,
    activity = pbmc_meta$Activity,
    responder = pbmc_meta$Responder, 
    responder12 = pbmc_meta$Response12,
    responder26 = pbmc_meta$Response26, 
    responder52 = pbmc_meta$Response52, 
    first_biopsy = pbmc_meta$First_biop
)

pbmc_cna_meta <- pbmc_cna_meta[!(is.na(pbmc_cna_meta$activity)), ]

# One Hot Encoder
oneHotify <- function(meta, variable) {
    dummy <- dummyVars(" ~ .", data = meta[variable])
    oneHotified <- predict(dummy, newdata = meta[variable])
    return(cbind(meta, oneHotified))
}

pbmc_cna_meta <- oneHotify(pbmc_cna_meta, "ISN")
pbmc_cna_meta <- oneHotify(pbmc_cna_meta, "responder")
pbmc_cna_meta <- oneHotify(pbmc_cna_meta, "responder12")
pbmc_cna_meta <- oneHotify(pbmc_cna_meta, "responder26")
pbmc_cna_meta <- oneHotify(pbmc_cna_meta, "responder52")

pbmc_hpcs <- pbmc_hpcs[rownames(pbmc_hpcs) %in% pbmc_cna_meta$barcode, c(1:20)]
colnames(pbmc_hpcs) <- paste0("hPC", seq(1,20))

pbmc_meta <- pbmc_meta[pbmc_meta$Cell %in% pbmc_cna_meta$barcode, ]
pbmc_UMAP <- data.frame(UMAP_1 = pbmc_meta$UMAP_1, UMAP_2 = pbmc_meta$UMAP_2)

tissue_meta <- tissue_meta[tissue_meta$cell %in% tissue_cna_meta$barcode, ]
tissue_hpcs <- tissue_meta[, paste0("hPC", seq(1,20))]
tissue_UMAP <- data.frame(UMAP_1 = tissue_meta$hUMAP1, UMAP_2 = tissue_meta$hUMAP2)


dir.create("metaFiles/")
write.csv(tissue_cna_meta, paste0("metaFiles/", name, "_tissue_cna_meta.csv"), row.names = FALSE)
write.csv(pbmc_cna_meta, paste0("metaFiles/", name, "_pbmc_cna_meta.csv"), row.names = FALSE)

write.csv(tissue_hpcs, paste0("metaFiles/", name, "_tissue_cna_hpcs.csv"), row.names = FALSE)
write.csv(pbmc_hpcs, paste0("metaFiles/", name, "_pbmc_cna_hpcs.csv"), row.names = FALSE)

write.csv(tissue_UMAP, paste0("metaFiles/", name, "_tissue_cna_umap.csv"), row.names = FALSE)
write.csv(pbmc_UMAP, paste0("metaFiles/", name, "_pbmc_cna_umap.csv"), row.names = FALSE)