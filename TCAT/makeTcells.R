suppressPackageStartupMessages({
    library(tidyverse)
    library(Matrix)
})

t_raw <- readRDS("/data/srlab1/Yu/tissue_blood_eQTL/AMPII_SLE/AMPII_blood/data/pheno/AllSample/GEX_raw.rds")
finalCells <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects_v3/t_nk_pbmc_norm.Rds")

normalize_channel <- function(x) {
  x %>%
    gsub("channel_", "channel", .) %>%
    gsub("_", "-", .) %>%
    gsub("channel([0-9]+)-", "channel\\1-", .)
}

colnames(finalCells) <- normalize_channel(colnames(finalCells))
tcells <- t_raw[, colnames(t_raw) %in% colnames(finalCells)]

writeMM(tcells, file = "onlyT_cells.mtx")

write.table(rownames(tcells),
            file = "onlyT_rows.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(colnames(tcells),
            file = "onlyT_cols.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)