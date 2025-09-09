suppressPackageStartupMessages({
    library(tidyverse)
    library(dplyr)
})

args = commandArgs(trailingOnly=TRUE)

cellType <- as.character(args[1])

ptM <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', cellType, '/injured_pt/sc_meta.csv'))
ptH <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', cellType, '/injured_pt/sc_harmony.csv'))
ptU <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', cellType, '/injured_pt/sc_umap.csv'))

line <- paste0(cellType, ",", 
    nrow(ptM), ",", nrow(ptM[!is.na(ptM$Final_Chronicity), ]))

write(line,file="cellCounts.csv",append=TRUE)

dir.create(paste0("ptChronicities/", cellType))

write.csv(ptM[!is.na(ptM$Final_Chronicity), ], paste0("ptChronicities/", cellType, "/sc_meta.csv"), row.names = FALSE)
write.csv(ptH[!is.na(ptM$Final_Chronicity), ], paste0("ptChronicities/", cellType, "/sc_harmony.csv"), row.names = FALSE)
write.csv(ptU[!is.na(ptM$Final_Chronicity), ], paste0("ptChronicities/", cellType, "/sc_umap.csv"), row.names = FALSE)