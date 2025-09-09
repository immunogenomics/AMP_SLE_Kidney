suppressPackageStartupMessages({
    library(tidyverse)
})

cleanUp <- function(name, type, tName) {
    pMeta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/pbmc/', name, '/', type, '/meta.csv'))
    pHarmony <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/pbmc/', name, '/', type, '/harmony.csv'))
    pUmap <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/pbmc/', name, '/', type, '/umap.csv'))
    tMeta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/', tName, '/', type, '/sc_meta.csv'))

    clinData <- readRDS('/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects/clinicalData.rds') %>% 
            select(AMP.ID, AMP.Subject_ID, sample)

    tMeta <- merge(tMeta, clinData, by = "sample", all.x = TRUE)

    cellsToKeep <- pMeta$Unified_Visit %in% unique(tMeta$AMP.ID)

    suppressWarnings({
        dir.create("cna_results")
        dir.create(paste0("cna_results/", name))
        dir.create(paste0("cna_results/", name, "/", type))
    })

    write.csv(pMeta[cellsToKeep, ], paste0("cna_results/", name, "/", type, "/meta.csv"), row.names = FALSE)
    write.csv(pHarmony[cellsToKeep, ], paste0("cna_results/", name, "/", type, "/harmony.csv"), row.names = FALSE)
    write.csv(pUmap[cellsToKeep, ], paste0("cna_results/", name, "/", type, "/umap.csv"), row.names = FALSE)

    oldSamples <- length(unique(pMeta$Unified_Visit))
    newSamples <- length(unique(pMeta$Unified_Visit[cellsToKeep]))
    oldCells <- nrow(pMeta)
    newCells <- nrow(pMeta[cellsToKeep, ])

    newRow <- data.frame(
        type = type, 
        samplesBefore = oldSamples, 
        samplesAfter = newSamples, 
        cellsBefore = oldCells,
        cellsAfter = newCells
    )

    return(newRow)
}

for (cellTypes in list(list("mono_dc", "myeloid"), list("t_nk", "t_nk"), list("b_cell", "b_plasma"))) {
    dir.create("sampleCounts/")
    powerCounts <- data.frame(
        type = character(),  # Empty character column
        samplesBefore = numeric(),   # Empty numeric column
        samplesAfter = numeric(),   # Empty numeric column
        cellsBefore = numeric(),   # Empty numeric column
        cellsAfter = numeric()   # Empty numeric column
    )

    for (type in c("case_control", "chronicity", "activity")) {
        if (cellTypes[[1]] == "b_cell" && type == "case_control") {next}
        powerCounts <<- rbind(powerCounts, cleanUp(cellTypes[[1]], type, cellTypes[[2]]))
    }

    write.csv(powerCounts, paste0("sampleCounts/", cellTypes[[1]], "_powerCounts.csv"), row.names = FALSE)
}