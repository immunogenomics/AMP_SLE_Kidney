
mergeMetas <- function(type) {
    meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/pbmc/', type, '/case_control/meta.csv'))
    clinData <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/pbmc/sex/meta.csv'))
}
