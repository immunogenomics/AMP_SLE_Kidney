suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
})

args = commandArgs(trailingOnly=TRUE)

name <- args[1]
type <- args[2]

sc_meta <- read.csv(paste0('/data/srlab/ssg34/SLE_kidney_v2/data/pbmc/', name, '/', type, '/meta.csv'))

meta <- read_excel("bloodMeta.xlsx")
meta <- meta %>% select("Subject ID", "Sex", "Age at Baseline") 
colnames(meta) <- c("Sample", "Sex", "Age")
meta <- meta %>%
  filter(!is.na(Sex), !is.na(Age)) %>%
  mutate(Sex_numeric = ifelse(Sex == "Male", 1, 0))

sc_meta <- sc_meta %>% left_join(meta,
            by = c("Donor" = "Sample"))

metaDir <- paste0("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/CNA_blood/new_metas/", name, "/", type)
dir.create(metaDir)
write.csv(sc_meta, paste0(metaDir, "/meta.csv"))