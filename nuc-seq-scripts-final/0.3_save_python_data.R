#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")


library(RcppCNPy)

#-------------------------------------------------------

filepath1 <- "/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/"
filepath2 <- "/outs/filtered_feature_bc_matrix/"
samples <- read.delim("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/list_of_nuc_samples.txt",
                            sep = "\n", header = F)
samples <- as.character(samples[,1])

#-------------------------------------------------------

for (i in 1:length(samples)) {

	data <- readRDS(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/filtered_data/",
                      samples[i],
                      "_data_filtered_500nGene_0.01pctmtUMI7500uppernGene_1000nUMI_40000uppernUMI.rds",
                      sep = "")
)
	dim(data)
	npySave(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/python_data/", 
		samples[i], 
		"_data_filtered_500nGene_0.01pctmtUMI7500uppernGene_1000nUMI_40000uppernUMI.py",
		sep = ""),
		as.matrix(data))
	rm(data)
	gc()
}
